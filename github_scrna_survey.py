#!/usr/bin/env python3
"""
github_scrna_survey.py  —  v1.0

GitHub REST API scout for single-cell reproducibility repositories
attached to top-journal single-cell RNA-seq publications.

Features
--------
* Multi-query search with deduplication across queries
* Language composition via /repos/{owner}/{repo}/languages
* Recursive git-tree scan for notebook / Rmd / workflow / env / Docker files
* Journal detection with abbreviation support (longer / more-specific names preferred
  so "Nature Communications" is never swallowed by "Nature")
* Conservative paper-year extraction (only when year appears in publication context)
* Study vs. method scoring (study_score / method_score / repo_type)
* Fork exclusion by default (--include-forks to override)
* Rate-limit aware with exponential back-off (up to 5 retries per request)
* Outputs both CSV and JSONL with confidence / score fields

Usage
-----
    export GITHUB_TOKEN=ghp_...
    python github_scrna_survey.py

    # Custom queries
    python github_scrna_survey.py \
        --queries "scrna-seq analysis Nature 2023" \
                  "single cell reproducibility 2024" \
        --output-prefix my_survey \
        --max-repos 500 \
        --min-stars 5

    # Include forks; skip slow tree scan
    python github_scrna_survey.py --include-forks --no-tree-scan

Options
-------
    --token TOKEN          GitHub PAT (or GITHUB_TOKEN env var)
    --queries QUERY [...]  One or more search query strings
    --output-prefix PREFIX File prefix for .csv/.jsonl outputs [scrna_survey]
    --max-repos N          Maximum repositories to collect per query [200]
    --include-forks        Include forked repos (excluded by default)
    --no-tree-scan         Skip recursive tree scan
    --min-stars N          Minimum star count filter [0]
    --verbose              Enable debug logging
"""

import argparse
import base64
import csv
import json
import logging
import os
import re
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

import requests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

API_BASE = "https://api.github.com"
USER_AGENT = "scrna-survey/1.0 (github.com/r2ascension/scrna_py_2025)"

# Earliest publication year we consider plausible for scRNA-seq papers
MIN_VALID_YEAR = 2010

DEFAULT_QUERIES: List[str] = [
    "single-cell RNA-seq reproducibility code",
    "scRNA-seq analysis Nature Cell Science Genome",
    "single cell transcriptomics clinical data analysis",
    "scrna snakemake nextflow reproducible",
    "single-cell atlas code paper",
]

# Journal patterns: list of (regex_pattern, canonical_name, is_abbreviation) tuples.
# IMPORTANT: longer / more-specific entries MUST come before shorter ones
# so that "Nature Communications" is never swallowed by bare "Nature".
# Abbreviations come after full names; single-word journals come last.
#
# is_abbreviation values:
#   False — full multi-word name           → confidence 0.9
#   True  — recognised abbreviation OR     → confidence 0.7
#           single-word name               → confidence 0.5 (also in _SINGLE_WORD_JOURNALS)
JOURNAL_PATTERNS: List[Tuple[str, str, bool]] = [
    # Multi-word full names (most specific first)
    (r"nature\s+communications",            "Nature Communications",         False),
    (r"nature\s+biotechnology",             "Nature Biotechnology",           False),
    (r"nature\s+methods",                   "Nature Methods",                 False),
    (r"nature\s+medicine",                  "Nature Medicine",                False),
    (r"nature\s+cell\s+biology",            "Nature Cell Biology",            False),
    (r"nature\s+genetics",                  "Nature Genetics",                False),
    (r"nature\s+aging",                     "Nature Aging",                   False),
    (r"nature\s+immunology",                "Nature Immunology",              False),
    (r"genome\s+biology",                   "Genome Biology",                 False),
    (r"cell\s+reports\s+methods",           "Cell Reports Methods",           False),
    (r"cell\s+reports",                     "Cell Reports",                   False),
    (r"cell\s+systems",                     "Cell Systems",                   False),
    (r"cell\s+stem\s+cell",                 "Cell Stem Cell",                 False),
    (r"developmental\s+cell",               "Developmental Cell",             False),
    (r"cancer\s+cell",                      "Cancer Cell",                    False),
    (r"molecular\s+cell",                   "Molecular Cell",                 False),
    (r"current\s+biology",                  "Current Biology",                False),
    (r"science\s+advances",                 "Science Advances",               False),
    (r"science\s+translational\s+medicine", "Science Translational Medicine", False),
    (r"nucleic\s+acids\s+research",         "Nucleic Acids Research",         False),
    (r"plos\s+computational\s+biology",     "PLoS Computational Biology",     False),
    (r"plos\s+biology",                     "PLoS Biology",                   False),
    (r"briefings\s+in\s+bioinformatics",    "Briefings in Bioinformatics",    False),
    (r"frontiers\s+in\s+genetics",          "Frontiers in Genetics",          False),
    (r"bmc\s+genomics",                     "BMC Genomics",                   False),
    (r"bmc\s+bioinformatics",               "BMC Bioinformatics",             False),
    # Recognised abbreviations
    (r"nat\.?\s*commun",                    "Nature Communications",          True),
    (r"nat\.?\s*biotechnol",                "Nature Biotechnology",           True),
    (r"nat\.?\s*methods",                   "Nature Methods",                 True),
    (r"nat\.?\s*med\b",                     "Nature Medicine",                True),
    (r"nat\.?\s*cell\s+biol",               "Nature Cell Biology",            True),
    (r"nat\.?\s*genet",                     "Nature Genetics",                True),
    (r"nat\.?\s*immunol",                   "Nature Immunology",              True),
    (r"genome\s+biol\b",                    "Genome Biology",                 True),
    (r"nucleic\s+acids\s+res\b",            "Nucleic Acids Research",         True),
    (r"sci\.\s*adv\b",                      "Science Advances",               True),
    # Single-word / short names (least specific — must be last)
    (r"\bimmunity\b",                       "Immunity",                       True),
    (r"\belife\b",                          "eLife",                          True),
    (r"\bbioinformatics\b",                 "Bioinformatics",                 True),
    (r"\bnature\b",                         "Nature",                         True),
    (r"\bcell\b",                           "Cell",                           True),
    (r"\bscience\b",                        "Science",                        True),
]

# Pre-compiled patterns (compiled once at import time)
_JOURNAL_COMPILED: List[Tuple[re.Pattern, str, bool]] = [
    (re.compile(pat, re.IGNORECASE), name, is_abbrev)
    for pat, name, is_abbrev in JOURNAL_PATTERNS
]
_SINGLE_WORD_JOURNALS = {"Nature", "Cell", "Science", "Immunity", "eLife", "Bioinformatics"}


# Year pattern: 4-digit year adjacent to a publication/DOI context.
# We only return a year when it appears near a recognised anchor keyword
# so that arbitrary README dates (e.g. "Last updated: 2023-01-01") are ignored.
_YEAR_IN_CONTEXT_RE = re.compile(
    r"""
    (?:
        \b(?:published|year|doi|biorxiv|arxiv|preprint|citation|submitted|accepted)\b
        [^\d]{0,40}
        (20\d{2}|19\d{2})
    |
        (20\d{2}|19\d{2})
        [^\d]{0,40}
        \b(?:doi|biorxiv|preprint|published|
            nature|cell|science|genome\s+biol|
            nat\s+commun|nat\s+biotechnol|nat\s+methods|
            nucleic\s+acids)\b
    )
    """,
    re.IGNORECASE | re.VERBOSE,
)

# Study vs. method keyword sets
_STUDY_KEYWORDS: List[str] = [
    "scrna-seq", "scrna", "single-cell", "single cell",
    "snrna", "scatac", "multiome",
    "cell atlas", "cell type", "transcriptome", "transcriptomic",
    "tissue", "patient", "clinical", "disease", "cohort",
    "reproduce", "reproducibility", "data", "dataset",
    "analysis of", "profiling",
]

_METHOD_KEYWORDS: List[str] = [
    "tool", "package", "library", "software", "framework",
    "pipeline", "workflow", "method", "algorithm",
    "benchmarking", "benchmark", "utility",
]

# File classification helpers
_NOTEBOOK_EXTS = frozenset({".ipynb"})
_RMD_EXTS = frozenset({".rmd"})
_SNAKEMAKE_BASENAMES = frozenset({"snakefile"})
_SNAKEMAKE_EXTS = frozenset({".smk", ".snakefile"})
_NEXTFLOW_EXTS = frozenset({".nf"})
_NEXTFLOW_BASENAMES = frozenset({"nextflow.config"})
_DOCKER_BASENAMES = frozenset({
    "dockerfile", "docker-compose.yml", "docker-compose.yaml",
})
_CONDA_BASENAMES = frozenset({
    "environment.yml", "environment.yaml", "conda.yml", "conda.yaml",
})
_CONDA_PREFIX_RE = re.compile(r"^conda[_.\-].*\.ya?ml$", re.IGNORECASE)
_REQUIREMENTS_BASENAMES = frozenset({
    "requirements.txt", "requirements-dev.txt", "requirements_dev.txt",
})

# CSV output field order
CSV_FIELDS: List[str] = [
    "full_name", "html_url", "description",
    "stars", "forks_count", "open_issues", "watchers",
    "created_at", "updated_at", "pushed_at",
    "is_fork", "is_archived",
    "top_language", "languages_json", "topics",
    "repo_type", "study_score", "method_score",
    "journal_detected", "journal_confidence", "paper_year",
    "has_notebook", "notebook_count",
    "has_rmd", "rmd_count",
    "has_snakemake", "snakemake_count",
    "has_nextflow", "nextflow_count",
    "has_docker", "docker_count",
    "has_conda_env", "conda_env_count",
    "has_requirements", "requirements_count",
    "readme_snippet",
]


# ---------------------------------------------------------------------------
# HTTP / API helpers
# ---------------------------------------------------------------------------

def make_session(token: Optional[str] = None) -> requests.Session:
    """Create a requests.Session with auth and User-Agent pre-set."""
    session = requests.Session()
    session.headers.update({
        "User-Agent": USER_AGENT,
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    })
    if token:
        session.headers["Authorization"] = f"Bearer {token}"
    return session


def _wait_for_rate_limit(response: requests.Response) -> None:
    """Sleep until rate-limit resets if budget is exhausted."""
    remaining = response.headers.get("X-RateLimit-Remaining")
    reset_ts = response.headers.get("X-RateLimit-Reset")
    if remaining is not None and int(remaining) < 5:
        if reset_ts:
            wait = max(1.0, int(reset_ts) - time.time() + 2.0)
            logging.warning(
                "Rate limit exhausted (%s remaining). Sleeping %.0fs until reset.",
                remaining, wait,
            )
            time.sleep(wait)
        else:
            logging.warning(
                "Rate limit nearly exhausted (%s remaining).", remaining
            )


def api_get(
    session: requests.Session,
    url: str,
    params: Optional[Dict[str, Any]] = None,
    max_retries: int = 5,
) -> Optional[Any]:
    """
    GET a GitHub API URL with retry / exponential back-off.

    Returns parsed JSON on success, or None on permanent failure.
    Handles 202 (async compute), 403/429 (rate limits), and 5xx errors.
    """
    delay = 2.0
    for attempt in range(max_retries):
        try:
            resp = session.get(url, params=params, timeout=30)
        except requests.RequestException as exc:
            logging.warning(
                "Request error (attempt %d/%d): %s", attempt + 1, max_retries, exc
            )
            time.sleep(delay)
            delay = min(delay * 2, 60.0)
            continue

        _wait_for_rate_limit(resp)

        if resp.status_code == 200:
            return resp.json()

        if resp.status_code == 202:
            # GitHub is computing stats asynchronously — retry after a pause
            logging.debug("202 Accepted for %s — waiting %.0fs", url, delay)
            time.sleep(delay)
            delay = min(delay * 2, 30.0)
            continue

        if resp.status_code == 404:
            logging.debug("404 Not Found: %s", url)
            return None

        if resp.status_code in (403, 429):
            retry_after = resp.headers.get("Retry-After")
            reset_ts = resp.headers.get("X-RateLimit-Reset")
            if retry_after:
                wait = float(retry_after) + 1.0
            elif reset_ts:
                wait = max(2.0, int(reset_ts) - time.time() + 2.0)
            else:
                wait = delay
            logging.warning(
                "HTTP %d on %s — sleeping %.0fs before retry",
                resp.status_code, url, wait,
            )
            time.sleep(wait)
            delay = min(delay * 2, 120.0)
            continue

        if resp.status_code in (500, 502, 503, 504):
            logging.warning(
                "HTTP %d on %s (attempt %d/%d) — retrying in %.0fs",
                resp.status_code, url, attempt + 1, max_retries, delay,
            )
            time.sleep(delay)
            delay = min(delay * 2, 120.0)
            continue

        logging.warning("HTTP %d on %s — not retrying", resp.status_code, url)
        return None

    logging.error("Exhausted %d retries for %s", max_retries, url)
    return None


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def search_repositories(
    session: requests.Session,
    query: str,
    include_forks: bool = False,
    max_repos: int = 200,
) -> List[Dict[str, Any]]:
    """
    Paginate /search/repositories for *query*.

    Returns a list of raw GitHub repo objects.  By default forks are
    excluded by appending ``fork:false`` to the query.
    """
    if not include_forks and "fork:" not in query:
        query = query + " fork:false"

    results: List[Dict[str, Any]] = []
    page = 1
    per_page = 30  # keep well below GitHub search page limits

    while len(results) < max_repos:
        data = api_get(
            session,
            f"{API_BASE}/search/repositories",
            params={
                "q": query,
                "sort": "stars",
                "order": "desc",
                "per_page": per_page,
                "page": page,
            },
        )
        if data is None:
            break

        items = data.get("items", [])
        if not items:
            break

        results.extend(items)
        logging.info("  Page %d — %d repos accumulated", page, len(results))

        # GitHub search allows 30 requests/min; 2 s pause keeps us safe
        time.sleep(2)

        if len(items) < per_page:
            break  # reached last page

        page += 1

    return results[:max_repos]


# ---------------------------------------------------------------------------
# Language enrichment
# ---------------------------------------------------------------------------

def fetch_languages(
    session: requests.Session,
    full_name: str,
) -> Dict[str, int]:
    """Return byte counts per language via /repos/{owner}/{repo}/languages."""
    data = api_get(session, f"{API_BASE}/repos/{full_name}/languages")
    return data if isinstance(data, dict) else {}


# ---------------------------------------------------------------------------
# Recursive tree scan
# ---------------------------------------------------------------------------

def fetch_tree(
    session: requests.Session,
    full_name: str,
) -> Optional[List[Dict[str, str]]]:
    """
    Fetch the full recursive git tree for the default branch HEAD.

    Returns a list of tree entry dicts, or None on failure.
    """
    repo_data = api_get(session, f"{API_BASE}/repos/{full_name}")
    if not repo_data:
        return None
    branch = repo_data.get("default_branch", "HEAD")

    data = api_get(
        session,
        f"{API_BASE}/repos/{full_name}/git/trees/{branch}",
        params={"recursive": "1"},
    )
    if data is None:
        return None

    if data.get("truncated"):
        logging.warning("Git tree truncated for %s (very large repo)", full_name)

    return data.get("tree", [])


def classify_tree(tree: List[Dict[str, str]]) -> Dict[str, Any]:
    """
    Walk tree entries and count artifact types.

    Returns a dict with ``{has_X, X_count}`` keys for each artifact type.
    """
    counts: Dict[str, int] = {
        "notebook_count": 0,
        "rmd_count": 0,
        "snakemake_count": 0,
        "nextflow_count": 0,
        "docker_count": 0,
        "conda_env_count": 0,
        "requirements_count": 0,
    }

    for entry in tree:
        if entry.get("type") != "blob":
            continue
        path = entry.get("path", "")
        basename = os.path.basename(path).lower()
        _, ext = os.path.splitext(basename)

        if ext in _NOTEBOOK_EXTS:
            counts["notebook_count"] += 1
        elif ext in _RMD_EXTS:
            counts["rmd_count"] += 1
        elif basename in _SNAKEMAKE_BASENAMES or ext in _SNAKEMAKE_EXTS:
            counts["snakemake_count"] += 1
        elif basename in _NEXTFLOW_BASENAMES or ext in _NEXTFLOW_EXTS:
            counts["nextflow_count"] += 1
        elif basename in _DOCKER_BASENAMES:
            counts["docker_count"] += 1
        elif basename in _CONDA_BASENAMES or _CONDA_PREFIX_RE.match(basename):
            counts["conda_env_count"] += 1
        elif basename in _REQUIREMENTS_BASENAMES:
            counts["requirements_count"] += 1

    result: Dict[str, Any] = dict(counts)
    result["has_notebook"] = counts["notebook_count"] > 0
    result["has_rmd"] = counts["rmd_count"] > 0
    result["has_snakemake"] = counts["snakemake_count"] > 0
    result["has_nextflow"] = counts["nextflow_count"] > 0
    result["has_docker"] = counts["docker_count"] > 0
    result["has_conda_env"] = counts["conda_env_count"] > 0
    result["has_requirements"] = counts["requirements_count"] > 0
    return result


# ---------------------------------------------------------------------------
# README fetch
# ---------------------------------------------------------------------------

def fetch_readme(session: requests.Session, full_name: str) -> str:
    """
    Fetch and decode the repository README (first 4 000 characters).

    Returns an empty string if the README is missing or cannot be decoded.
    """
    data = api_get(session, f"{API_BASE}/repos/{full_name}/readme")
    if not data:
        return ""
    try:
        content = data.get("content", "")
        encoding = data.get("encoding", "base64")
        if encoding == "base64":
            text = base64.b64decode(
                content.replace("\n", "")
            ).decode("utf-8", errors="replace")
        else:
            text = content
        return text[:4000]
    except Exception as exc:
        logging.debug("README decode error for %s: %s", full_name, exc)
        return ""


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

def detect_journal(text: str) -> Tuple[Optional[str], float]:
    """
    Detect the most likely journal name in *text*.

    Patterns are checked in order from most-specific to least-specific, so
    "Nature Communications" is always preferred over bare "Nature".

    Returns ``(journal_name, confidence)`` where confidence ∈ {0.5, 0.7, 0.9, 0.0}.
    """
    for regex, name, is_abbrev in _JOURNAL_COMPILED:
        if regex.search(text):
            if name in _SINGLE_WORD_JOURNALS:
                conf = 0.5
            elif is_abbrev:
                conf = 0.7
            else:
                conf = 0.9
            return name, conf
    return None, 0.0


def infer_year(text: str) -> Optional[int]:
    """
    Conservatively infer publication year from *text*.

    A year is only returned when it appears adjacent to a recognised
    publication anchor (DOI, bioRxiv, journal name, "published", etc.).
    Standalone dates in README prose are ignored.
    """
    current_year = datetime.now(tz=timezone.utc).year
    for match in _YEAR_IN_CONTEXT_RE.finditer(text):
        raw = match.group(1) or match.group(2)
        if raw:
            year = int(raw)
            if MIN_VALID_YEAR <= year <= current_year + 1:
                return year
    return None


def score_repo(repo: Dict[str, Any], readme: str) -> Tuple[int, int]:
    """
    Compute ``(study_score, method_score)`` from repo metadata and README.

    Scores are raw keyword-frequency counts used to distinguish study
    reproducibility repositories from general software/tool repositories.
    """
    corpus = " ".join([
        repo.get("name", ""),
        repo.get("description") or "",
        " ".join(repo.get("topics", [])),
        readme[:1000],
    ]).lower()

    study_score = sum(corpus.count(kw.lower()) for kw in _STUDY_KEYWORDS)
    method_score = sum(corpus.count(kw.lower()) for kw in _METHOD_KEYWORDS)
    return study_score, method_score


def classify_repo_type(study_score: int, method_score: int) -> str:
    """
    Map ``(study_score, method_score)`` to a repo-type label.

    Labels:
      * ``"study"``   — ratio ≥ 0.65 study
      * ``"method"``  — ratio ≤ 0.35 study  (≥ 0.65 method)
      * ``"mixed"``   — between those thresholds
      * ``"unknown"`` — both scores are zero
    """
    if study_score == 0 and method_score == 0:
        return "unknown"
    ratio = study_score / (study_score + method_score)
    if ratio >= 0.65:
        return "study"
    if ratio <= 0.35:
        return "method"
    return "mixed"


# ---------------------------------------------------------------------------
# Enrichment pipeline
# ---------------------------------------------------------------------------

_TREE_SCAN_DEFAULTS: Dict[str, Any] = {
    "notebook_count": 0,
    "rmd_count": 0,
    "snakemake_count": 0,
    "nextflow_count": 0,
    "docker_count": 0,
    "conda_env_count": 0,
    "requirements_count": 0,
    "has_notebook": False,
    "has_rmd": False,
    "has_snakemake": False,
    "has_nextflow": False,
    "has_docker": False,
    "has_conda_env": False,
    "has_requirements": False,
}


def enrich_repo(
    session: requests.Session,
    repo: Dict[str, Any],
    do_tree_scan: bool = True,
) -> Dict[str, Any]:
    """
    Enrich a raw GitHub search repo item with extra metadata.

    Returns a flat dict ready for CSV/JSONL serialisation.
    """
    full_name: str = repo["full_name"]
    logging.info("  Enriching %s", full_name)

    readme = fetch_readme(session, full_name)
    languages = fetch_languages(session, full_name)

    top_language: Optional[str] = (
        max(languages, key=lambda k: languages[k]) if languages
        else repo.get("language")
    )

    tree_info: Dict[str, Any]
    if do_tree_scan:
        tree = fetch_tree(session, full_name)
        tree_info = classify_tree(tree) if tree is not None else dict(_TREE_SCAN_DEFAULTS)
    else:
        tree_info = dict(_TREE_SCAN_DEFAULTS)

    combined_text = (repo.get("description") or "") + " " + readme
    journal, journal_conf = detect_journal(combined_text)
    year = infer_year(combined_text)
    study_score, method_score = score_repo(repo, readme)
    repo_type = classify_repo_type(study_score, method_score)

    record: Dict[str, Any] = {
        # Identity
        "full_name": full_name,
        "html_url": repo.get("html_url", ""),
        "description": repo.get("description") or "",
        # Popularity
        "stars": repo.get("stargazers_count", 0),
        "forks_count": repo.get("forks_count", 0),
        "open_issues": repo.get("open_issues_count", 0),
        "watchers": repo.get("watchers_count", 0),
        # Dates
        "created_at": repo.get("created_at", ""),
        "updated_at": repo.get("updated_at", ""),
        "pushed_at": repo.get("pushed_at", ""),
        # Fork / archive status
        "is_fork": repo.get("fork", False),
        "is_archived": repo.get("archived", False),
        # Language
        "top_language": top_language or "",
        "languages_json": json.dumps(languages) if languages else "{}",
        # Topics
        "topics": ",".join(repo.get("topics", [])),
        # Repo classification
        "repo_type": repo_type,
        "study_score": study_score,
        "method_score": method_score,
        # Journal / paper
        "journal_detected": journal or "",
        "journal_confidence": round(journal_conf, 2),
        "paper_year": year or "",
        # README snippet
        "readme_snippet": readme[:500].replace("\n", " "),
    }
    record.update(tree_info)
    return record


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_outputs(records: List[Dict[str, Any]], prefix: str) -> None:
    """Write *records* to ``{prefix}.csv`` and ``{prefix}.jsonl``."""
    csv_path = f"{prefix}.csv"
    jsonl_path = f"{prefix}.jsonl"

    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    logging.info("Wrote %d rows → %s", len(records), csv_path)

    with open(jsonl_path, "w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    logging.info("Wrote %d lines → %s", len(records), jsonl_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--token",
        default=os.environ.get("GITHUB_TOKEN", ""),
        help="GitHub personal access token (or GITHUB_TOKEN env var)",
    )
    parser.add_argument(
        "--queries",
        nargs="+",
        default=DEFAULT_QUERIES,
        metavar="QUERY",
        help="One or more GitHub search query strings",
    )
    parser.add_argument(
        "--output-prefix",
        default="scrna_survey",
        metavar="PREFIX",
        help="Prefix for output .csv/.jsonl files [scrna_survey]",
    )
    parser.add_argument(
        "--max-repos",
        type=int,
        default=200,
        metavar="N",
        help="Maximum repositories per query [200]",
    )
    parser.add_argument(
        "--include-forks",
        action="store_true",
        default=False,
        help="Include forked repositories (excluded by default)",
    )
    parser.add_argument(
        "--no-tree-scan",
        action="store_true",
        default=False,
        help="Skip recursive git-tree scan (faster but less detail)",
    )
    parser.add_argument(
        "--min-stars",
        type=int,
        default=0,
        metavar="N",
        help="Minimum star count [0]",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Enable debug logging",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    if not args.token:
        logging.warning(
            "No GitHub token provided. Unauthenticated requests are "
            "rate-limited to 10 search requests/min and 60 requests/hour."
        )

    session = make_session(args.token or None)

    # --- Collect candidates from all queries, deduplicated by full_name ---
    seen: Dict[str, Dict[str, Any]] = {}
    for query in args.queries:
        logging.info("Searching: %r", query)
        items = search_repositories(
            session,
            query,
            include_forks=args.include_forks,
            max_repos=args.max_repos,
        )
        for item in items:
            seen.setdefault(item["full_name"], item)
        logging.info("Unique candidates so far: %d", len(seen))
        time.sleep(2)

    # --- Star filter ---
    candidates = [
        r for r in seen.values()
        if r.get("stargazers_count", 0) >= args.min_stars
    ]
    logging.info(
        "After star filter (>=%d stars): %d candidates",
        args.min_stars, len(candidates),
    )

    # --- Enrich each repo ---
    do_tree = not args.no_tree_scan
    records: List[Dict[str, Any]] = []
    for i, repo in enumerate(candidates, start=1):
        logging.info("[%d/%d] %s", i, len(candidates), repo["full_name"])
        try:
            record = enrich_repo(session, repo, do_tree_scan=do_tree)
            records.append(record)
        except Exception as exc:
            logging.error(
                "Error enriching %s: %s",
                repo["full_name"], exc,
                exc_info=args.verbose,
            )
        # Polite pause between repo enrichments
        time.sleep(1)

    # --- Write outputs ---
    if records:
        write_outputs(records, args.output_prefix)
    else:
        logging.warning("No records to write.")

    logging.info("Done. %d repositories saved.", len(records))
    return 0


if __name__ == "__main__":
    sys.exit(main())
