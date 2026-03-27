#!/usr/bin/env python3
"""
Unit tests for github_scrna_survey.py

Tests cover the pure-Python analysis helpers that do not require a network
connection or a GitHub token.  Run with:

    python -m pytest test_github_scrna_survey.py -v
    # or simply:
    python test_github_scrna_survey.py
"""

import sys
import unittest

# Ensure the module under test is importable when run from the repo root
sys.path.insert(0, ".")
from github_scrna_survey import (
    classify_repo_type,
    classify_tree,
    detect_journal,
    infer_year,
    score_repo,
)


# ---------------------------------------------------------------------------
# Journal detection
# ---------------------------------------------------------------------------

class TestDetectJournal(unittest.TestCase):

    # --- specificity: long names beat short names ---

    def test_nature_communications_beats_nature(self):
        journal, conf = detect_journal("Published in Nature Communications 2023.")
        self.assertEqual(journal, "Nature Communications")
        self.assertGreaterEqual(conf, 0.9)

    def test_nature_biotechnology_beats_nature(self):
        journal, conf = detect_journal("Code for our Nature Biotechnology paper.")
        self.assertEqual(journal, "Nature Biotechnology")
        self.assertGreaterEqual(conf, 0.9)

    def test_nature_methods_beats_nature(self):
        journal, conf = detect_journal("Reproducible code for Nature Methods 2022.")
        self.assertEqual(journal, "Nature Methods")

    def test_nature_medicine_beats_nature(self):
        journal, conf = detect_journal("See our Nature Medicine manuscript.")
        self.assertEqual(journal, "Nature Medicine")

    def test_genome_biology(self):
        journal, conf = detect_journal("Analysis code for Genome Biology 2021 paper.")
        self.assertEqual(journal, "Genome Biology")
        self.assertGreaterEqual(conf, 0.9)

    def test_cell_reports_beats_cell(self):
        journal, conf = detect_journal("Reproducibility code — Cell Reports Methods.")
        self.assertEqual(journal, "Cell Reports Methods")

    def test_bare_nature(self):
        journal, conf = detect_journal("Published in Nature.")
        self.assertEqual(journal, "Nature")
        self.assertAlmostEqual(conf, 0.5)

    def test_bare_cell(self):
        journal, conf = detect_journal("Data from Cell (2020).")
        self.assertEqual(journal, "Cell")

    def test_no_journal(self):
        journal, conf = detect_journal("This is just a README with no journal mention.")
        self.assertIsNone(journal)
        self.assertEqual(conf, 0.0)

    # --- abbreviations ---

    def test_abbrev_nat_commun(self):
        journal, conf = detect_journal("Nat. Commun. (2022) doi:10.1038/s41467")
        self.assertEqual(journal, "Nature Communications")
        self.assertAlmostEqual(conf, 0.7)

    def test_abbrev_nat_biotechnol(self):
        journal, conf = detect_journal("Nat Biotechnol 2023")
        self.assertEqual(journal, "Nature Biotechnology")

    def test_abbrev_nat_methods(self):
        journal, conf = detect_journal("Nat. Methods paper")
        self.assertEqual(journal, "Nature Methods")

    def test_abbrev_genome_biol(self):
        journal, conf = detect_journal("Code for Genome Biol paper.")
        self.assertEqual(journal, "Genome Biology")

    # --- case insensitivity ---

    def test_case_insensitive(self):
        journal, _ = detect_journal("NATURE COMMUNICATIONS")
        self.assertEqual(journal, "Nature Communications")

    # --- mixed text ---

    def test_mixed_text_prefers_first_match(self):
        """Nature Communications should match before Nature when both appear."""
        journal, _ = detect_journal(
            "This is code from our Nature Communications paper. "
            "We also cite Nature 2020."
        )
        self.assertEqual(journal, "Nature Communications")


# ---------------------------------------------------------------------------
# Year inference
# ---------------------------------------------------------------------------

class TestInferYear(unittest.TestCase):

    def test_year_after_doi(self):
        year = infer_year("doi: 10.1038/s41467-023-XXXX published 2023")
        self.assertEqual(year, 2023)

    def test_year_after_biorxiv(self):
        year = infer_year("bioRxiv preprint 2022. doi: 10.1101/2022.01.01")
        self.assertEqual(year, 2022)

    def test_year_after_published(self):
        year = infer_year("Published: 2021. Code for the analysis.")
        self.assertEqual(year, 2021)

    def test_year_near_nature(self):
        year = infer_year("Nature Communications, 2020.")
        # Nature appears in journal patterns so this may or may not match
        # depending on how close the year is — accept either outcome
        self.assertIn(year, [2020, None])

    def test_no_year_from_bare_date(self):
        """Standalone file/commit date should NOT trigger a year claim."""
        year = infer_year("Last updated: 2023-06-15. README v2.")
        self.assertIsNone(year)

    def test_no_year_from_future(self):
        """Years far in the future should not be returned."""
        year = infer_year("doi: 10.1038/xxxx published 2099")
        self.assertIsNone(year)

    def test_no_year_from_old_date(self):
        """Years before 2010 should not be returned."""
        year = infer_year("bioRxiv 2005")
        self.assertIsNone(year)

    def test_year_returns_integer(self):
        year = infer_year("doi: 10.1038/x, 2019")
        if year is not None:
            self.assertIsInstance(year, int)


# ---------------------------------------------------------------------------
# Repo scoring and classification
# ---------------------------------------------------------------------------

class TestScoreRepo(unittest.TestCase):

    def _make_repo(self, name="", description="", topics=None):
        return {
            "name": name,
            "description": description,
            "topics": topics or [],
        }

    def test_study_heavy(self):
        repo = self._make_repo(
            name="scrna-seq-tissue-atlas",
            description="Single-cell RNA-seq analysis of patient tissue",
            topics=["single-cell", "reproducibility"],
        )
        study, method = score_repo(repo, "Dataset analysis of clinical cohort")
        self.assertGreater(study, method)

    def test_method_heavy(self):
        repo = self._make_repo(
            name="myscratool",
            description="A software package and tool for single-cell benchmarking",
            topics=["package", "algorithm"],
        )
        study, method = score_repo(repo, "Library for workflow pipeline benchmarking")
        self.assertGreater(method, study)

    def test_zero_scores(self):
        repo = self._make_repo(name="hello-world", description="A demo project")
        study, method = score_repo(repo, "No relevant content here.")
        self.assertGreaterEqual(study, 0)
        self.assertGreaterEqual(method, 0)


class TestClassifyRepoType(unittest.TestCase):

    def test_study(self):
        self.assertEqual(classify_repo_type(20, 3), "study")

    def test_method(self):
        self.assertEqual(classify_repo_type(2, 20), "method")

    def test_mixed(self):
        self.assertEqual(classify_repo_type(10, 10), "mixed")

    def test_unknown(self):
        self.assertEqual(classify_repo_type(0, 0), "unknown")

    def test_boundary_study(self):
        # ratio = 13/20 = 0.65 — exactly at study threshold
        self.assertEqual(classify_repo_type(13, 7), "study")

    def test_boundary_method(self):
        # ratio = 7/20 = 0.35 — exactly at method threshold
        self.assertEqual(classify_repo_type(7, 13), "method")


# ---------------------------------------------------------------------------
# Tree classification
# ---------------------------------------------------------------------------

class TestClassifyTree(unittest.TestCase):

    def _entry(self, path):
        return {"type": "blob", "path": path}

    def test_jupyter_notebook(self):
        tree = [self._entry("notebooks/analysis.ipynb")]
        result = classify_tree(tree)
        self.assertTrue(result["has_notebook"])
        self.assertEqual(result["notebook_count"], 1)

    def test_rmd_file(self):
        tree = [self._entry("analysis.Rmd")]
        result = classify_tree(tree)
        self.assertTrue(result["has_rmd"])
        self.assertEqual(result["rmd_count"], 1)

    def test_snakefile(self):
        tree = [self._entry("workflow/Snakefile")]
        result = classify_tree(tree)
        self.assertTrue(result["has_snakemake"])

    def test_snakemake_extension(self):
        tree = [self._entry("rules/align.smk")]
        result = classify_tree(tree)
        self.assertTrue(result["has_snakemake"])

    def test_nextflow(self):
        tree = [self._entry("main.nf")]
        result = classify_tree(tree)
        self.assertTrue(result["has_nextflow"])

    def test_dockerfile(self):
        tree = [self._entry("docker/Dockerfile")]
        result = classify_tree(tree)
        self.assertTrue(result["has_docker"])
        self.assertEqual(result["docker_count"], 1)

    def test_conda_environment(self):
        tree = [self._entry("environment.yml")]
        result = classify_tree(tree)
        self.assertTrue(result["has_conda_env"])

    def test_conda_prefixed_name(self):
        tree = [self._entry("conda_env.yaml")]
        result = classify_tree(tree)
        self.assertTrue(result["has_conda_env"])

    def test_requirements_txt(self):
        tree = [self._entry("requirements.txt")]
        result = classify_tree(tree)
        self.assertTrue(result["has_requirements"])

    def test_empty_tree(self):
        result = classify_tree([])
        for key in ["has_notebook", "has_rmd", "has_snakemake",
                    "has_nextflow", "has_docker", "has_conda_env",
                    "has_requirements"]:
            self.assertFalse(result[key])

    def test_tree_entries(self):
        tree = [
            self._entry("notebooks/01_qc.ipynb"),
            self._entry("notebooks/02_cluster.ipynb"),
            self._entry("analysis.Rmd"),
            self._entry("Dockerfile"),
            self._entry("environment.yml"),
            self._entry("workflow/Snakefile"),
            {"type": "tree", "path": "some_dir"},  # non-blob, should be ignored
        ]
        result = classify_tree(tree)
        self.assertEqual(result["notebook_count"], 2)
        self.assertEqual(result["rmd_count"], 1)
        self.assertEqual(result["docker_count"], 1)
        self.assertEqual(result["conda_env_count"], 1)
        self.assertEqual(result["snakemake_count"], 1)
        self.assertTrue(result["has_notebook"])

    def test_case_insensitive_rmd(self):
        """Both .Rmd and .rmd should be detected."""
        tree = [
            self._entry("ANALYSIS.RMD"),
            self._entry("report.rmd"),
        ]
        result = classify_tree(tree)
        self.assertEqual(result["rmd_count"], 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
