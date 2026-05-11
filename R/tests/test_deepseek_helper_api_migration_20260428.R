#!/usr/bin/env Rscript

source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R")

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

assert_identical <- function(x, y, msg) {
  if (!identical(x, y)) {
    stop(sprintf("%s\nExpected: %s\nActual: %s", msg, paste(capture.output(str(y)), collapse = " "), paste(capture.output(str(x)), collapse = " ")), call. = FALSE)
  }
}

chat_cfg <- tc_deepseek_resolve_model("deepseek-chat")
assert_identical(chat_cfg$model, "deepseek-v4-flash", "deepseek-chat should map to deepseek-v4-flash")
assert_true(!isTRUE(chat_cfg$thinking_enabled), "deepseek-chat alias should not enable thinking")

reasoner_cfg <- tc_deepseek_resolve_model("deepseek-reasoner")
assert_identical(reasoner_cfg$model, "deepseek-v4-flash", "deepseek-reasoner should map to deepseek-v4-flash")
assert_true(isTRUE(reasoner_cfg$thinking_enabled), "deepseek-reasoner alias should enable thinking")

explicit_cfg <- tc_deepseek_resolve_model("deepseek-v4-pro")
assert_identical(explicit_cfg$model, "deepseek-v4-pro", "explicit v4 model should be preserved")

endpoint <- tc_deepseek_chat_endpoint("https://api.deepseek.com/")
assert_identical(endpoint, "https://api.deepseek.com/chat/completions", "chat endpoint should use official OpenAI-compatible path")

payload_reasoner <- tc_deepseek_build_payload(
  prompt = "hello",
  model = "deepseek-reasoner",
  system_prompt = "You are a helpful assistant."
)
assert_identical(payload_reasoner$model, "deepseek-v4-flash", "reasoner payload should use mapped v4 model")
assert_true(is.list(payload_reasoner$thinking), "reasoner payload should include thinking config")
assert_identical(payload_reasoner$thinking$type, "enabled", "thinking config should be enabled")
assert_true(length(payload_reasoner$messages) == 2L, "payload should contain system and user messages")
assert_identical(payload_reasoner$messages[[1]]$role, "system", "first message should be system")
assert_identical(payload_reasoner$messages[[2]]$role, "user", "second message should be user")

payload_chat <- tc_deepseek_build_payload(
  prompt = "hello",
  model = "deepseek-chat"
)
assert_true(is.null(payload_chat$thinking), "non-thinking payload should not include thinking config")
assert_identical(payload_chat$model, "deepseek-v4-flash", "chat payload should also use mapped v4 model")

cfg <- tc_deepseek_request_config(
  prompt = "hello",
  model = "deepseek-reasoner",
  api_key = "sk-test-1234567890",
  base_url = "https://api.deepseek.com"
)
assert_identical(cfg$url, "https://api.deepseek.com/chat/completions", "request config should build correct URL")
assert_identical(cfg$headers[["Authorization"]], "Bearer sk-test-1234567890", "request config should carry bearer token")
assert_identical(cfg$payload$model, "deepseek-v4-flash", "request config payload should use mapped model")

cat("DeepSeek helper API migration tests passed.\n")
