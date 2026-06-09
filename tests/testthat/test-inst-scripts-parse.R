# Guard: every canonical supervised usage script under inst/scripts/ must parse
# (syntax-only). These are templates with <PATH_TO_...> placeholders and a
# RUN <- FALSE guard, so we only check they PARSE — we never source/run them.
test_that("all inst/scripts/*.R parse without syntax error", {
  scripts_dir <- system.file("scripts", package = "OtsuFire")
  # During pkgload::load_all() inst/ is mapped, so system.file resolves it.
  if (!nzchar(scripts_dir) || !dir.exists(scripts_dir)) {
    skip("inst/scripts not resolvable in this build context")
  }
  rfiles <- list.files(scripts_dir, pattern = "[.]R$", full.names = TRUE)
  expect_gte(length(rfiles), 6L)
  for (f in rfiles) {
    expect_error(parse(file = f), NA,
                 info = paste("parse failed:", basename(f)))
  }
})
