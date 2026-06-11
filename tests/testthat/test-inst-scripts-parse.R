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
  # 4 canonical example scripts after the protocol-pair demos
  # (03_supervised_protocol_legacy.R / 04_supervised_protocol_nested_refit.R)
  # were removed: there is one training procedure, no training_protocol choice.
  expect_gte(length(rfiles), 4L)
  for (f in rfiles) {
    expect_error(parse(file = f), NA,
                 info = paste("parse failed:", basename(f)))
  }
})
