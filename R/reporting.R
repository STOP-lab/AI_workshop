render_pdf_report <- function(cfg, sample_table, results) {
  require_package("rmarkdown")

  template <- resolve_path("reports/pipeline_report.Rmd", getwd(), must_work = TRUE)
  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "report"))
  report_path <- file.path(out_dir, "point_seq_splicing_report.pdf")

  render_env <- new.env(parent = globalenv())
  render_env$cfg <- cfg
  render_env$sample_table <- sample_table
  render_env$results <- results

  rmarkdown::render(
    input = template,
    output_format = "pdf_document",
    output_file = basename(report_path),
    output_dir = out_dir,
    params = list(
      output_dir = project_output_dir(cfg),
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    envir = render_env,
    quiet = TRUE
  )

  report_path
}
