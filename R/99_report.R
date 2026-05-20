# 99_report.R
# Quarto/RMarkdown stub. Render summary report after the benchmark and final
# model fits complete. Replace path to the .qmd once you author the template.

source(here::here("R", "00_config.R"))

render_report <- function(input = file.path(PATHS$reports, "summary.qmd"),
                          output_dir = PATHS$reports) {
  if (!file.exists(input)) {
    log_msg("No report template at ", input, " - skipping.")
    return(invisible(NULL))
  }
  quarto::quarto_render(input, output_dir = output_dir)
}

if (sys.nframe() == 0) render_report()
