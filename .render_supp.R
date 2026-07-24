suppressMessages(library(pdftools))
outdir <- file.path(getwd(), ".paper_pages")
dir.create(outdir, showWarnings = FALSE)
info <- pdf_info("paper_original.pdf")
cat("pages:", info$pages, "\n")
for (pg in 19:min(info$pages, 45)) {
  out <- file.path(outdir, sprintf("page_%02d.png", pg))
  pdf_convert("paper_original.pdf", format = "png", pages = pg,
              dpi = 170, filenames = out, verbose = FALSE)
}
cat("done\n")
