# Makefile for Mass Gap Paper
# Usage:
#   make          - Build the PDF
#   make clean    - Remove generated files
#   make watch    - Continuously rebuild on changes (requires latexmk)

PAPER = mass_gap
LATEX = pdflatex
BIBTEX = bibtex

# Default target
all: $(PAPER).pdf

TABLE_DEPS = $(wildcard tables/*.tex)
FIGURE_DEPS = $(wildcard figures/*.png figures/*.pdf)

# Build PDF (run twice for references, with bibtex in between)
$(PAPER).pdf: $(PAPER).tex $(PAPER).bib aastex701.cls aasjournalv7.bst $(TABLE_DEPS) $(FIGURE_DEPS)
	$(LATEX) $(PAPER)
	$(BIBTEX) $(PAPER)
	$(LATEX) $(PAPER)
	$(LATEX) $(PAPER)

# Quick build (single pass, no bibliography update)
quick:
	$(LATEX) $(PAPER)

# Clean generated files
clean:
	rm -f $(PAPER).aux $(PAPER).bbl $(PAPER).blg $(PAPER).log \
	      $(PAPER).out $(PAPER).toc $(PAPER).fls $(PAPER).fdb_latexmk \
	      $(PAPER).synctex.gz \
	      figures/* tables/*

# Clean everything including PDF
distclean: clean
	rm -f $(PAPER).pdf

# Continuous build (requires latexmk)
watch:
	latexmk -pdf -pvc $(PAPER)

# Generate tables from data (placeholder - to be expanded)
tables:
	@echo "Table generation scripts not yet implemented"
	@echo "Run: python scripts/tables/tab_summary.py"

.PHONY: all quick clean distclean watch tables
