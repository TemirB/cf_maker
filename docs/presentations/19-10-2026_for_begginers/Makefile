XELATEX ?= xelatex
PYTHON ?= python3
STEM := 19-10-2026_for_begginers
SLIDES := presentation_$(STEM)
NOTES := speaker_notes_$(STEM)

.PHONY: all notes figures
all: notes
	mkdir -p build
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build $(SLIDES).tex
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build $(SLIDES).tex
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build $(NOTES).tex
	cp build/$(SLIDES).pdf $(SLIDES).pdf
	cp build/$(NOTES).pdf $(NOTES).pdf

notes:
	$(PYTHON) build_notes.py

figures:
	$(PYTHON) make_figures.py
	$(PYTHON) extract_research_figures.py
