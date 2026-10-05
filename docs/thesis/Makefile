LATEXMK ?= latexmk
GO ?= go
CONFIG ?= config/kt.json
RESULTS ?=
PRESENTATION_PLOTS ?= 12
BUILD_DIR ?= build
STRICT_ASSETS ?= 0

# abspath treats spaces as separators; preserve the output directory as one path.
DOCS_BUILD_DIR := $(if $(filter /%,$(firstword $(BUILD_DIR))),$(BUILD_DIR),$(CURDIR)/$(BUILD_DIR))
REPORT_BUILD_DIR := $(DOCS_BUILD_DIR)/report
PRESENTATION_BUILD_DIR := $(DOCS_BUILD_DIR)/presentation
PRESENTATION_BUILDER ?= $(PRESENTATION_BUILD_DIR)/presentation-builder
LATEXMK_FLAGS := -cd -xelatex -interaction=nonstopmode -halt-on-error -file-line-error

# Force a new TeX pass: cached draft output must not bypass strict asset checks.
ifeq ($(STRICT_ASSETS),1)
STRICT_FLAGS := -g -usepretex='\def\ThesisStrictAssets{1}'
endif

.DEFAULT_GOAL := all
.NOTPARALLEL:
.PHONY: all report presentation clean doctor help check-tools check-latexmk \
        check-go check-presentation-tools

all: report presentation

help:
	@printf '%s\n' \
	  'make                     Собрать текст работы и презентацию' \
	  'make report              Собрать build/report/thesis.pdf' \
	  'make presentation        Собрать PDF со свежими результатами CONFIG' \
	  'make clean               Удалить PDF и вспомогательные файлы через latexmk' \
	  'make doctor              Проверить Go, latexmk, XeLaTeX и Biber' \
	  '' \
	  'Настройки: BUILD_DIR=/path/to/output LATEXMK=/path/to/latexmk' \
	  'CONFIG=config/kt.json     Конфиг расчёта; относительный путь от корня репы' \
	  'RESULTS=/path/to/results  Готовый расчёт вместо выбора через CONFIG' \
	  'GO=go                     Go compiler для presentation-builder' \
	  'PRESENTATION_PLOTS=12     Максимум графиков в PDF; 0 — все' \
	  'STRICT_ASSETS=1           Прервать сборку, если отсутствуют исходные рисунки'

check-latexmk:
	@command -v "$(LATEXMK)" >/dev/null 2>&1 || { \
	  printf '%s\n' 'Не найден latexmk. Установите TeX Live и latexmk или задайте LATEXMK=/path/to/latexmk.' >&2; \
	  exit 1; \
	}

check-tools: check-latexmk
	@status=0; \
	for tool in xelatex biber; do \
	  if ! command -v "$$tool" >/dev/null 2>&1; then \
	    printf 'Не найден %s. Требуются XeLaTeX и Biber; см. docs/thesis/README.md.\n' "$$tool" >&2; \
	    status=1; \
	  fi; \
	done; \
	exit "$$status"

check-go:
	@command -v "$(GO)" >/dev/null 2>&1 || { \
	  printf '%s\n' 'Не найден Go. Установите Go или задайте GO=/path/to/go.' >&2; \
	  exit 1; \
	}

check-presentation-tools: check-tools check-go

report: check-tools
	mkdir -p "$(REPORT_BUILD_DIR)"
	"$(LATEXMK)" $(LATEXMK_FLAGS) $(STRICT_FLAGS) \
	  -outdir="$(REPORT_BUILD_DIR)" -jobname=thesis report/compile.tex

# XeTeX's image loader also needs the output tree in its explicit search path.
$(PRESENTATION_BUILDER): ../../tools/presentation_builder/main.go
	mkdir -p "$(dir $@)"
	"$(GO)" build -o "$@" "$<"

presentation: check-presentation-tools $(PRESENTATION_BUILDER)
	mkdir -p "$(PRESENTATION_BUILD_DIR)"
	"$(PRESENTATION_BUILDER)" -config "$(CONFIG)" -repo-root "$(CURDIR)/../.." \
	  -source-dir "$(CURDIR)" -build-dir "$(PRESENTATION_BUILD_DIR)" \
	  -latexmk "$(LATEXMK)" -plots "$(PRESENTATION_PLOTS)" \
	  -strict-assets="$(if $(filter 1,$(STRICT_ASSETS)),true,false)" \
	  $(if $(strip $(RESULTS)),-results "$(RESULTS)")

clean: check-latexmk
	"$(LATEXMK)" -cd -C -outdir="$(REPORT_BUILD_DIR)" -jobname=thesis report/compile.tex
	"$(LATEXMK)" -cd -C -outdir="$(PRESENTATION_BUILD_DIR)" -jobname=presentation presentation/main.tex

doctor: check-presentation-tools
	"$(GO)" version
	"$(LATEXMK)" -v
	xelatex --version
	biber --version
