# Linux defaults. A compiler override must be explicit: make CXX=clang++.
ifneq ($(origin CXX),command line)
CXX := g++
endif

CMAKE ?= cmake
CTEST ?= ctest
BUILD_DIR ?= build
BUILD_TYPE ?= Release
JOBS ?= 2
CONFIG ?= config/kt.json
ROOT_DIR ?=
CMAKE_ARGS ?=
CTEST_ARGS ?=
SANITIZE_BUILD_DIR ?= build-sanitize
FIXTURE_DIR ?= input/gaussian_make
CLANG_FORMAT ?= clang-format
CLANG_TIDY ?= clang-tidy
CLANG_TIDY_ARGS ?=
LATEXMK ?= latexmk
GO ?= go
HOST_OS := $(shell uname -s)
HOST_ARCH := $(shell uname -m)
ifeq ($(HOST_OS),Linux)
  ifeq ($(HOST_ARCH),x86_64)
    DEFAULT_PRESENTATION_BUILDER := tools/bin/presentation-builder-linux-amd64
  else ifneq ($(filter aarch64 arm64,$(HOST_ARCH)),)
    DEFAULT_PRESENTATION_BUILDER := tools/bin/presentation-builder-linux-arm64
  endif
endif
PRESENTATION_BUILDER ?= $(if $(DEFAULT_PRESENTATION_BUILDER),$(DEFAULT_PRESENTATION_BUILDER),$(BUILD_DIR)/presentation-builder)
RESULTS ?=
PRESENTATION_PLOTS ?= 12
DOCS_BUILD_DIR ?= docs/thesis/build
STRICT_ASSETS ?= 0
DOCS_OUTPUT_DIR := $(if $(filter /%,$(firstword $(DOCS_BUILD_DIR))),$(DOCS_BUILD_DIR),$(CURDIR)/$(DOCS_BUILD_DIR))

ROOT_ARGS = $(if $(ROOT_DIR),-DROOT_DIR="$(ROOT_DIR)")

.DEFAULT_GOAL := build
# Serialize wrapper targets even with make -j; CMake builds still use JOBS.
.NOTPARALLEL:
.PHONY: help configure build test run check rebuild clean doctor format format-check \
        lint sanitize-configure sanitize sanitize-test fixture demo docs thesis presentation \
        run-presentation check-presentation-tools docs-clean

help:
	@printf '%s\n' \
	  'make                     Release-сборка через g++ (по умолчанию JOBS=2)' \
	  'make test                Собрать и выполнить CTest' \
	  'make run                 Собрать и запустить CONFIG=config/kt.json' \
	  'make check               Собрать, проверить тесты, запустить CONFIG' \
	  'make configure           Восстановить файлы CMake, включая VerifyGlobs.cmake' \
	  'make rebuild             Заново собрать с --clean-first' \
	  'make clean               Удалить результаты компиляции через CMake' \
	  'make doctor              Показать компилятор, CMake и активный ROOT' \
	  'make format              Отформатировать src/tests/tools' \
	  'make format-check        Проверить формат без изменений' \
	  'make lint                Запустить clang-tidy по compile_commands.json' \
	  'make sanitize            ASan/UBSan-сборка приложения и тестов в build-sanitize' \
	  'make sanitize-test       Собрать с ASan/UBSan и выполнить CTest' \
	  'make fixture             Создать Gaussian-вход в FIXTURE_DIR (без перезаписи)' \
	  'make demo                Создать/использовать Gaussian-вход и запустить анализ' \
	  'make docs                Собрать PDF текста диплома и презентации' \
	  'make thesis              Собрать PDF текста диплома' \
	  'make presentation        Обновить PDF из результатов CONFIG=config/kt.json' \
	  'make run-presentation    Собрать программу, выполнить CONFIG и получить PDF' \
	  'make docs-clean          Очистить артефакты сборки LaTeX' \
	  '' \
	  'Настройки: JOBS=4 CXX=clang++ CONFIG=config/rapidity.json BUILD_DIR=build-linux' \
	  'ROOT_DIR=/path/to/root/cmake CMAKE_ARGS="..." CTEST_ARGS="-R ratios"' \
	  'make run-presentation CONFIG=config/y.json' \
	  'make presentation RESULTS=/path/to/results' \
	  'PRESENTATION_PLOTS=12     Максимум графиков в PDF; 0 — все' \
	  'Перед сборкой активируйте установленный ROOT: source /path/to/root/bin/thisroot.sh'

# Always configure first: a stale build tree may have lost VerifyGlobs.cmake.
configure:
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" \
	  -DCMAKE_CXX_COMPILER="$(CXX)" -DCF_MAKER_SANITIZE=OFF $(ROOT_ARGS) $(CMAKE_ARGS)

build: configure
	$(CMAKE) --build "$(BUILD_DIR)" --parallel "$(JOBS)"

test: build
	$(CTEST) --test-dir "$(BUILD_DIR)" --output-on-failure $(CTEST_ARGS)

run: build
	"$(BUILD_DIR)/main" "$(CONFIG)"

check: test
	"$(BUILD_DIR)/main" "$(CONFIG)"

rebuild: configure
	$(CMAKE) --build "$(BUILD_DIR)" --clean-first --parallel "$(JOBS)"

clean:
	$(CMAKE) --build "$(BUILD_DIR)" --target clean

doctor:
	"$(CXX)" --version
	$(CMAKE) --version
	root-config --version
	root-config --prefix

format:
	find src tests tools -type f \( -name '*.cpp' -o -name '*.h' -o -name '*.C' \) \
	  -exec "$(CLANG_FORMAT)" -i {} +

format-check:
	find src tests tools -type f \( -name '*.cpp' -o -name '*.h' -o -name '*.C' \) \
	  -exec "$(CLANG_FORMAT)" --dry-run --Werror {} +

lint: configure
	find src -type f -name '*.cpp' \
	  -exec "$(CLANG_TIDY)" -p "$(BUILD_DIR)" $(CLANG_TIDY_ARGS) {} +

# Use a separate, ignored directory; build-asan is tracked in the existing repo.
# Global flags instrument tests and the generator as well as the main program.
sanitize-configure:
	$(CMAKE) -S . -B "$(SANITIZE_BUILD_DIR)" -DCMAKE_BUILD_TYPE=Debug \
	  -DCMAKE_CXX_COMPILER="$(CXX)" -DCF_MAKER_SANITIZE=ON \
	  '-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer' \
	  '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined' $(ROOT_ARGS) $(CMAKE_ARGS)

sanitize: sanitize-configure
	$(CMAKE) --build "$(SANITIZE_BUILD_DIR)" --parallel "$(JOBS)"

sanitize-test: sanitize
	$(CTEST) --test-dir "$(SANITIZE_BUILD_DIR)" --output-on-failure $(CTEST_ARGS)

fixture: build
	"$(BUILD_DIR)/generate_gaussian" "$(FIXTURE_DIR)"

demo: build
	@if [ ! -f "$(FIXTURE_DIR)/config.json" ] || [ ! -f "$(FIXTURE_DIR)/input.root" ]; then \
	  "$(BUILD_DIR)/generate_gaussian" "$(FIXTURE_DIR)"; \
	fi
	"$(BUILD_DIR)/main" "$(FIXTURE_DIR)/config.json"

# Documents have their own dependencies and do not need a C++/ROOT build.
docs:
	$(MAKE) -C docs/thesis all LATEXMK="$(LATEXMK)" \
	  BUILD_DIR="$(DOCS_OUTPUT_DIR)" STRICT_ASSETS="$(STRICT_ASSETS)" \
	  GO="$(GO)" PRESENTATION_BUILDER="$(abspath $(PRESENTATION_BUILDER))" \
	  CONFIG="$(CONFIG)" RESULTS="$(RESULTS)" \
	  PRESENTATION_PLOTS="$(PRESENTATION_PLOTS)"

thesis:
	$(MAKE) -C docs/thesis report LATEXMK="$(LATEXMK)" \
	  BUILD_DIR="$(DOCS_OUTPUT_DIR)" STRICT_ASSETS="$(STRICT_ASSETS)"

presentation:
	$(MAKE) -C docs/thesis presentation LATEXMK="$(LATEXMK)" GO="$(GO)" \
	  BUILD_DIR="$(DOCS_OUTPUT_DIR)" STRICT_ASSETS="$(STRICT_ASSETS)" \
	  PRESENTATION_BUILDER="$(abspath $(PRESENTATION_BUILDER))" \
	  CONFIG="$(CONFIG)" RESULTS="$(RESULTS)" \
	  PRESENTATION_PLOTS="$(PRESENTATION_PLOTS)"

check-presentation-tools:
	$(MAKE) -C docs/thesis check-presentation-tools LATEXMK="$(LATEXMK)" GO="$(GO)" \
	  PRESENTATION_BUILDER="$(abspath $(PRESENTATION_BUILDER))"

# Diagnose missing document tools before starting the analysis. Exit code 2
# still permits a diagnostic presentation; ordinary run/check keep their status.
run-presentation: check-presentation-tools build $(PRESENTATION_BUILDER)
	"$(abspath $(PRESENTATION_BUILDER))" -config "$(CONFIG)" \
	  -repo-root "$(CURDIR)" -source-dir "$(CURDIR)/docs/thesis" \
	  -build-dir "$(DOCS_OUTPUT_DIR)/presentation" -latexmk "$(LATEXMK)" \
	  -plots "$(PRESENTATION_PLOTS)" \
	  -strict-assets="$(if $(filter 1,$(STRICT_ASSETS)),true,false)" \
	  -executable "$(abspath $(BUILD_DIR)/main)"

ifeq ($(wildcard $(PRESENTATION_BUILDER)),)
$(PRESENTATION_BUILDER): tools/presentation_builder/main.go
	mkdir -p "$(dir $@)"
	"$(GO)" build -o "$@" "$<"
endif

docs-clean:
	$(MAKE) -C docs/thesis clean LATEXMK="$(LATEXMK)" \
	  BUILD_DIR="$(DOCS_OUTPUT_DIR)"
