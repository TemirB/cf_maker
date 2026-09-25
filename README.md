# cf_maker_new

Утилита для построения и фита 3D корреляционных функций (HBT) с использованием ROOT.

Проект предназначен для batch-обработки больших ROOT-файлов и сохранения результатов
в виде отдельных выходных файлов.

---

# Структура

```
src/
├── main.cpp            — точка входа
├── config/             — загрузка и валидация JSON-конфигурации
├── core/               — общие примитивы: бининг, LCMS-оси, утилиты файловой системы
├── io/                 — чтение входных гистограмм
├── fit/                — модель фита, FitResult, начальные параметры
├── draw/               — отрисовка и стилизация
└── analysis/           — стадии анализа: 3D CF, 1D/2D проекции, отношения, зависимости
```

---

# Установка зависимостей

Общие зависимости: C++17 компилятор, CMake ≥ 3.20, git, **ROOT 6**,
OpenSSL (headers), опционально ccache / clang-format / clang-tidy.

`nlohmann/json` подключается как git-сабмодуль — после клонирования:

```bash
git submodule update --init --recursive
```

## Linux

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y build-essential cmake git libssl-dev \
    root-system ccache clang-format clang-tidy
```

Если `root-system` в репозитории слишком старый — поставьте прекомпилированный
бинарь с <https://root.cern/install/#linux-packages-and-pre-compiled-binaries>
и активируйте окружение перед сборкой:

```bash
source /path/to/root-install/bin/thisroot.sh
```

Альтернатива — conda:

```bash
conda install -c conda-forge root
```

### Fedora / RHEL

```bash
sudo dnf install -y gcc-c++ cmake git openssl-devel root ccache clang clang-tools-extra
```

### Arch Linux

```bash
sudo pacman -S gcc cmake make git openssl root ccache clang
```

## macOS

```bash
brew install cmake root openssl clang-format ccache
```

`clang-tidy` на macOS доступен только в keg-only пакете `llvm` (см. AGENTS.md).

---

# Сборка

Из корня проекта:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

После сборки исполняемый файл `main` будет находиться в директории `build/`.

# Запуск

Программа принимает один аргумент — путь к JSON-конфигурации:

```bash
./build/main config/kt.json
```

Запускать нужно из корня репозитория. Результат работы будет находиться в директории
`<base_output>/<dir>` из конфигурации.

# Генератор тестовых данных

```bash
./build/generate_gaussian input/gaussian
./build/main input/gaussian/config.json
```

Генератор создаёт `input.root` и готовый `config.json` в указанной директории.
Существующие файлы не перезаписываются. Все зарядовые группы, центральности и
бины получают одинаковую гауссову корреляционную функцию:

```text
C(q) = 1 + lambda * exp(-(R_out²*q_out² + R_side²*q_side² + R_long²*q_long²) / 0.197²)
```

Радиусы задаются в фм, q — в ГэВ/c. Cross-термы равны нулю.
Имена гистограмм совместимы с анализом: `bp_CH_CENT_num_BIN` — постоянный
знаменатель, `bp_CH_CENT_num_wei_BIN` — числитель `den * C(q)`.
Значения вычисляются в центрах бинов. Параметры генерации записываются
в объект `generator_truth` внутри ROOT-файла.

```bash
./build/generate_gaussian input/gaussian_noisy --type rapidity \
    --r-out 4 --r-side 5 --r-long 6 --lambda 0.7 \
    --bins 40 --q-max 0.3 --counts 10000 --noise --seed 42
```

По умолчанию используется `kt`, радиусы 4/5/6 фм, λ=0.7, сетка 40³
в диапазоне ±0.3 ГэВ/c и знаменатель 10000 на бин. Без `--noise` CF точная;
с `--noise` числитель получает пуассоновские флуктуации с заданным seed.
Знаменатель фиксирован, его ошибка нулевая; ошибка числителя равна √num.
Это искусственные данные для проверки фита и пайплайна. Применяемая анализом
биномиальная формула ошибок отношения не является моделью независимых
пуассоновских счётчиков, поэтому шумовой режим не предназначен для калибровки
ошибок параметров и p-value.

Генерируемый конфиг запускает все стадии в одном потоке без сохранения картинок;
результаты попадают в `<директория>/results`. Полный список опций:
`./build/generate_gaussian --help`.

# Конфигурация

Все секции, кроме `machine` и `vars`, опциональны — у них есть дефолты в коде,
поэтому минимальный конфиг выглядит так:

```json
{
    "machine": { "base_input": "/data", "base_output": "/results" },
    "vars": {
        "input": { "file": "merged.root", "type": "kt" },
        "output": { "dir": "run1" }
    }
}
```

Полный пример — `config/kt.json`. Секции:

| Секция | Что настраивает |
|---|---|
| `machine` | корневые пути машины (`base_input` + `vars.input.file`, `base_output` + `vars.output.dir`) |
| `vars` | входной файл и тип анализа (`kt` / `rapidity`) |
| `fit` | `q_max` — диапазон фита по q; `use_default_ip` — дефолтные начальные параметры вместо табличных; `options` — опции ROOT `Fit()` (внимание: `M` принудительно включает TMinuit+MigradImproved, для Minuit2 не указывайте `M`); `minimizer` — `Minuit`/`Minuit2` (по умолчанию `Minuit2`); `retry_with_defaults` — при плохом фите перефит с дефолтными начальными параметрами (лучший по chi2); `use_integral` — интеграл функции по бину вместо центра (опция `I`); `minos_errors` — Minos-ошибки (опция `E`, медленно); `limits` — лимиты параметров (`radius_sq`, `cross`, `lambda` — пары `[min, max]`); `freeze` — зафиксированные параметры: `{"r_os": 0.0, ...}` (значение — фикс; `r_out`/`r_side`/`r_long` задаются в фм, остальные как в параметризации; отсутствующий ключ = свободен). Секция `freeze` при наличии заменяет дефолт целиком: `"freeze": {}` освобождает cross-термы |
| `projections` | ширины срезов: `slice_1d`, `slice_2d` (freeze), `crop_2d`, `slice_ratio`, `fit_over_cf_range`; диапазоны осей 1D-канвасов: `axis_range_1d` (± по q), `cf_y_range_1d`, `fit_over_cf_y_range_1d` (пары `[min, max]`) |
| `binning` | **необязательно**: свои границы `values` (edges) + `names` (на каждый бин, т.е. `values.size() - 1`), опционально `file_names` для имён файлов. Если секции нет — встроенные таблицы kt/rapidity |
| `images` | `need: false` — не сохранять картинки (batch); `format` — `pdf`/`png`/... |
| `stages` | вкл/выкл стадии пайплайна: `cf3d`, `dependency`, `projections_1d`, `projections_2d`, `ratios`. Зависимые стадии требуют, чтобы `cf3d` отрабатывал (или `cf3d.root` уже был создан прошлым прогоном) |
| `selection` | `charges: [0, 1]`, `centralities: [0, 1, 2, 3]` — подмножество для быстрых прогонов |
| `logging` | `level`: `debug`/`info`/`warning`/`error`; `file`: имя лог-файла в выходной директории (пустая строка — только stderr) |
| `threads` | число потоков для параллельного фитирования, `0` — автоматически (`hardware_concurrency`) |

Неизвестные ключи отклоняются с ошибкой, неверные значения валидируются
(например, `binning.names` должно содержать `values.size() - 1` элементов).

При каждом запуске в выходную директорию пишется снапшот эффективного конфига —
`run_config.json` (со всеми дефолтами и полными путями), чтобы всегда было видно,
чем именно получены результаты.

# Тесты

```bash
ctest --test-dir build --output-on-failure
```

# Тулчейн

- **clang-format** — конфиг `.clang-format`, отформатировать всё:

  ```bash
  clang-format -i $(find src tests -name '*.cpp' -o -name '*.h')
  ```

- **clang-tidy** — конфиг `.clang-tidy`, прогнать линтер:

  Linux:

  ```bash
  clang-tidy -p build $(find src -name '*.cpp')
  ```

  macOS (keg-only `llvm` + `-isysroot`, подробнее в AGENTS.md):

  ```bash
  /opt/homebrew/opt/llvm/bin/clang-tidy -p build \
      --extra-arg="-isysroot" --extra-arg="$(xcrun --show-sdk-path)" \
      $(find src -name '*.cpp')
  ```

- **ccache** — подхватывается CMake автоматически, если установлен
- **Sanitizers** — `cmake -S . -B build-asan -DCF_MAKER_SANITIZE=ON`

# Данные для анализа в лабе:

y-bining = 3609 файлов на ядро всего 10 ядер
`../../sthbtmaker-repo/py_output/y_run_3609_fpc/merged_y_run_3609_fpc.root`

kt-bining = 2707 файлов на ядро всего 16 ядер
`../../sthbtmaker-repo/py_output/kt_run_2707_fpc/merged_kt_run_2707_fpc.root`
