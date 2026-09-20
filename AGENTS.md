# Naming

- Пространства имён: `snake_case` (`log`, `draw`, `lcms`, `charge`)
- Типы (классы, структуры, enum): `PascalCase` (`FitResult`, `LCMSAxis`)
- Значения `enum class`: `PascalCase` (`LCMSAxis::Out`, `Binning::Kt`)
- Функции и методы: `snake_case` (`fit_cf_3d`, `bin_center`)
- Переменные: `snake_case` (`fit_results`, `cf_hist`)
- Приватные члены класса: с суффиксом `_` (`pData_`)
- Константы (`constexpr`): `k` + `PascalCase` (`kCount`, `kHc2`)
- Макросы: `SCREAMING_SNAKE`

# Build

```bash
cmake -S . -B build
cmake --build build -j
```

Sanitizers:

```bash
cmake -S . -B build-asan -DCF_MAKER_SANITIZE=ON
cmake --build build-asan -j
```

# Test

```bash
ctest --test-dir build --output-on-failure
```

# Format

```bash
clang-format -i $(find src tests -name '*.cpp' -o -name '*.h')
```

# Lint

Linux (`apt install clang-tidy`):

```bash
clang-tidy -p build $(find src -name '*.cpp')
```

macOS: clang-tidy доступен только в keg-only пакете `llvm` (brew),
и ему нужно указать SDK через `-isysroot`:

```bash
/opt/homebrew/opt/llvm/bin/clang-tidy -p build \
    --extra-arg="-isysroot" --extra-arg="$(xcrun --show-sdk-path)" \
    $(find src -name '*.cpp')
```
