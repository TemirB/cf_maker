# Текст работы и презентация

Материалы скопированы из `../../tex/analysis-thesis`, версия `b4109e1`.
Исходная директория не изменяется и не нужна для сборки из этой репы.
`report/main.tex` сохранён без изменений, включая исходный титульный лист отчёта
по НИР. Сборка использует отдельный файл `report/compile.tex`.

Состав директории:

- `report/` — текст работы и библиография; рисунки ожидаются в `report/images/`.
- `presentation/` — презентация с обновлениями по проверке и исправлению `cf_maker`;
  рисунки ожидаются в `presentation/fig/`.
- `reference/report-original.pdf` — готовая историческая версия отчёта из исходного репозитория.
- `build/` — PDF и вспомогательные файлы новой сборки; в Git не добавляются.

SHA-256 исходного `report/main.tex`:
`6635f35ffc4ec9d27d9adeae0e170573ebdc6cd8696f30f4625b64278c3936c9`.

## Сборка из корня репозитория

```bash
make docs          # текст работы и презентация
make thesis        # только текст работы
make presentation  # только презентация
make docs-clean    # очистить результаты сборки документов
```

Готовые файлы:

- `docs/thesis/build/report/thesis.pdf`
- `docs/thesis/build/presentation/presentation.pdf`

Эти цели не требуют ROOT, CMake или компиляции программы. Для другого каталога
результатов задайте `DOCS_BUILD_DIR`, для другого latexmk — `LATEXMK`:

```bash
make docs DOCS_BUILD_DIR=/tmp/cf-maker-docs LATEXMK=/path/to/latexmk
```

Можно также собирать внутри директории материалов:

```bash
make -C docs/thesis
make -C docs/thesis report
make -C docs/thesis presentation
make -C docs/thesis doctor
make -C docs/thesis clean
```

Здесь каталог результатов задаётся переменной `BUILD_DIR`, по умолчанию `build`.
Сборка выполняется через `latexmk -xelatex`; latexmk запускает Biber автоматически
для библиографии и повторяет TeX-проходы до разрешения ссылок.
Для строгой проверки используется штатная опция `-usepretex`;
её описание приведено в [руководстве latexmk](https://www.cantab.net/users/johncollins/latexmk/latexmk-488.pdf).

## Зависимости

Нужны XeLaTeX, latexmk и Biber. Минимальный набор пакетов для Ubuntu:

```bash
sudo apt install texlive-xetex texlive-latex-extra texlive-fonts-extra \
  texlive-bibtex-extra texlive-lang-cyrillic latexmk biber fonts-freefont-otf
```

Альтернатива — полный набор TeX Live:

```bash
sudo apt install texlive-full latexmk biber fonts-freefont-otf
```

Документы используют кириллицу, `polyglossia`, `fontspec`, `biblatex-gost`,
`algorithm2e`, а презентация — тему Metropolis. В macOS подходит установленный
MacTeX с `latexmk`, XeLaTeX и Biber в `PATH`.

Эта копия проверена 5 октября 2026 г. на XeLaTeX из TeX Live 2026,
latexmk 4.88 и Biber 2.22: собраны текст работы на 15 страницах и презентация
на 33 слайдах. Новые слайды 17–21 описывают статистические ошибки, точность
гистограмм, проекции, контроль аппроксимации и результаты проверок программы.
Строгий режим отдельно проверен для обоих документов.

## Недостающие рисунки

В исходной копии нет 14 рисунков отчёта и 25 рисунков презентации. При обычной
сборке их места занимают явно обозначенные рамки с именами отсутствующих файлов.
Это позволяет проверять текст и обновлённые слайды; рамки не являются данными
анализа. Такой PDF нужно дополнить исходными рисунками перед защитой.

Скопируйте настоящие файлы в `report/images/` и `presentation/fig/`, сохранив
имена из `\includegraphics`. Для проверки полноты материалов включите строгую
сборку:

```bash
make docs STRICT_ASSETS=1
# или из директории материалов:
make -C docs/thesis STRICT_ASSETS=1
```

Она выполняет новый TeX-проход и завершится ошибкой при отсутствии рисунка,
даже если черновой PDF уже собран.
