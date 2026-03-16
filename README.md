# cf_maker_new

Утилита для построения и фита 3D корреляционных функций (HBT) с использованием ROOT.

Проект предназначен для batch-обработки больших ROOT-файлов и сохранения результатов
в виде отдельных выходных файлов.

---
# Сборка

Из корня проекта, 

```bash
mkdir -p build
cd build
cmake ..
make
```

После сборки исполняемый файл `main` будет находиться в директории `build/`.

# Запуск

Программа принимает два аргумента командной строки:

```bash
./main <input.root> <output_dir>
```

Результат работы будет находиться в директории `results/<output_dir>`.

# Данные для анализа в лабе:
При запуске из build:
y-bining = 480 файлов на ядро всего 10 ядер
`../..sthbtmaker-repo/py_output/y_bining_480_fpc/merged_y_bining_480_fpc.root`
kt-bining = 896 файлов на ядро всего 16 ядер
`../..sthbtmaker-repo/py_output/kt_bining_896_fpc/merged_kt_bining_896_fpc.root`