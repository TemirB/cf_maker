#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"
if [[ "${1:-}" == "--help" ]]; then
    echo 'Usage: bash proof/run_statistics_diff.sh [NEW_OUTPUT_DIR] [CONFIG ...]'
    echo 'Defaults: statistics_diff_YYYYmmdd_HHMMSS, config/kt.json config/y.json'
    echo 'Requires the repository, CMake, a C++17 compiler and activated ROOT.'
    exit 0
fi
output_dir="${1:-statistics_diff_$(date +%Y%m%d_%H%M%S)}"
if [[ $# -gt 0 ]]; then shift; fi
if [[ $# -eq 0 ]]; then set -- config/kt.json config/y.json; fi
build_dir="${STATISTICS_BUILD_DIR:-build-statistics-diff}"
cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_dir" --target statistics_diff --parallel "${JOBS:-2}"
"$build_dir/statistics_diff" "$output_dir" "$@"
