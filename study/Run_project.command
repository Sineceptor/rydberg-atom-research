#!/bin/zsh
set -eu
cd -- "${0:A:h}"
task_python="$PWD/.venv/bin/python"
if [[ ! -x "$task_python" ]]; then
  print 'This folder needs its own Python environment. No packages were installed.'
  print 'Follow START_HERE.md: create .venv with Python 3.12 and install requirements.txt.'
  print 'The launcher never uses a Codex cache or another project environment.'
  exit 1
fi
export MPLBACKEND=Agg
"$task_python" run_project.py "$@"
