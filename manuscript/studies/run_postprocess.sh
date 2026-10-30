#!/bin/bash
# Post-processing of all manuscript studies (after run_reruns.sh).
cd "$(dirname "$0")"
python -I postprocess.py
python -I make_companion_data.py
python -I quadrature_table.py
