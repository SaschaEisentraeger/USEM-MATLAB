#!/bin/bash
# Post-processing of all manuscript studies (after run_reruns.sh).
cd "$(dirname "$0")"
CYLRES=../../studies/thick_cylinder/results
CYL_RAW=$CYLRES/raw_calibrated CYL_TAB=$CYLRES/tables_calibrated \
  CYL_GEN=$CYLRES/generated_calibrated python -I ../../studies/thick_cylinder/postprocess.py
python -I postprocess.py
python -I make_companion_data.py
python -I quadrature_table.py
