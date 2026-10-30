#!/bin/bash
# Calibration of the Gauss rules of all manuscript studies (MATLAB R2026a).
# Writes manuscript/data/quadrature/<family>.dat; six jobs in parallel.
# ms_setup disables the automatic parallel pool of each session (the code
# base contains parfor loops).
cd "$(dirname "$0")"
mkdir -p logs
run() { matlab -batch "$1" > "logs/$2.log" 2>&1; }
cyl() { for f in "$@"; do echo "ms_quadrature_calibration('$f', 2:8);"; done; }
run "ms_quadrature_calibration('conv010', 2:6); ms_quadrature_calibration('conv020', 2:6); ms_quadrature_calibration('conv030', 2:6);" calib_conv &
run "ms_quadrature_calibration('curved15', 1:10); ms_quadrature_calibration('curved20', 1:10); ms_quadrature_calibration('combined30', 1:10); ms_quadrature_calibration('rect3x2', 12);" calib_3x2 &
run "$(cyl cyl_a000_N1 cyl_a000_N2 cyl_a000_N4 cyl_a000_N8 cyl_a015_N1 cyl_a015_N2 cyl_a015_N4 cyl_a015_N8)" calib_cyl1 &
run "$(cyl cyl_a030_N1 cyl_a030_N2 cyl_a030_N4 cyl_a030_N8)" calib_cyl2 &
run "$(cyl cyl_a045_N1 cyl_a045_N2 cyl_a045_N4 cyl_a045_N8 cyl_a060_N1 cyl_a060_N2)" calib_cyl3 &
run "$(cyl cyl_a060_N4 cyl_a060_N8 cyl_a090_N1 cyl_a090_N2 cyl_a090_N4 cyl_a090_N8)" calib_cyl4 &
wait
# restore the automatic parallel pool (disabled by ms_setup in batch sessions)
matlab -batch "ps = parallel.Settings; ps.Pool.AutoCreate = true;"
echo CALIBRATION DONE
