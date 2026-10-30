#!/bin/bash
# Studies of the manuscript with the calibrated Gauss rules (MATLAB R2026a).
# Run run_calibration.sh first; then run_postprocess.sh.
cd "$(dirname "$0")"
mkdir -p logs
run() { matlab -batch "$1" > "logs/$2.log" 2>&1; }
run "ms_convergence(0.1, 2:6, [1 2 4 8 16], 'd010'); ms_convergence(0.2, 2:6, [1 2 4 8 16], 'd020');" conv &
run "ms_convergence(0.3, 2:6, [1 2 4 8 16], 'd030'); ms_spectrum(1:10);" spectrum &
run "ms_growth(2:8); for p = [4 5 6 7 8], ms_transient(p, 'curved20', 50, 20); ms_transient(p, 'combined30', 50, 20); end" dynamics &
run "ms_cylinder([0 90]);" cyl_a000_a090 &
run "ms_cylinder([15 60]);" cyl_a015_a060 &
run "ms_cylinder(30);" cyl_a030 &
run "ms_cylinder(45);" cyl_a045 &
wait
# restore the automatic parallel pool (disabled by ms_setup in batch sessions)
matlab -batch "ps = parallel.Settings; ps.Pool.AutoCreate = true;"
echo RERUNS DONE
