#!/bin/bash
# Element-level studies and panel patch tests of the manuscript (MATLAB R2026a).
# Runs five batch sessions side by side, then the timing study alone on an
# otherwise idle machine; afterwards run run_postprocess.sh.
cd "$(dirname "$0")"
mkdir -p logs
run() { matlab -batch "$1" > "logs/$2.log" 2>&1; }
run "ms_element_study();" element &
run "ms_distortion_sweep(); ms_mesh_figures();" sweep &
run "ms_perturbation(1:12); ms_scaling_split();" perturbation &
run "addpath('../../studies/rajendran_panel/matlab'); rp_run('$(pwd -W)/../data/panel');" panel &
run "unisolvence_check([0 8 8; 15 8 8; 30 8 6; 30 8 7; 30 8 8; 45 8 6; 45 8 8; 60 4 8; 60 8 6; 60 8 7; 60 8 8; 90 1 7; 90 1 8; 90 8 6; 90 8 7; 90 8 8], '$(pwd -W)/../data/unisolvence.txt');" unisolvence &
wait
matlab -singleCompThread -batch "ms_timing();" > logs/timing.log 2>&1
# restore the automatic parallel pool (disabled by ms_setup in batch sessions)
matlab -batch "ps = parallel.Settings; ps.Pool.AutoCreate = true;"
echo ELEMENT STUDIES DONE
