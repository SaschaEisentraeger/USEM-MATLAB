#!/bin/bash
# Element-level studies and panel patch tests of the manuscript (MATLAB R2026a).
# Runs six batch sessions side by side, then the timing study alone on an
# otherwise idle machine; afterwards run run_postprocess.sh.
cd "$(dirname "$0")"
mkdir -p logs
run() { matlab -batch "$1" > "logs/$2.log" 2>&1; }
run "ms_element_study(); ms_deviation_check();" element &
run "ms_distortion_sweep(); ms_mesh_figures();" sweep &
run "ms_perturbation(1:12); ms_scaling_split();" perturbation &
run "addpath('../../studies/rajendran_panel/matlab'); rp_run('$(pwd -W)/../data/panel');" panel &
run "unisolvence_check([0 8 8; 15 8 8; 30 8 6; 30 8 7; 30 8 8; 45 8 6; 45 8 8; 60 4 8; 60 8 6; 60 8 7; 60 8 8; 90 1 7; 90 1 8; 90 8 6; 90 8 7; 90 8 8], '$(pwd -W)/../data/unisolvence.txt');" unisolvence &
run "frame_check();" frames &
wait
matlab -singleCompThread -batch "ms_timing();" > logs/timing.log 2>&1
echo ELEMENT STUDIES DONE
