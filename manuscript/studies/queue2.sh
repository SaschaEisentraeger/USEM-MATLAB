#!/bin/bash
# second queue: remaining studies, started as cores become free
cd "$(dirname "$0")"
while pgrep -f "ms_spectrum" > /dev/null; do sleep 30; done
octave --no-gui --quiet --eval "addpath('../../studies/rajendran_panel/matlab'); rp_run('$(pwd)/../data/panel');" > logs/panel.log 2>&1 &
while pgrep -f "'trapezoid'" > /dev/null; do sleep 30; done
octave --no-gui --quiet --eval "for p = [4 5 6 7 8], ms_transient(p,'curved20',50,20); ms_transient(p,'combined30',50,20); end" > logs/transient.log 2>&1 &
while pgrep -f "'curved','general'" > /dev/null; do sleep 30; done
octave --no-gui --quiet --eval "ms_perturbation(1:12);" > logs/perturbation.log 2>&1
wait
echo QUEUE2 DONE
