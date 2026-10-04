#!/bin/bash
# alp2/ (3 Oct 2026): rate comparison with the adopted severities (fus cc_in 0.95, sac cc_in 0.9).
# 1. search run for the stopping time (untimed), then make_stop.py writes the timed decks;
# 2. untimed runs: time-step study (5), rate comparison (6), the two logs of Table 3;
# 3. timed runs, one thread, alone: stop/timing_stop.py and timing_fair2.py.
# Executable: copy ~/FEBio-stage1. Progress in progress.txt.
cd "$(dirname "$0")"
F="${FEBIO:-febio4}"
run() { (cd "$1" && OMP_NUM_THREADS=$3 MKL_NUM_THREADS=$3 "$F" -i "$2" -silent > "${2%.feb}.stdout" 2>&1; echo "$1/$2 exit=$? $(date +%H:%M)" >> "$OLDPWD/progress.txt"); }
echo "start $(date +%H:%M)" > progress.txt
run speedup-fair/stop/find fus-sl100-2st-ds5-500.feb 4
python3 make_stop.py >> progress.txt 2>&1 || { echo "make_stop failed" >> progress.txt; exit 1; }
D=full-KG/apex/ds-aneu-t1
run $D fus_ds001.feb 4 &
run $D sac_ds001.feb 4 &
run $D fus_ds005.feb 4 & p1=$!
run $D fus_ds01.feb 4 & p2=$!
run $D sac_ds01.feb 4 & p3=$!
wait $p1 $p2 $p3
for f in rate-fusiform-MBE rate-saccular-MBE rate-saccular-full-sload10; do run full-KG/apex $f.feb 2 & done
for f in rate-fusiform-full-sload10 rate-fusiform-full-sload100 rate-saccular-full-sload100; do run full-KG/apex/to5 $f.feb 2 & done
run timing a_ds1.0.feb 1 &
run timing conv-aneu-MBE.feb 1 &
wait
echo "UNTIMED DONE $(date +%H:%M)" >> progress.txt
pgrep -x febio4 > /dev/null && { echo "another febio4 is running: timing refused" >> progress.txt; exit 1; }
(cd speedup-fair/stop && python3 timing_stop.py > timing_stop.stdout 2>&1; echo "timing_stop exit=$? $(date +%H:%M)" >> ../../progress.txt)
(cd speedup-fair && python3 timing_fair2.py > timing_fair2.stdout 2>&1; echo "timing_fair2 exit=$? $(date +%H:%M)" >> ../progress.txt)
echo ALL DONE >> progress.txt
