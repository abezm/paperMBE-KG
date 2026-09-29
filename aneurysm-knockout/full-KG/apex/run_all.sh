#!/bin/bash
# Rate comparison (Figs. 6-7) re-run with the probe at the apex: elem 10 / node 20.
# Copies of the decks in ../ (only the two primary loggers changed), outputs stay in apex/.
cd "$(dirname "$0")"
F="${FEBIO:-febio4}"
run() { d=$(dirname "$1"); b=$(basename "$1"); (cd "$d" && OMP_NUM_THREADS=$2 MKL_NUM_THREADS=$2 "$F" -i "$b" -silent > "${b%.feb}.stdout" 2>&1; echo "$1 exit $?" >> "$OLDPWD/progress.txt"); }
: > progress.txt
run ds-aneu/sac_ds001.feb 4 &
run ds-aneu/fus_ds001.feb 4 &
for f in rate-fusiform-MBE rate-fusiform-full-sload10 rate-fusiform-full-sload100 rate-saccular-MBE rate-saccular-full-sload10 rate-saccular-full-sload100; do run $f.feb 2 & done
wait
echo ALL DONE >> progress.txt
