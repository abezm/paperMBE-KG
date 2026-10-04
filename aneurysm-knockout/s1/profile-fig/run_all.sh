#!/bin/bash
# s1/profile-fig: the 12 decks one at a time, 2 threads (as run_untimed.sh). Progress in progress.txt.
cd "$(dirname "$0")"
F="${FEBIO:-febio4}"
echo "start $(date +%H:%M)" > progress.txt
for f in prof-*.feb; do
  OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 "$F" -i "$f" -silent > "${f%.feb}.stdout" 2>&1
  echo "$f exit=$? $(date +%H:%M)" >> progress.txt
done
echo "ALL DONE $(date +%H:%M)" >> progress.txt
