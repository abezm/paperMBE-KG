#!/bin/bash
# runs every deck of gate/, z2/, L30/ that has no log yet; 4 threads; wall time in run.out
cd "$(dirname "$0")"
for d in gate z2 L30; do for f in $d/*.feb; do
  [ -e "${f%.feb}.log" ] && continue
  s=$(date +%s); (cd $d && OMP_NUM_THREADS=4 ~/FEBio-stage1/build/bin/febio4 -i $(basename $f) -silent > $(basename ${f%.feb}).stdout 2>&1); rc=$?
  echo "$f rc=$rc $(( $(date +%s)-s ))s $(tail -2 ${f%.feb}.log | head -1 | tr -s ' ')" >> run.out
done; done
echo ALLDONE >> run.out
