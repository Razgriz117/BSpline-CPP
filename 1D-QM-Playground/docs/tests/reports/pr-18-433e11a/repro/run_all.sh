#!/usr/bin/env bash
# Re-run every reproduction from the PR-18 verification report.
# Usage: from 1D-QM-Playground/, after building TISE/build/tise_solver:
#   bash docs/tests/reports/pr-18-433e11a/repro/run_all.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; S="${TISE_SOLVER:-TISE/build/tise_solver}"
for t in control two-sep dup near1e-3 near1e-6 near1e-9 near1e-11 near1e-13 npts1 npts0 nulldelta nulldelta2; do
  echo "=== $t"; "$S" --config "$HERE/$t.yaml" --output-dir "/tmp/pr18-out-$t" 2>&1 | grep -v "below E=0.0" | head -2; echo "exit=${PIPESTATUS[0]}"
  [ -f "/tmp/pr18-out-$t/eigenvalues.dat" ] && sed -n 2,3p "/tmp/pr18-out-$t/eigenvalues.dat"
done
for t in npts0_full npts1_full visnull nulldelta_full; do echo "=== controller.py $t"; python3 controller.py --config "$HERE/$t.yaml" 2>&1 | tail -2; echo "exit=${PIPESTATUS[0]}"; done
echo "=== Eigen / LAPACK standalone"; g++ -O1 -std=c++17 -I/usr/include/eigen3 "$HERE/eig.cpp" -o /tmp/pr18-eig -llapack -lblas && /tmp/pr18-eig
echo "=== LLT on the solver's own overlap.dat"; g++ -O1 -std=c++17 -I/usr/include/eigen3 "$HERE/llt.cpp" -o /tmp/pr18-llt && for t in control near1e-9 near1e-13; do /tmp/pr18-llt /tmp/pr18-out-$t/overlap.dat; done
