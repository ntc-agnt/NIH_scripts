#!/usr/bin/env python3
# Functional-identity comparator for make_heatmap matrices. Compares the numeric
# cells of two output matrices and reports the max absolute and max RELATIVE
# difference. Exit 0 if max relative diff <= tolerance (default 2e-5, i.e. about one
# unit in the 6-significant-figure printed output), else exit 1. Header lines and
# the six leading feature columns are ignored; row sets must match.
#   compare_matrix.py A B [tol]
import sys

def rows(fn):
    out = {}
    for L in open(fn):
        if "\t" not in L:
            continue
        f = L.rstrip("\n").split("\t")
        if f[0] == "Gene ID":            # column header row
            continue
        try:
            out[f[0]] = [float(x) for x in f[6:]]
        except ValueError:
            continue
    return out

a = rows(sys.argv[1]); b = rows(sys.argv[2])
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 2e-5
if set(a) != set(b):
    print("ROWSET DIFFERS (%d vs %d features)" % (len(a), len(b))); sys.exit(2)
maxabs = maxrel = 0.0
for k in a:
    if len(a[k]) != len(b[k]):
        print("COLCOUNT DIFFERS for %s" % k); sys.exit(2)
    for x, y in zip(a[k], b[k]):
        d = abs(x - y); maxabs = max(maxabs, d)
        den = max(abs(x), abs(y))
        if den > 0:
            maxrel = max(maxrel, d / den)
ok = maxrel <= tol
print("max_abs=%.3e max_rel=%.3e tol=%.1e %s" % (maxabs, maxrel, tol, "PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
