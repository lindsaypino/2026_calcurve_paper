"""Regenerate every figure. Scripts that read only committed figures-of-merit run
by default; pass --with-raw to also run the example-curve / model figures that need
the raw MS data pointed at by config.yaml.

    python run_all.py            # FOM-only figures (no raw data needed)
    python run_all.py --with-raw # everything
"""
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))

FOM_ONLY = [
    "figures/fig01b_distributions.py",
    "figures/fig02_triptychs.py",
    "figures/supp_before_after.py",
    "figures/supp_tiers.py",
]
NEEDS_RAW = [
    "figures/fig01b_examples.py",
    "figures/fig03_retentiontime.py",
    "figures/fig_uloq_examples_linear.py",
]

scripts = FOM_ONLY + (NEEDS_RAW if "--with-raw" in sys.argv else [])
os.makedirs(os.path.join(HERE, "output"), exist_ok=True)

failed = []
for s in scripts:
    path = os.path.join(HERE, s)
    if not os.path.exists(path):
        print(f"SKIP  {s} (not present yet)"); continue
    print(f"RUN   {s}")
    r = subprocess.run([sys.executable, path], cwd=HERE)
    if r.returncode != 0:
        failed.append(s)

print("\nDONE." + (f" FAILED: {failed}" if failed else " all figures regenerated."))
sys.exit(1 if failed else 0)
