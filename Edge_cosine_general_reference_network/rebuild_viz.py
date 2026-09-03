import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent

commands = [
    # run_network --viz-only rebuilds the overview and comparison views (both
    # embedded into robust_network_{condition}.html) from cache — no separate
    # run_comparison call needed for output; run it manually if you just want
    # the shared/only-in-network/only-in-reference stats printed to console.
    [sys.executable, "-m", "Edge_cosine_general_reference_network.run_network", "--condition", "healthy", "--viz-only"],
    [sys.executable, "-m", "Edge_cosine_general_reference_network.run_network", "--condition", "obese",   "--viz-only"],
]

for cmd in commands:
    print(f"\n{'='*60}")
    print(f"Running: {' '.join(cmd[1:])}")
    print('='*60)
    result = subprocess.run(cmd, cwd=ROOT)
    if result.returncode != 0:
        print(f"\n[!] Command failed with exit code {result.returncode}. Stopping.")
        sys.exit(result.returncode)

print("\n[ok] All done.")
