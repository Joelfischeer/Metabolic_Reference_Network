import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent

commands = [
    # run_network builds the overview and comparison views (both embedded
    # into robust_network_{condition}.html) as its final step — no separate
    # run_comparison call needed for output.
    [sys.executable, "-m", "Edge_cosine_met_reference_network.run_network", "--condition", "healthy"],
    [sys.executable, "-m", "Edge_cosine_met_reference_network.run_network", "--condition", "obese"],
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
