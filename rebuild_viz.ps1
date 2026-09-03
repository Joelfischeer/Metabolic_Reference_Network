Set-Location $PSScriptRoot

uv run -m Edge_metabolism_reference_network.run_network --condition healthy --viz-only
uv run -m Edge_metabolism_reference_network.run_network --condition obese --viz-only
uv run -m Edge_metabolism_reference_network.run_comparison --condition healthy
uv run -m Edge_metabolism_reference_network.run_comparison --condition obese
