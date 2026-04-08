# Plot LModeA Project Guide

This directory contains an installable Python tool for post-processing LModeA output.

Package name:
- plotlmodea

Installed CLI commands:
- lmodea-extract
- lmodea-plot

Compatibility wrappers are still present:
- extractor.py
- plotLmodes.py

Core logic modules live in:
- src/plotlmodea/extractor.py
- src/plotlmodea/plot_lmodes.py

## 1) What each script does

### extractor.py

Purpose:
- Reads an LModeA output text file (for example bodipy.out).
- Extracts three datasets:
  - Analysis of local modes table.
  - Decomposition matrix (normal modes into local modes).
  - Normal mode properties (frequency, reduced mass, force constant, IR intensity).
- Writes those datasets as CSV files used by plotting and downstream analysis.

Main outputs:
- analysis_of_local_modes.csv
- local_mode_properties.csv
- normal_mode_properties.csv

### plotLmodes.py

Purpose:
- Reads local_mode_properties.csv matrix and normalizes local-mode contributions per normal mode to 100%.
- Reads analysis_of_local_modes.csv labels (for legend and grouping metadata).
- Optionally reads a TOML file (default: bodipy.toml) with pop/group rules.
- Produces a stacked bar chart image with:
  - x-axis: selected normal modes.
  - y-axis: local mode character (%).
  - legend on a dedicated right panel (40% of figure width).

Main output:
- local_mode_character.png (or your selected output name)

## 2) Requirements

Python:
- Python 3.11+ recommended.
  - plotLmodes.py uses tomllib for TOML parsing (standard in Python 3.11+).

Packages:
- numpy
- matplotlib

Install and sync with uv:
- uv sync

Run commands with uv:
- uv run lmodea-extract --help
- uv run lmodea-plot --help

## 3) Typical workflow

Step A: Extract CSV tables from LModeA output
- uv run lmodea-extract -i bodipy.out

Step B: Plot contributions
- uv run lmodea-plot

Step C: Plot with grouping TOML
- uv run lmodea-plot --group-toml bodipy.toml --output bodipy_grouped.png

## 4) extractor.py details

### Inputs

Required by content:
- LModeA output file containing these sections:
  - Cartesian coordinates
  - Analysis of Local Modes
  - Decomposition of normal modes into local modes
  - Results of vibrations

CLI option:
- -i, --input
  - Path to input .out file
  - Default: bodipy.out

### Outputs and columns

1) analysis_of_local_modes.csv
- Columns:
  - No, i, j, k, l, q_n, Name, AtomSymbols
- AtomSymbols is generated from the atom index-to-element mapping parsed from Cartesian coordinates.
  - Example for bond: C6-H4
  - Example for angle: C6-C1-C2
  - Example for dihedral: C6-C1-C2-C5

2) local_mode_properties.csv
- Header row is normal mode indices.
- Each following row corresponds to one q_n row from the decomposition matrix.

3) normal_mode_properties.csv
- Columns:
  - Mode
  - Irrep
  - Frequency_cm-1
  - ReducedMass_AMU
  - ForceConstant_mDyn_per_A
  - IRIntensity_km_per_mol

### Full CLI

- uv run lmodea-extract
- uv run lmodea-extract -i bodipy.out
- uv run lmodea-extract --analysis-csv analysis_of_local_modes.csv --local-props-csv local_mode_properties.csv --normal-modes-csv normal_mode_properties.csv

Options:
- --analysis-csv
  - Output path for local mode parameter table
  - Default: analysis_of_local_modes.csv
- --local-props-csv
  - Output path for decomposition matrix
  - Default: local_mode_properties.csv
- --normal-modes-csv
  - Output path for normal mode properties
  - Default: normal_mode_properties.csv

### Error behavior

extractor.py raises descriptive ValueError messages when required sections are missing or malformed, for example:
- Could not find section: Analysis of Local Modes
- Could not find decomposition section for local mode matrix
- Parsed empty decomposition matrix

## 5) plotLmodes.py details

### Inputs

Required files:
- local_mode_properties.csv
- analysis_of_local_modes.csv

Optional file:
- bodipy.toml (or any TOML provided with --group-toml)

### Core plotting logic

1. Load matrix and header mode labels.
2. Normalize each normal-mode column so each column sums to 100%.
3. Select normal modes:
- Either first N modes via --max-modes
- Or explicit list via --mode-list
4. Build local mode labels and metadata from analysis CSV.
5. Apply optional grouping rules from TOML.
6. Rank rows by total contribution and keep top N via --top-n.
7. Plot stacked bars and render legend in right panel.

### Grouping algorithm summary

For each local mode row:
- mode_size is inferred from nonzero i,j,k,l
  - 2 for stretch-like rows
  - 3 for bend-like rows
  - 4 for dihedral-like rows
- AtomSymbols is split by dash into a token set.

For each TOML group rule:
- popMode limits eligible mode_size values.
- popElement provides one or more atom token tuples.
- Match condition:
  - row mode_size is in popMode
  - and at least one popElement tuple is a subset of that row atom token set.

All matched rows for a group are summed into one grouped row.
First matching group wins for a row (rows are assigned once).

### Group colors

- Grouped rows share a color per group name.
- Ungrouped rows use gray.

### Figure layout

- Two-panel figure:
  - Left panel: stacked bar chart (60% width)
  - Right panel: legend area (40% width)
- This prevents chart compression and keeps legend readable.

### Full CLI

- uv run lmodea-plot
- uv run lmodea-plot --output local_mode_character_all.png
- uv run lmodea-plot --max-modes 48 --top-n 48
- uv run lmodea-plot --mode-list 1,2,5,10
- uv run lmodea-plot --group-toml bodipy.toml --output bodipy_grouped.png

Options:
- --matrix-csv
  - Input matrix CSV
  - Default: local_mode_properties.csv
- --analysis-csv
  - Input analysis CSV
  - Default: analysis_of_local_modes.csv
- --output
  - Output image path
  - Default: local_mode_character.png
- --max-modes
  - Number of normal modes from start of matrix header
  - Default: 48
- --mode-list
  - Comma-separated explicit normal mode list
  - Overrides --max-modes
- --top-n
  - Number of local-mode rows (after grouping) to keep in legend/chart
  - Remaining rows merged to Others
  - Default: 48
- --group-toml
  - TOML path for optional grouping rules
  - Default: bodipy.toml

## 6) TOML format reference

Example:

[[ch-methyl]]
popElement = [["C6","H4"],["C6","H6"],["C6","H5"]]
popMode = [2]

[[ch-ring]]
popElement = [["C2","H13"],["C10","H17"]]

Notes:
- Table array name is used as the group name in legend and color mapping.
- popElement entries must match AtomSymbols token style from analysis_of_local_modes.csv.
- If popMode is omitted, default allowed modes are [2, 3, 4].

## 7) Practical tips

- If grouped colors are not visible:
  - Ensure rows are actually grouped (matching token style and popMode).
  - Check that the TOML file path is correct.
- If requested modes fail:
  - Verify mode indices in --mode-list exist in matrix header.
- If no TOML is provided or file is missing:
  - Script runs without grouping.

## 8) Reproducible command set

From this project directory:

1. Extract tables
- uv run lmodea-extract -i bodipy.out

2. Plot baseline
- uv run lmodea-plot --output local_mode_character.png

3. Plot grouped
- uv run lmodea-plot --group-toml bodipy.toml --output bodipy_grouped.png

4. Plot selected normal modes
- uv run lmodea-plot --mode-list 1,2,5,10 --group-toml bodipy.toml --output local_mode_character_selected.png

## 9) File map in this folder

- pyproject.toml
- uv.lock
- src/plotlmodea/extractor.py
- src/plotlmodea/plot_lmodes.py
- extractor.py
- plotLmodes.py
- bodipy.out
- bodipy.toml
- analysis_of_local_modes.csv
- local_mode_properties.csv
- normal_mode_properties.csv
- local_mode_character.png
- bodipy_grouped.png
