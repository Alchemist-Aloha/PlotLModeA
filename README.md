# Extract and Plot LModeA Output

This project is a uv-managed, installable Python CLI package for parsing LModeA v3.0 output and plotting local mode contributions.

## Project components

### lmodea-extract:
- Reads an LModeA output text file (for example bodipy.out).
- Extracts three datasets:
  - Analysis of local modes table.
  - Decomposition matrix (normal modes into local modes).
  - Normal mode properties (frequency, reduced mass, force constant, IR intensity).
- Writes those datasets as CSV files used by plotting and downstream analysis.

### lmodea-plot:
- Reads local_mode_properties.csv matrix and normalizes local-mode contributions per normal mode to 100%.
- Reads analysis_of_local_modes.csv labels (for legend and grouping metadata).
- Optionally reads a TOML file with pop/group rules.
- Produces a stacked bar chart image with:
  - x-axis: selected normal modes.
  - y-axis: local mode character (%).
  - legend on a dedicated right panel (40% of figure width).

Main output:
- local_mode_character.png (or your selected output name)

## Setup

Requirements:
- Python >= 3.11
- uv

Install dependencies and create environment:
- `uv sync`

Check command help:
- `uv run lmodea-extract --help`
- `uv run lmodea-plot --help`

## Typical workflow 

Because default filenames are relative, run commands from the dataset folder unless you pass explicit paths.

Option A: run from `bodipy` directory

1. `cd bodipy`
2. `uv run --project .. lmodea-extract -i bodipy.out`
3. `uv run --project .. lmodea-plot --group-toml bodipy.toml --output bodipy_grouped.png`

Option B: run from project root with explicit paths

1. `uv run lmodea-extract -i bodipy/bodipy.out --analysis-csv bodipy/analysis_of_local_modes.csv --local-props-csv bodipy/local_mode_properties.csv --normal-modes-csv bodipy/normal_mode_properties.csv`
2. `uv run lmodea-plot --matrix-csv bodipy/local_mode_properties.csv --analysis-csv bodipy/analysis_of_local_modes.csv --group-toml bodipy/bodipy.toml --output bodipy/bodipy_grouped.png`

## Command reference

### lmodea-extract

Purpose:
- Parse an LModeA `.out` file and export:
  - `analysis_of_local_modes.csv`
  - `local_mode_properties.csv`
  - `normal_mode_properties.csv`

Options:
- `-i`, `--input` (Path to LModeA output file)
- `--analysis-csv` (Output CSV path for Analysis of Local Modes table, default: `analysis_of_local_modes.csv`)
- `--local-props-csv` (Output CSV path for Local mode properties table, default: `local_mode_properties.csv`)
- `--normal-modes-csv` (Output CSV path for normal mode properties, default: `normal_mode_properties.csv`)

### lmodea-plot

Purpose:
- Read extracted CSVs.
- Normalize local-mode contribution per normal mode to 100%.
- Optionally group modes from TOML rules.
- Produce stacked bar plot with legend panel on the right.

#### Inputs

Required files:
- local_mode_properties.csv
- analysis_of_local_modes.csv

Optional file:
- bodipy.toml (or any TOML provided with --group-toml)

#### Core plotting logic

1. Load matrix and header mode labels.
2. Normalize each normal-mode column so each column sums to 100%.
3. Select normal modes:
- Either first N modes via --max-modes
- Or explicit list via --mode-list
4. Build local mode labels and metadata from analysis CSV.
5. Apply optional grouping rules from TOML.
6. Rank rows by total contribution and keep top N via --top-n.
7. Plot stacked bars and render legend in right panel.

#### Grouping algorithm summary

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

#### Group colors

- Grouped rows share a color per group name.
- Ungrouped rows use gray.

#### Figure layout

- Two-panel figure:
  - Left panel: stacked bar chart (60% width)
  - Right panel: legend area (40% width)
- This prevents chart compression and keeps legend readable.

Options:
- `--matrix-csv` (CSV file with local mode properties matrix (q_n vs normal modes), default: `local_mode_properties.csv`)
- `--analysis-csv` (CSV file with local mode property analysis (No, i, j, k, l, q_n, Name, AtomSymbols), default: `analysis_of_local_modes.csv`)
- `--output` (Output PNG file for the plot, default: `local_mode_character.png`)
- `--max-modes` (Number of normal modes (columns) to plot, default: `888`)
- `--mode-list` (comma-delimited mode ids, overrides `--max-modes`)
- `--top-n` (Number of local modes (rows) to show in legend, default: `888`)
- `--group-toml` (Optional TOML file with pop/group rules)

## TOML grouping format

Example:

```toml
[[ch-methyl]]
popElement = [["C6","H4"],["C6","H6"],["C6","H5"]]
popMode = [2]

[[ch-ring]]
popElement = [["C2","H13"],["C10","H17"]]
```

Rules:
- Group name comes from table array name (for example `ch-methyl`).
- `popElement` tuple tokens must match `AtomSymbols` tokens in analysis CSV.
- `popMode` limits matching mode sizes (`2`, `3`, `4`).
- If `popMode` is omitted, defaults are `[2, 3, 4]`.
## Project structure

```
project/
  pyproject.toml
  uv.lock
  README.md
  src/
    plotlmodea/
      __init__.py
      __main__.py
      extractor.py
      plot_lmodes.py
  example/
    example.out
    example.toml
    analysis_of_local_modes.csv
    local_mode_properties.csv
    normal_mode_properties.csv
```

Notes:
- The installable package code lives under `src/plotlmodea`.
## Source files

- `src/plotlmodea/extractor.py`: extraction engine and CLI main
- `src/plotlmodea/plot_lmodes.py`: plotting engine and CLI main
- `src/plotlmodea/__main__.py`: module execution entry (`python -m plotlmodea`)
- `pyproject.toml`: package metadata, dependencies, and scripts
- `uv.lock`: locked dependency set
