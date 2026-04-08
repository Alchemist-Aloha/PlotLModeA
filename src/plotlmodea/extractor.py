import argparse
import csv
import re
from pathlib import Path


def _find_line_index(lines, needle, start=0):
	needle_lower = needle.lower()
	for idx in range(start, len(lines)):
		if needle_lower in lines[idx].lower():
			return idx
	return -1


def _parse_atom_symbols(lines):
	atom_symbols = {}
	start = _find_line_index(lines, "Cartesian coordinates")
	if start == -1:
		return atom_symbols

	for idx in range(start + 1, len(lines)):
		parts = lines[idx].split()
		if len(parts) >= 7 and parts[0].isdigit():
			atom_symbols[int(parts[0])] = parts[1]
			continue

		# Stop once coordinate table has finished.
		if atom_symbols and lines[idx].strip().startswith("---"):
			break

	return atom_symbols


def _build_atom_symbol_label(name, i, j, k, l, atom_symbols):
	indices = [i, j, k, l]
	tags = []
	for n in indices:
		if n <= 0:
			continue
		sym = atom_symbols.get(n, "A")
		tags.append(f"{sym}{n}")

	if not tags:
		return ""

	if name == "Bond length" and len(tags) >= 2:
		return f"{tags[0]}-{tags[1]}"
	if name == "Bond angle" and len(tags) >= 3:
		return f"{tags[0]}-{tags[1]}-{tags[2]}"
	if name == "Dihedral angle" and len(tags) >= 4:
		return f"{tags[0]}-{tags[1]}-{tags[2]}-{tags[3]}"

	return "-".join(tags)


def _extract_decomposition_matrix(lines):
	start = _find_line_index(lines, "Decomposition of normal modes into local modes")
	if start == -1:
		raise ValueError("Could not find decomposition section for local mode matrix")

	end = _find_line_index(lines, "<<< ACS >>>", start=start)
	if end == -1:
		end = len(lines)

	mode_order = []
	matrix = {}
	current_modes = []

	for idx in range(start, end):
		line = lines[idx].rstrip("\n")
		stripped = line.strip()

		if "Vib. Mode" in line:
			current_modes = [int(x) for x in re.findall(r"\d+", line)]
			for mode in current_modes:
				if mode not in mode_order:
					mode_order.append(mode)
			continue

		if not current_modes:
			continue

		tokens = stripped.split()
		if not tokens:
			continue

		# First data row in each block starts with "q_n:", remaining rows start with q_n index.
		if tokens[0] == "q_n:":
			if len(tokens) < 2:
				continue
			q_n = int(tokens[1])
			values = tokens[2:]
		elif tokens[0].isdigit():
			q_n = int(tokens[0])
			values = tokens[1:]
		else:
			continue

		if len(values) < len(current_modes):
			continue

		row = matrix.setdefault(q_n, {})
		for mode, value in zip(current_modes, values[: len(current_modes)]):
			row[mode] = value

	mode_order = sorted(mode_order)
	q_n_order = sorted(matrix)

	if not mode_order or not q_n_order:
		raise ValueError("Parsed empty decomposition matrix")

	return mode_order, q_n_order, matrix


def _extract_float_values(line):
	return [float(x) for x in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)]


def _extract_normal_mode_properties(lines):
	start = _find_line_index(lines, "Results of vibrations:")
	if start == -1:
		raise ValueError("Could not find normal mode vibration results section")

	end = _find_line_index(lines, "Results of translations and rotations:", start=start)
	if end == -1:
		raise ValueError("Could not find end of normal mode vibration section")

	normal_modes = []
	block_start = None
	block_size = 0

	for idx in range(start, end):
		line = lines[idx].strip()
		if not line:
			continue

		if line.startswith("Irreps"):
			irreps = lines[idx].split()[1:]
			if not irreps:
				continue
			block_start = len(normal_modes)
			block_size = len(irreps)
			for ir in irreps:
				normal_modes.append(
					{
						"Mode": len(normal_modes) + 1,
						"Irrep": ir,
						"Frequency_cm-1": "",
						"ReducedMass_AMU": "",
						"ForceConstant_mDyn_per_A": "",
						"IRIntensity_km_per_mol": "",
					}
				)
			continue

		if block_start is None:
			continue

		if line.startswith("Frequencies"):
			vals = _extract_float_values(lines[idx])
			for j, v in enumerate(vals[:block_size]):
				normal_modes[block_start + j]["Frequency_cm-1"] = v
			continue

		if line.startswith("Reduced masses"):
			vals = _extract_float_values(lines[idx])
			for j, v in enumerate(vals[:block_size]):
				normal_modes[block_start + j]["ReducedMass_AMU"] = v
			continue

		if line.startswith("Force constants"):
			vals = _extract_float_values(lines[idx])
			for j, v in enumerate(vals[:block_size]):
				normal_modes[block_start + j]["ForceConstant_mDyn_per_A"] = v
			continue

		if line.startswith("IR intensities"):
			vals = _extract_float_values(lines[idx])
			for j, v in enumerate(vals[:block_size]):
				normal_modes[block_start + j]["IRIntensity_km_per_mol"] = v

	if not normal_modes:
		raise ValueError("Parsed empty normal mode properties table")

	return normal_modes


def extract_local_mode_tables(input_path, analysis_csv, local_props_csv, normal_modes_csv):
	lines = Path(input_path).read_text(encoding="utf-8", errors="replace").splitlines(True)
	atom_symbols = _parse_atom_symbols(lines)

	analysis_idx = _find_line_index(lines, "Analysis of Local Modes")
	if analysis_idx == -1:
		raise ValueError("Could not find section: Analysis of Local Modes")

	params_table_idx = _find_line_index(lines, "No.  IB(", start=analysis_idx)
	if params_table_idx == -1:
		raise ValueError("Could not find parameter table under Analysis of Local Modes")

	params_rows = []
	# Parse: No.  IB( i j k l ) q_n Name
	sep_below = _find_line_index(lines, "----", start=params_table_idx + 1)
	if sep_below == -1:
		raise ValueError("Malformed Analysis of Local Modes table")

	for idx in range(sep_below + 1, len(lines)):
		raw = lines[idx].rstrip("\n")
		stripped = raw.strip()
		if not stripped:
			continue
		if stripped.startswith("---"):
			break

		parts = raw.split()
		if not parts:
			continue
		if not parts[0].isdigit():
			continue
		if len(parts) < 7:
			continue

		params_rows.append(
			{
				"No": parts[0],
				"i": parts[1],
				"j": parts[2],
				"k": parts[3],
				"l": parts[4],
				"q_n": parts[5],
				"Name": " ".join(parts[6:]),
				"AtomSymbols": _build_atom_symbol_label(
					" ".join(parts[6:]),
					int(parts[1]),
					int(parts[2]),
					int(parts[3]),
					int(parts[4]),
					atom_symbols,
				),
			}
		)

	mode_order, q_n_order, matrix = _extract_decomposition_matrix(lines)
	normal_modes = _extract_normal_mode_properties(lines)

	with Path(analysis_csv).open("w", newline="", encoding="utf-8") as f:
		writer = csv.DictWriter(f, fieldnames=["No", "i", "j", "k", "l", "q_n", "Name", "AtomSymbols"])
		writer.writeheader()
		writer.writerows(params_rows)

	with Path(local_props_csv).open("w", newline="", encoding="utf-8") as f:
		writer = csv.writer(f)
		header = [str(mode) for mode in mode_order]
		writer.writerow(header)
		for q_n in q_n_order:
			writer.writerow([matrix[q_n].get(mode, "") for mode in mode_order])

	with Path(normal_modes_csv).open("w", newline="", encoding="utf-8") as f:
		writer = csv.DictWriter(
			f,
			fieldnames=[
				"Mode",
				"Irrep",
				"Frequency_cm-1",
				"ReducedMass_AMU",
				"ForceConstant_mDyn_per_A",
				"IRIntensity_km_per_mol",
			],
		)
		writer.writeheader()
		writer.writerows(normal_modes)

	return len(params_rows), len(q_n_order), len(mode_order), len(normal_modes)


def main():
	parser = argparse.ArgumentParser(
		description="Extract 'Analysis of Local Modes' and 'Local mode properties' tables from an LModeA .out file."
	)
	parser.add_argument("-i", "--input", help="Path to LModeA output file")
	parser.add_argument(
		"--analysis-csv",
		default="analysis_of_local_modes.csv",
		help="Output CSV path for Analysis of Local Modes table",
	)
	parser.add_argument(
		"--local-props-csv",
		default="local_mode_properties.csv",
		help="Output CSV path for Local mode properties table",
	)
	parser.add_argument(
		"--normal-modes-csv",
		default="normal_mode_properties.csv",
		help="Output CSV path for normal mode properties",
	)

	args = parser.parse_args()

	n_analysis, n_rows, n_cols, n_normal = extract_local_mode_tables(
		args.input,
		args.analysis_csv,
		args.local_props_csv,
		args.normal_modes_csv,
	)

	print(f"Wrote {n_analysis} rows to {args.analysis_csv}")
	print(f"Wrote {n_rows}x{n_cols} matrix to {args.local_props_csv}")
	print(f"Wrote {n_normal} rows to {args.normal_modes_csv}")


if __name__ == "__main__":
	main()

