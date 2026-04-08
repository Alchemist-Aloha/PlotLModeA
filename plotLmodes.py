import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    tomllib = None


def load_matrix(csv_path):
    with Path(csv_path).open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [row for row in reader]

    mode_labels = [int(x) for x in header]
    data = np.array(rows, dtype=float)
    return mode_labels, data


def normalize_columns_to_percent(data):
    col_sum = data.sum(axis=0, keepdims=True)
    safe = np.where(col_sum == 0.0, 1.0, col_sum)
    return data / safe * 100.0


def parse_mode_list(mode_list_text):
    # Accept comma-delimited positive integers, e.g. "1,2,5,10".
    selected = []
    seen = set()
    for token in mode_list_text.split(","):
        token = token.strip()
        if not token:
            continue
        if not token.isdigit() or int(token) <= 0:
            raise ValueError(f"Invalid normal mode '{token}'. Use a comma-delimited list like 1,2,5,10")
        mode = int(token)
        if mode not in seen:
            selected.append(mode)
            seen.add(mode)
    if not selected:
        raise ValueError("No valid normal modes were provided")
    return selected


def load_analysis_rows(analysis_csv):
    rows = []
    with Path(analysis_csv).open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = (row.get("AtomSymbols") or "").strip()
            if not label:
                label = f"q{row['No']}"
            row["label"] = label
            row["atoms"] = extract_atoms_from_label(label)
            row["mode_size"] = infer_mode_size(row)
            rows.append(row)
    return rows


def infer_mode_size(row):
    fields = [row.get("i", "0"), row.get("j", "0"), row.get("k", "0"), row.get("l", "0")]
    count = 0
    for value in fields:
        try:
            if int(value) != 0:
                count += 1
        except (TypeError, ValueError):
            continue
    return count


def extract_atoms_from_label(label):
    tokens = [tok.strip() for tok in label.split("-") if tok.strip()]
    return set(tokens)


def parse_group_config(config_path):
    config_file = Path(config_path)
    if not config_file.exists():
        return []
    if tomllib is None:
        raise RuntimeError("TOML parsing requires Python 3.11+ (tomllib)")

    with config_file.open("rb") as f:
        raw = tomllib.load(f)

    groups = []
    for group_name, group_block in raw.items():
        entries = group_block if isinstance(group_block, list) else [group_block]
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            elements = entry.get("popElement", [])
            pop_modes = entry.get("popMode", [2, 3, 4])

            normalized_patterns = []
            for pattern in elements:
                if isinstance(pattern, list):
                    atoms = {str(atom).strip() for atom in pattern if str(atom).strip()}
                    if atoms:
                        normalized_patterns.append(atoms)

            normalized_modes = set()
            for mode in pop_modes:
                try:
                    normalized_modes.add(int(mode))
                except (TypeError, ValueError):
                    continue

            if normalized_patterns and normalized_modes:
                groups.append(
                    {
                        "name": str(group_name),
                        "patterns": normalized_patterns,
                        "modes": normalized_modes,
                    }
                )

    return groups


def mode_matches_group(analysis_row, group_rule):
    if analysis_row["mode_size"] not in group_rule["modes"]:
        return False

    atoms = analysis_row["atoms"]
    for pattern in group_rule["patterns"]:
        if pattern.issubset(atoms):
            return True
    return False


def apply_grouping(norm_data, analysis_rows, group_rules):
    n_rows = norm_data.shape[0]
    row_labels = [
        analysis_rows[i]["label"] if i < len(analysis_rows) else f"q_n {i + 1}"
        for i in range(n_rows)
    ]
    row_groups = ["ungrouped"] * n_rows

    if not group_rules:
        return norm_data, row_labels, row_groups

    assigned = set()
    grouped_rows = []
    grouped_labels = []
    grouped_names = []

    for rule in group_rules:
        match_idx = []
        for idx in range(min(n_rows, len(analysis_rows))):
            if idx in assigned:
                continue
            if mode_matches_group(analysis_rows[idx], rule):
                match_idx.append(idx)

        if not match_idx:
            continue

        group_sum = norm_data[match_idx, :].sum(axis=0)
        grouped_rows.append(group_sum)
        grouped_labels.append(rule["name"])
        grouped_names.append(rule["name"])
        assigned.update(match_idx)

    for idx in range(n_rows):
        if idx in assigned:
            continue
        grouped_rows.append(norm_data[idx, :])
        grouped_labels.append(row_labels[idx])
        grouped_names.append("ungrouped")

    if not grouped_rows:
        return norm_data, row_labels, row_groups

    grouped_data = np.vstack(grouped_rows)
    return grouped_data, grouped_labels, grouped_names


def select_rows_for_plot(norm_data, row_labels, row_groups, top_n):
    # Rank by total contribution across all normal modes.
    scores = norm_data.sum(axis=1)
    order = np.argsort(scores)[::-1]
    keep = order[:top_n]

    kept_data = norm_data[keep, :]
    kept_labels = [row_labels[i] if i < len(row_labels) else f"q_n {i + 1}" for i in keep]
    kept_groups = [row_groups[i] if i < len(row_groups) else "ungrouped" for i in keep]

    if top_n < norm_data.shape[0]:
        other = norm_data.sum(axis=0) - kept_data.sum(axis=0)
        kept_data = np.vstack([kept_data, other])
        kept_labels.append("Others")
        kept_groups.append("ungrouped")

    return kept_data, kept_labels, kept_groups


def build_group_colors(row_groups):
    unique_groups = []
    for group in row_groups:
        if group not in unique_groups:
            unique_groups.append(group)

    palette = plt.cm.tab20(np.linspace(0, 1, max(len(unique_groups), 1)))
    colors = {}
    for idx, group in enumerate(unique_groups):
        colors[group] = palette[idx]
    if "ungrouped" in colors:
        colors["ungrouped"] = (0.72, 0.72, 0.72, 1.0)
    return [colors[group] for group in row_groups]


def make_plot(mode_labels, stacked_data, row_labels, row_groups, out_png):
    x = np.arange(len(mode_labels))
    fig, (ax, ax_leg) = plt.subplots(
        1,
        2,
        figsize=(16, 7),
        dpi=150,
        gridspec_kw={"width_ratios": [3, 2]},  # 60% plot, 40% legend
    )

    colors = build_group_colors(row_groups)
    bottom = np.zeros(len(mode_labels))

    for idx, (row, label) in enumerate(zip(stacked_data, row_labels)):
        ax.bar(
            x,
            row,
            bottom=bottom,
            color=colors[idx],
            edgecolor="black",
            linewidth=0.2,
            width=0.9,
            label=label,
        )
        bottom += row

    ax.set_xlim(-0.6, len(mode_labels) - 0.4)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Local Mode Character (%)")
    ax.set_xlabel("Normal Mode u")
    ax.set_xticks(x)
    ax.set_xticklabels(mode_labels, rotation=45, ha="right", fontsize=7)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)

    handles, labels = ax.get_legend_handles_labels()
    ax_leg.axis("off")
    ax_leg.legend(
        handles,
        labels,
        loc="upper left",
        fontsize=7,
        ncol=4,
        frameon=True,
        borderaxespad=0.0,
        handlelength=1.0,
        handletextpad=0.3,
        columnspacing=0.6,
    )

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot local mode character chart from local_mode_properties.csv")
    parser.add_argument("--matrix-csv", default="local_mode_properties.csv")
    parser.add_argument("--analysis-csv", default="analysis_of_local_modes.csv")
    parser.add_argument("--output", default="local_mode_character.png")
    parser.add_argument("--max-modes", type=int, default=48, help="Number of normal modes (columns) to plot")
    parser.add_argument(
        "--mode-list",
        default="",
        help="Comma-delimited normal modes to plot on x-axis, e.g. 1,2,5,10. Overrides --max-modes.",
    )
    parser.add_argument("--top-n", type=int, default=48, help="Number of local modes (rows) to show in legend")
    parser.add_argument(
        "--group-toml",
        default="bodipy.toml",
        help="Optional TOML file with pop/group rules (default: bodipy.toml)",
    )
    args = parser.parse_args()

    mode_labels, data = load_matrix(args.matrix_csv)
    norm_data = normalize_columns_to_percent(data)

    if args.mode_list.strip():
        requested_modes = parse_mode_list(args.mode_list)
        mode_to_index = {mode: idx for idx, mode in enumerate(mode_labels)}
        missing = [m for m in requested_modes if m not in mode_to_index]
        if missing:
            raise ValueError(f"Requested normal modes not found in matrix: {missing}")
        keep_idx = [mode_to_index[m] for m in requested_modes]
        mode_labels = [mode_labels[i] for i in keep_idx]
        norm_data = norm_data[:, keep_idx]
    else:
        max_modes = min(args.max_modes, norm_data.shape[1])
        mode_labels = mode_labels[:max_modes]
        norm_data = norm_data[:, :max_modes]

    analysis_rows = load_analysis_rows(args.analysis_csv)
    group_rules = parse_group_config(args.group_toml)
    grouped_data, grouped_labels, grouped_names = apply_grouping(norm_data, analysis_rows, group_rules)

    stacked_data, row_labels, row_groups = select_rows_for_plot(
        grouped_data,
        grouped_labels,
        grouped_names,
        args.top_n,
    )

    make_plot(mode_labels, stacked_data, row_labels, row_groups, args.output)
    print(f"Saved plot to {args.output}")


if __name__ == "__main__":
    main()
