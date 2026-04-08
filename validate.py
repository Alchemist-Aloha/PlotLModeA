import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


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


def load_qn_labels(analysis_csv):
    qn_labels = []

    with Path(analysis_csv).open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = (row.get("AtomSymbols") or "").strip()
            if not label:
                label = f"q{row['No']}"

            qn_labels.append(label)

    return qn_labels


def select_rows_for_plot(norm_data, qn_labels, top_n):
    # Rank by total contribution across all normal modes.
    scores = norm_data.sum(axis=1)
    order = np.argsort(scores)[::-1]
    keep = order[:top_n]

    kept_data = norm_data[keep, :]
    kept_labels = [qn_labels[i] if i < len(qn_labels) else f"q_n {i + 1}" for i in keep]

    if top_n < norm_data.shape[0]:
        other = norm_data.sum(axis=0) - kept_data.sum(axis=0)
        kept_data = np.vstack([kept_data, other])
        kept_labels.append("Others")

    return kept_data, kept_labels


def make_plot(mode_labels, stacked_data, row_labels, out_png):
    x = np.arange(len(mode_labels))
    fig, ax = plt.subplots(figsize=(16, 7), dpi=150)

    colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(row_labels)))
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

    # Put legend to the right like the reference image.
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        fontsize=7,
        ncol=2,
        frameon=True,
        borderaxespad=0.0,
        handlelength=1.0,
        handletextpad=0.3,
        columnspacing=0.6,
    )

    fig.tight_layout(rect=[0.0, 0.0, 0.82, 1.0])
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

    qn_labels = load_qn_labels(args.analysis_csv)
    stacked_data, row_labels = select_rows_for_plot(norm_data, qn_labels, args.top_n)

    make_plot(mode_labels, stacked_data, row_labels, args.output)
    print(f"Saved plot to {args.output}")


if __name__ == "__main__":
    main()
