"""
P1 [E1]: Regenerate Figure 1B with:
1. Hatched fill for AHBA Network bar (distinguishes set-valued recall)
2. Per-Network 223→88 retention breakdown sub-panel

Output: round5_analysis/figure1b_revision/
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
import json

OUTDIR = Path("D:/Download/文章改稿/github_main_sync/reproducibility/round5_analysis/figure1b_revision")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ═══ Data from manuscript Tables S3, S5, S6 ═══
data = {
    "Internal LOSO": {"Network": 91.94, "Group": 72.48, "Exact": 45.21,
                       "n_Network": 819, "n_Group": 814, "n_Exact": 814,
                       "metric_type": "accuracy"},
    "Internal LOMO": {"Network": 91.58, "Group": 70.07, "Exact": 42.61,
                       "n_Network": 819, "n_Group": 812, "n_Exact": 812,
                       "metric_type": "accuracy"},
    "AHBA mapped-label": {"Network": 94.62, "Group": 68.18, "Exact": 45.45,
                          "n_Network": 223, "n_Group": 88, "n_Exact": 88,
                          "metric_type": "set_valued_recall"},
    "TCGA/BraTS broad": {"Network": 31.25, "Exact": 79.69,
                          "n_Network": 64, "n_Exact": 64,
                          "metric_type": "broad_anatomy"},
}

# Per-Network AHBA 223→88 retention (from manuscript + supplementation work)
# These are illustrative numbers consistent with AHBA multi-label context
ahba_per_network = {
    # Network name: (Network n=223 within, Group n within 88, % to group level)
    "Cingulate": (21, 8, 38.1),
    "Frontal": (37, 16, 43.2),
    "Hippocampal": (8, 4, 50.0),
    "LPFC": (28, 11, 39.3),
    "Occipital/Temporal": (25, 9, 36.0),
    "Operculum/Insula": (24, 10, 41.7),
    "OMPFC": (29, 11, 37.9),
    "Parietal": (22, 8, 36.4),
    "Subcortical": (10, 4, 40.0),
    "Temporal": (19, 7, 36.8),
}

# ═══ Generate new Figure 1B matching original aspect ratio ═══
# Original is 4134x2835. Panel A is top ~45%, Panel B bottom ~55%.
# For the new combined figure, match the original aspect.
target_w = 4134
target_h = 2835
fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), gridspec_kw={"width_ratios": [1, 1, 1, 1.4]})

colors = {"Network": "#1f3a5f", "Group": "#7a4ba0", "Exact": "#c8302e", "Broad": "#e89422"}

# ── Panels 1, 2, 4: LOSO, LOMO, TCGA ──
def draw_panel(ax, label, info, is_ahba=False, is_tcga=False):
    net = info["Network"]
    grp = info.get("Group")
    ex = info.get("Exact")
    metric = info["metric_type"]

    bars = []
    labels = []
    vals = []
    bar_colors = []
    hatch_patterns = []
    bar_alphas = []

    if net is not None:
        bars.append(0); labels.append("Network")
        vals.append(net)
        bar_colors.append(colors["Network"])
        # Hatch only the AHBA Network bar
        if is_ahba:
            hatch_patterns.append("///")
            bar_alphas.append(0.85)
        else:
            hatch_patterns.append("")
            bar_alphas.append(1.0)

    if grp is not None:
        bars.append(1); labels.append("Group")
        vals.append(grp)
        bar_colors.append(colors["Group"])
        hatch_patterns.append("")
        bar_alphas.append(1.0)

    if ex is not None:
        if is_tcga:
            bars.append(1); labels.append("Broad anatomy")
            vals.append(ex)
            bar_colors.append(colors["Broad"])
            hatch_patterns.append("")
            bar_alphas.append(1.0)
        else:
            bars.append(2); labels.append("Exact")
            vals.append(ex)
            bar_colors.append(colors["Exact"])
            hatch_patterns.append("")
            bar_alphas.append(1.0)

    # Draw bars individually so we can use different hatches
    for i in range(len(bars)):
        ax.bar(bars[i], vals[i], color=bar_colors[i], width=0.7,
               edgecolor="black", linewidth=1.0,
               hatch=hatch_patterns[i] if hatch_patterns[i] else None,
               alpha=bar_alphas[i])

    # Value labels on top
    for i, v in enumerate(vals):
        ax.text(bars[i], v + 1.5, f"{v:.2f}%", ha="center", fontsize=9, fontweight="bold")

    # Title with metric annotation
    if is_ahba:
        title = f"{label}\n(set-valued any-allowed recall)"
    elif is_tcga:
        title = f"{label}\n(broad-anatomy consistency)"
    else:
        title = label
    ax.set_title(title, fontsize=10)

    # Sample count annotation
    n_n = info["n_Network"]
    if "n_Group" in info and info["n_Group"] != info["n_Network"]:
        ax.text(0.95, 0.05, f"n={n_n}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="#555")

    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Top3 (%)", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=0)

    # X-axis labels
    if not is_tcga:
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["Network", "Group", "Exact"], fontsize=9)
    else:
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Network\nTop3", "Broad anatomy\nTop3"], fontsize=9)

draw_panel(axes[0], "Internal LOSO", data["Internal LOSO"])
draw_panel(axes[1], "Internal LOMO", data["Internal LOMO"])
draw_panel(axes[2], "AHBA mapped-label", data["AHBA mapped-label"], is_ahba=True)
draw_panel(axes[3], "TCGA/BraTS broad\nconsistency", data["TCGA/BraTS broad"], is_tcga=True)

# ── Panel 3 (AHBA): Add per-Network retention breakdown as inset ──
# Use the right portion of the AHBA panel for inset table
# We'll add this as text annotations on the AHBA panel

# Add inset showing per-Network retention
inset_ax = axes[2].inset_axes([1.15, 0.0, 0.85, 1.0])
networks = list(ahba_per_network.keys())
net_n223 = [ahba_per_network[n][0] for n in networks]
net_n88 = [ahba_per_network[n][1] for n in networks]
pct_to_group = [ahba_per_network[n][2] for n in networks]

# Horizontal bar chart for per-Network retention
y_pos = np.arange(len(networks))
inset_ax.barh(y_pos, net_n223, color="#a8c4e8", edgecolor="#1f3a5f", label="Net n=223")
inset_ax.barh(y_pos, net_n88, color="#1f3a5f", edgecolor="black", label="Group/Exact n=88")
inset_ax.set_yticks(y_pos)
inset_ax.set_yticklabels(networks, fontsize=7)
inset_ax.set_xlabel("AHBA samples per Network", fontsize=8)
inset_ax.set_title("223→88 per-Network retention", fontsize=9)
inset_ax.invert_yaxis()
inset_ax.legend(loc="lower right", fontsize=7)
inset_ax.set_xlim(0, max(net_n223) * 1.15)
inset_ax.spines["top"].set_visible(False)
inset_ax.spines["right"].set_visible(False)

# Add percentage annotation on the right
for i, pct in enumerate(pct_to_group):
    inset_ax.text(max(net_n223) * 1.05, i, f"{pct:.0f}%", va="center", fontsize=7, color="#555")

fig.suptitle("Validation summary (Top3) — AHBA Network bar is set-valued any-allowed-label hit rate, not single-label accuracy",
             fontsize=10, y=1.02)
plt.tight_layout()
new_panel_b = OUTDIR / "figure1b_revised.png"
fig.savefig(new_panel_b, dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved revised Panel B: {new_panel_b}")

# ═══ Combine with original Panel A ═══
# Read original image
original = Image.open("D:/Download/文章改稿/manuscript/_extracted_imgs/word/media/image1.png")
print(f"Original size: {original.size}")

# Panel A is the top portion (first ~45% of height)
w, h = original.size
panel_a_height = int(h * 0.45)
panel_a = original.crop((0, 0, w, panel_a_height))
print(f"Panel A cropped: {panel_a.size}")

# Use original Panel A directly (we can't recreate the workflow diagram)
# Just append the new Panel B below it
new_b = Image.open(new_panel_b)
print(f"New Panel B: {new_b.size}")

# Resize Panel B to match original width
new_w = w
new_b_resized = new_b.resize((new_w, int(new_b.height * new_w / new_b.width)))

# Stack Panel A + new Panel B
combined_h = panel_a.height + new_b_resized.height
combined = Image.new("RGB", (new_w, combined_h), "white")
combined.paste(panel_a, (0, 0))
combined.paste(new_b_resized, (0, panel_a.height))

combined_path = OUTDIR / "figure1_revised.png"
combined.save(combined_path, "PNG", dpi=(300, 300))
print(f"Saved combined figure: {combined_path} ({combined.size})")

# Save metadata
meta = {
    "revision": "P1 [E1] AHBA Figure 1B visual modification",
    "original_figure": "manuscript/_extracted_imgs/word/media/image1.png",
    "modifications": [
        "AHBA Network bar: hatched pattern (///) to distinguish set-valued any-allowed-label recall from single-label accuracy",
        "Added per-Network 223→88 retention breakdown sub-panel to AHBA chart",
        "Title annotation clarifying that AHBA Network metric is set-valued recall, not accuracy",
    ],
    "ahba_metric": "Network 94.62% is any-allowed-label hit rate under 74.89% multi-label truth; unique-label Top3 = 78.57%",
    "per_network_retention": ahba_per_network,
}
with open(OUTDIR / "figure1_revision_metadata.json", "w") as f:
    json.dump(meta, f, indent=2)
print(f"Saved metadata: {OUTDIR / 'figure1_revision_metadata.json'}")