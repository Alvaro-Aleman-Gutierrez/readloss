"""Figures for deciding what to trust (``readloss plot``).

Each figure answers one question a user has to decide on, and each writes, next to the
image, a TSV with exactly the values drawn, so the figure can be redrawn elsewhere.

  selftest   Before a screen: can the rule see each gene? Where do its reads go?
  power      After a screen: are absences and intermediates real, or short of reads?
  coverage   One run: do the candidate reads cover the whole protein, or one region?
  margins    One run: how many reads sit near the limits of the assignment rule?
  identity   One run: are candidate reads as similar to the queries as control reads?

Colours follow the Okabe-Ito palette, which stays distinct for colour-blind readers.
"""
import csv
import math
import random

COL = {"present": "#0072B2", "absent": "#D55E00", "intermediate": "#E69F00", "uninformative": "#9E9E9E",
       "kept": "#0072B2", "discarded": "#D55E00", "identity_only": "#CC79A7", "control": "#009E73",
       "control2": "#56B4E9", "other": "#999999", "nohit": "#DDDDDD"}
CONTROL_COLS = ["#009E73", "#56B4E9", "#F0E442", "#999999"]


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left"})
    return plt


def _write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _save(fig, out):
    fig.savefig(out, bbox_inches="tight", dpi=200)
    data = out.rsplit(".", 1)[0] + ".tsv"
    return data


# --- selftest -------------------------------------------------------------------------------

SELF_PARTS = [("assigned_self", "assigned to its gene", "kept"),
              ("assigned_candidate", "to the candidate", "discarded"),
              ("assigned_other_controls", "to another control", "control2"),
              ("assigned_other_genes", "to another family gene", "other"),
              ("discarded", "discarded by the rule", "identity_only"),
              ("no_hit", "no hit", "nohit")]


def selftest(rows, out, min_assignable=0.2):
    """Stacked bar per gene: where its simulated reads went under the rule."""
    plt = _plt()
    rows = [r for r in rows if r.get("gene")]
    fig, ax = plt.subplots(figsize=(6.5, 0.55 * len(rows) + 1.4))
    data = []
    for i, r in enumerate(reversed(rows)):
        n = _num(r["fragments"]) or 1
        left = 0.0
        for key, label, col in SELF_PARTS:
            v = _num(r.get(key, 0)) / n
            if v > 0:
                ax.barh(i, v, left=left, color=COL[col], edgecolor="white", linewidth=0.5,
                        label=label if i == 0 else None)
            left += v
            data.append({"gene": r["gene"], "role": r["role"], "part": key, "fraction": round(v, 4)})
        ax.text(1.01, i, f"{_num(r['assignable_fraction']):.2f}", va="center", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{r['gene']} ({r['role']})" for r in reversed(rows)])
    ax.axvline(min_assignable, color="#444444", lw=0.8, ls="--")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of simulated reads of the gene")
    ax.set_title("Can the rule assign each gene's reads to it?")
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.28))
    fig.text(0.01, 0.01, f"Number at right: assignable fraction. Dashed line: {min_assignable:g}; a candidate "
             "below it cannot be tested (calls are uninformative).", fontsize=7, color="#555555")
    return _write_tsv(_save(fig, out), data, ["gene", "role", "part", "fraction"])


# --- power ----------------------------------------------------------------------------------

def power(rows, out, absent_max=0.10, present_min=0.30, min_expected=10.0, label_max=30):
    """Expected reads against ratio: calls made with too few reads stand out on the left."""
    plt = _plt()
    pts, skipped = [], 0
    for r in rows:
        e, v = _num(r.get("expected")), _num(r.get("ratio"))
        if not (e == e and e > 0 and v == v):
            skipped += 1
            continue
        hi = _num(r.get("ci_high_cond") or r.get("ci_high"))
        lo = _num(r.get("ci_low_cond") or r.get("ci_low"))
        pts.append({"label": r.get("label") or r.get("run", ""), "expected": e, "ratio": v,
                    "ci_low": lo, "ci_high": hi, "status": r.get("status", ""), "flags": r.get("flags", "")})
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.axhspan(absent_max, present_min, color="#EEEEEE", lw=0, zorder=0)
    ax.axvspan(0.3, min_expected, color="#F6F6F6", lw=0, zorder=0)
    ax.axvline(min_expected, color="#666666", lw=0.8, ls="--")
    ax.text(min_expected, 0.02, f" {min_expected:g} expected reads", transform=ax.get_xaxis_transform(),
            fontsize=7, color="#555555", va="bottom")
    seen = set()
    for p in pts:
        c = COL.get(p["status"], "#000000")
        if p["ci_low"] == p["ci_low"] and p["ci_high"] == p["ci_high"] and p["ci_high"] != math.inf:
            ax.plot([p["expected"]] * 2, [p["ci_low"], p["ci_high"]], color=c, lw=0.8, alpha=0.5)
        ax.plot(p["expected"], p["ratio"], "o", ms=5, color=c, mec="white", mew=0.5,
                label=p["status"] if p["status"] not in seen else None)
        seen.add(p["status"])
        if len(pts) <= label_max:
            ax.annotate(p["label"], (p["expected"], p["ratio"]), xytext=(4, 3), textcoords="offset points",
                        fontsize=6.5, color="#333333")
    ax.set_xscale("log")
    xs = [p["expected"] for p in pts] or [1, 100]
    ax.set_xlim(max(0.3, min(min(xs), min_expected) / 2.5), max(max(xs), min_expected) * 2.5)
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    lo_, hi_ = ax.get_xlim()
    ticks = [t for t in (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000) if lo_ <= t <= hi_]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ys = [p["ci_high"] for p in pts if p["ci_high"] == p["ci_high"] and p["ci_high"] != math.inf] + \
        [p["ratio"] for p in pts] + [1.2]
    ax.set_ylim(-0.05, min(3, max(ys) * 1.08))
    ax.axhline(1, color="#AAAAAA", lw=0.6, ls=":")
    ax.set_xlabel("Expected candidate reads (log scale)")
    ax.set_ylabel("Observed / expected")
    ax.set_title("Which calls rest on enough reads?")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    note = "Left of the dashed line: download more spots before trusting the call."
    if skipped:
        note += f"  {skipped} run(s) without expected reads not shown."
    fig.text(0.01, -0.02, note, fontsize=7, color="#555555")
    return _write_tsv(_save(fig, out), pts, ["label", "expected", "ratio", "ci_low", "ci_high", "status", "flags"])


# --- per-run figures ------------------------------------------------------------------------

def _depth(spans, length, step=1):
    L = int(length)
    d = [0] * (L + 2)
    for a, b in spans:
        a, b = max(1, int(a)), min(L, int(b))
        if b >= a:
            d[a] += 1
            d[b + 1] -= 1
    out, run = [], 0
    for i in range(1, L + 1):
        run += d[i]
        out.append(run)
    return out


def coverage(counts, family, res, out, label=""):
    """Reads covering each position of the candidate's and the controls' reference proteins."""
    plt = _plt()
    genes = [family.candidate, *family.controls]
    fig, axes = plt.subplots(len(genes), 1, figsize=(6.5, 1.35 * len(genes) + 0.9), squeeze=False)
    data = []
    exp = _num(res.get("expected"))
    for k, (g, ax) in enumerate(zip(genes, axes[:, 0])):
        L = int(family.lengths[g])
        kept = _depth(counts.spans.get(g, []), L)
        x = list(range(1, L + 1))
        if g == family.candidate:
            disc = _depth(counts.discarded_spans, L)
            io = _depth(counts.identity_only_spans, L)
            ax.fill_between(x, disc, step="mid", color=COL["discarded"], alpha=0.25, lw=0,
                            label="discarded by the rule")
            ax.fill_between(x, io, step="mid", color=COL["identity_only"], alpha=0.55, lw=0,
                            label="discarded, low identity only")
            ax.fill_between(x, kept, step="mid", color=COL["kept"], alpha=0.85, lw=0, label="counted")
            spans = counts.spans.get(g, [])
            mean_len = sum(b - a + 1 for a, b in spans) / len(spans) if spans else 0
            if exp == exp and exp > 0:
                ml = mean_len or (sum(b - a + 1 for a, b in counts.discarded_spans) /
                                  max(1, len(counts.discarded_spans))) or 40
                ax.axhline(exp * ml / L, color="#333333", lw=0.8, ls="--", label="expected if present")
            for i in range(L):
                data.append({"gene": g, "position": i + 1, "counted": kept[i], "discarded": disc[i],
                             "identity_only": io[i]})
            ax.set_title(f"{g} (candidate) · {counts.by_gene.get(g, 0)} reads counted, "
                         f"{counts.discarded} discarded" + (f" · {label}" if label else ""), fontsize=9)
            ax.legend(frameon=False, fontsize=6.5, ncol=4, loc="upper right")
        else:
            col = CONTROL_COLS[(k - 1) % len(CONTROL_COLS)]
            ax.fill_between(x, kept, step="mid", color=col, alpha=0.85, lw=0)
            for i in range(L):
                data.append({"gene": g, "position": i + 1, "counted": kept[i], "discarded": "", "identity_only": ""})
            ax.set_title(f"{g} (control) · {counts.by_gene.get(g, 0)} reads", fontsize=9)
        ax.set_xlim(1, L)
        ax.set_ylabel("Reads")
        top = max([1] + [v for line in ax.get_lines() for v in line.get_ydata()] +
                  [max(kept or [0])] + ([max(_depth(counts.discarded_spans, L) or [0])] if g == family.candidate else []))
        ax.set_ylim(0, top * (1.75 if g == family.candidate else 1.15))
    axes[-1, 0].set_xlabel("Position on the reference protein of each gene (aa)")
    fig.suptitle("Do the candidate reads cover the whole protein?", x=0.01, ha="left", fontsize=10,
                 fontweight="bold")
    fig.tight_layout()
    return _write_tsv(_save(fig, out), data, ["gene", "position", "counted", "discarded", "identity_only"])


def margins(counts, family, rule, out, label="", clip=20.0):
    """Bit margin over the closest competitor against identity, with the limits of the rule."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.axvspan(-1, rule.margin, color="#F3E3DA", lw=0, zorder=0)
    ax.axhspan(0, rule.min_identity, color="#F3E3DA", lw=0, zorder=0)
    ax.axvline(rule.margin, color="#666666", lw=0.8, ls="--")
    ax.axhline(rule.min_identity, color="#666666", lw=0.8, ls="--")
    rng = random.Random(1)
    data = []
    groups = {}
    for g, m, pid, kept in counts.margins:
        role = "candidate" if g == family.candidate else "control"
        key = f"{role}, {'counted' if kept else 'discarded'}"
        x = clip + 1.5 + rng.uniform(-1.2, 1.2) if m >= clip else m + rng.uniform(-0.2, 0.2)
        groups.setdefault(key, []).append((x, pid))
        data.append({"gene": g, "role": role, "margin": round(m, 2), "identity": pid, "counted": kept})
    style = {"candidate, counted": (COL["kept"], "o", 18), "candidate, discarded": (COL["discarded"], "X", 26),
             "control, counted": (COL["control"], "o", 8), "control, discarded": (COL["identity_only"], "X", 18)}
    for key in ["control, counted", "control, discarded", "candidate, discarded", "candidate, counted"]:
        if key in groups:
            xs, ys = zip(*groups[key])
            c, mk, s = style[key]
            ax.scatter(xs, ys, s=s, c=c, marker=mk, alpha=0.7 if key.startswith("control") else 0.9,
                       linewidths=0.3, edgecolors="white", label=f"{key} ({len(xs)})")
    ax.axvspan(clip + 0.2, clip + 2.9, color="#F2F2F2", lw=0, zorder=0)
    ax.text(clip + 1.5, 1.0, "clear\nwinner", transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=7, color="#555555")
    ax.set_xlim(-1, clip + 3)
    ax.set_xticks([t for t in range(0, int(clip) + 1, 5)] + [clip + 1.5])
    ax.set_xticklabels([str(t) for t in range(0, int(clip) + 1, 5)] + [f"\u2265{clip:g}"])
    ax.set_ylim(max(0, min([p for _, _, p, _ in counts.margins] + [rule.min_identity]) - 5), 101)
    ax.set_xlabel("Margin over the closest competitor (bits)")
    ax.set_ylabel("Identity to the query (%)")
    ax.set_title("How many reads sit near the limits of the rule?" + (f"  {label}" if label else ""))
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.text(0.01, -0.09, f"Shaded: rejected for a candidate read (margin < {rule.margin:g} or identity < "
             f"{rule.min_identity:g}%). A cloud pressed against margin 0 means near-identical family members.\n"
             "Control reads of low identity may come from genes missing from the family (e.g. other kinases).",
             fontsize=7, color="#555555")
    return _write_tsv(_save(fig, out), data, ["gene", "role", "margin", "identity", "counted"])


def identity(counts, family, out, label=""):
    """Identity of counted candidate reads, discarded candidate reads and control reads."""
    plt = _plt()
    cats = [(f"{family.candidate}\ncounted", counts.candidate_identity, COL["kept"]),
            (f"{family.candidate}\ndiscarded", counts.discarded_identity, COL["discarded"])]
    for k, g in enumerate(family.controls):
        cats.append((f"{g}\ncontrol", counts.identity.get(g, []), CONTROL_COLS[k % len(CONTROL_COLS)]))
    fig, ax = plt.subplots(figsize=(1.3 * len(cats) + 2, 3.8))
    rng = random.Random(2)
    data = []
    for i, (name, vals, col) in enumerate(cats):
        xs = [i + rng.uniform(-0.18, 0.18) for _ in vals]
        ax.scatter(xs, vals, s=10, color=col, alpha=0.6, linewidths=0)
        if vals:
            med = sorted(vals)[len(vals) // 2]
            ax.plot([i - 0.28, i + 0.28], [med, med], color="#222222", lw=1.4)
        ax.text(i, 101.5, f"n={len(vals)}", ha="center", fontsize=7, color="#444444")
        for v in vals:
            data.append({"group": name.replace("\n", " "), "identity": v})
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([c[0] for c in cats], fontsize=8)
    ax.set_xlim(-0.6, len(cats) - 0.4)
    ax.set_ylim(min([v for _, vals, _ in cats for v in vals] + [60]) - 5, 104)
    ax.set_ylabel("Identity to the query (%)")
    ax.set_title("Are candidate reads as similar to the queries as control reads?" + (f"  {label}" if label else ""))
    fig.text(0.01, -0.06, "Lower than the controls: a divergent candidate. Near 100% with few reads: "
             "possible contamination. Very low and discarded: fragments of paralogs.", fontsize=7, color="#555555")
    return _write_tsv(_save(fig, out), data, ["group", "identity"])
