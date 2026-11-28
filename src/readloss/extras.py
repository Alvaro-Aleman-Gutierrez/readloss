"""Panel of single-copy genes, GC bias and sex inference (version 0.2.0).

Panel. The expected count of the family test rests on two or three controls of the same
family. A panel of unrelated single-copy genes (for example a subset of BUSCO genes of
the clade) gives a second, independent expectation. Its median read rate per aa is
robust to a few genes that are duplicated, lost or divergent in the lineage tested.
Comparing the control rate with the panel rate shows controls that are not single copy
(about twice the panel after a genome duplication) or that collect foreign reads.

GC. Illumina libraries made with PCR under-cover GC-rich DNA, the main reason why genes
are missing from short-read assemblies of birds. The same bias lowers the read count of
a GC-rich candidate. With a panel, the relation between read GC and coverage is
estimated in each run from the panel genes, and the coverage expected at the GC of the
candidate is compared with that at the GC of the controls.

Sex. Markers are genes of known location on the X or Z chromosome. Their reads relative
to the autosomal controls are about 1 in the homogametic sex and 0.5 in the heterogametic
sex, which tells the sex of a specimen whose sex was not reported.
"""
import statistics

PANEL_MIN_LEN = 100
GC_MIN_READS = 20
GC_BIN = 0.03
GC_MIN_GENES_BIN = 3


def _mean(v):
    return sum(v) / len(v) if v else float("nan")


def panel_rate(counts, family):
    """Median reads per aa across panel genes (genes without reads count as zero)."""
    rates = [counts.by_gene.get(g, 0) / L for g, L in family.panel.items() if L >= PANEL_MIN_LEN]
    return statistics.median(rates) if rates else float("nan"), len(rates)


def gc_curve(counts, family):
    """Binned median coverage (reads per aa, relative to the overall median) against read GC."""
    pts = []
    for g, L in family.panel.items():
        gc = counts.gc.get(g, [])
        if len(gc) >= GC_MIN_READS and L >= PANEL_MIN_LEN:
            pts.append((_mean(gc), counts.by_gene.get(g, 0) / L))
    if len(pts) < 2 * GC_MIN_GENES_BIN:
        return []
    med = statistics.median(r for _, r in pts)
    if med <= 0:
        return []
    bins = {}
    for gc, r in pts:
        bins.setdefault(int(gc / GC_BIN), []).append(r / med)
    return sorted(((k + 0.5) * GC_BIN, statistics.median(v)) for k, v in bins.items()
                  if len(v) >= GC_MIN_GENES_BIN)


def curve_at(curve, gc):
    """Linear interpolation of the curve; nan outside its range."""
    if not curve or gc != gc:
        return float("nan")
    if gc < curve[0][0] - GC_BIN / 2 or gc > curve[-1][0] + GC_BIN / 2:
        return float("nan")                     # no extrapolation beyond half a bin
    if gc <= curve[0][0]:
        return curve[0][1]
    if gc >= curve[-1][0]:
        return curve[-1][1]
    for (x0, y0), (x1, y1) in zip(curve, curve[1:]):
        if x0 <= gc <= x1:
            return y0 + (y1 - y0) * (gc - x0) / (x1 - x0) if x1 > x0 else y0
    return curve[-1][1] if gc == curve[-1][0] else float("nan")


def infer_sex(counts, family, ctrl_rate):
    """Sex from the X or Z markers: (sex, marker ratio) or ('unknown', nan)."""
    if not family.markers or not ctrl_rate or ctrl_rate != ctrl_rate or ctrl_rate <= 0:
        return "unknown", float("nan")
    obs = sum(counts.by_gene.get(m, 0) for m in family.markers)
    exp = ctrl_rate * sum(family.lengths[m] * family.assignable.get(m, 1.0) * family.copies.get(m, 1.0)
                          for m in family.markers)
    if exp < 10:
        return "unknown", obs / exp if exp else float("nan")
    r = obs / exp
    loc = family.location[family.markers[0]]
    homo, hetero = ("female", "male") if loc == "X" else ("male", "female")
    if r >= 0.75:
        return homo, r
    if r <= 0.65:
        return hetero, r
    return "unknown", r
