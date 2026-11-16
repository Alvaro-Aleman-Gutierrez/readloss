#!/usr/bin/env python3
"""Asignacion de lecturas (DIAMOND blastx vs familia DCX/DCLK) y cociente observado/esperado de DCX.
Regla fijada con controles (2026-11-07): alineamiento >= 30 aa; una lectura cuenta como DCX solo si su mejor hit es DCX,
supera al mejor hit no-DCX en >= 5 bits y tiene identidad >= 75 %. Controles monocopia: dclk1a + dclk2a
(dclk3 excluido: DIAMOND le asigna exceso de lecturas). Esperado = (lecturas control / aa control) x 357 aa.
uso: parse_dmnd.py <tsv> <etiqueta>"""
import sys, collections
from scipy.stats import poisson
DCX = {"DCX","dcx","LOC103178947","LOC113084161","LOC103035707","LOC102692282","LOC102366980"}
L = {"DCX": 357, "dclk1a": 770, "dclk2a": 826}
def gen(s):
    g = s.split("_")[1]
    if g in DCX: return "DCX"
    if s.startswith(("Danrer_dclk1a","Astmex_dclk1a","Lepocu_dclk1a","Calmil_dclk1a","Latcha_DCLK1A")): return "dclk1a"
    if "dclk2a" in s or "DCLK2A" in s: return "dclk2a"
    return "other"
f, tag = sys.argv[1], sys.argv[2]
hits = collections.defaultdict(list)
for l in open(f):
    q, s, pid, ln, bs = l.split()
    if int(ln) >= 30: hits[q].append((float(bs), gen(s), float(pid)))
cnt = collections.Counter(); pidd = []
for q, h in hits.items():
    h.sort(reverse=True); b = h[0]
    if b[1] == "DCX":
        otro = max((x[0] for x in h if x[1] != "DCX"), default=0)
        if b[0] - otro < 5 or b[2] < 75: cnt["DCX_descartada"] += 1; continue
        pidd.append(b[2])
    cnt[b[1]] += 1
esp = (cnt["dclk1a"] + cnt["dclk2a"]) / (L["dclk1a"] + L["dclk2a"]) * L["DCX"]; o = cnt["DCX"]
r = o / esp if esp else float("nan")
p = poisson.cdf(o, esp) if esp else float("nan")
lo = poisson.ppf(0.025, 1) if False else None
print(f"{tag}\tDCX={o}\tdclk1a={cnt['dclk1a']}\tdclk2a={cnt['dclk2a']}\tDCX_descartadas={cnt['DCX_descartada']}"
      f"\tesperado={esp:.1f}\tobs/esp={r:.3f}\tP(<=obs|presente)={p:.2g}\tpid_DCX={sum(pidd)/len(pidd) if pidd else 0:.1f}")
