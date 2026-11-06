#!/usr/bin/env python3
"""Compose the editable-text Fig 3 / Fig 4 SVGs from their panel SVGs.

The R composite scripts (make_fig3_composite.R / make_fig4_composite.R) only emit
raster PNG previews. This reassembles the same layouts from the panels' vector SVGs
so Fig3.svg / Fig4.svg keep editable text (panel letters, axis labels, numbers).
Only the atlas's anatomy icons remain raster (they are raster in the source panel).
Run after regenerating any panel:  python3 results/figures/final/compose_svgs.py
Requires: svgutils (pip install svgutils).
"""
import svgutils.transform as st
from lxml import etree as le
import xml.etree.ElementTree as ET
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.normpath(os.path.join(HERE, "..", "composite")) + "/"
OUT = HERE + "/"
MM2PT = 72 / 25.4

def intrinsic(stem):
    r = ET.parse(C + stem + ".svg").getroot()
    return float(r.get("width").replace("pt", "")), float(r.get("height").replace("pt", ""))

def asp(stem):
    w, h = intrinsic(stem)
    return w / h

def panel(stem, x_mm, y_mm, w_mm):
    iw, ih = intrinsic(stem)
    root = st.fromfile(C + stem + ".svg").getroot()
    root.moveto(x_mm * MM2PT, y_mm * MM2PT, scale_x=(w_mm * MM2PT) / iw)
    return root

def save(elems, W_mm, H_mm, out):
    fig = st.SVGFigure()
    fig.set_size((f"{W_mm*MM2PT:.2f}pt", f"{H_mm*MM2PT:.2f}pt"))
    fig.append(st.FigureElement(le.Element("rect", {
        "x": "0", "y": "0", "width": f"{W_mm*MM2PT}", "height": f"{H_mm*MM2PT}", "fill": "white"})))
    for e in elems:
        fig.append(e)
    fig.save(out)
    print(f"wrote {os.path.relpath(out)}  ({W_mm:.0f} x {H_mm:.0f} mm)")

# ---- Figure 3: a flow banner + 2x2 (b atlas | d dose / c violin | e alluvial) ----
W, gap = 234, 4
bh = W / asp("ICI_Fig3_flow")
wa, wc, wb, wd = 128, W - 4 - 128, 146, W - 4 - 146
hA, hC = wa / asp("ICI_Fig3_atlas"), wc / asp("ICI_Fig3c_doseresponse")
hB, hD = wb / asp("ICI_Fig3b_separation"), wd / asp("ICI_Fig_archetypes")
row1, row2 = max(hA, hC), max(hB, hD)
H = bh + gap + row1 + gap + row2
y1, y2 = bh + gap, bh + gap + row1 + gap
save([panel("ICI_Fig3_flow", 0, 0, W),
      panel("ICI_Fig3_atlas", 0, y1, wa),
      panel("ICI_Fig3c_doseresponse", W - wc, y1, wc),
      panel("ICI_Fig3b_separation", 0, y2, wb),
      panel("ICI_Fig_archetypes", W - wd, y2, wd)], W, H, OUT + "Fig3.svg")

# ---- Figure 4: stack a ranking / b heatmap / c grids, full width ----
W4, g4 = 174, 2
hA4, hB4, hC4 = (W4 / asp("ICI_Fig4_a_ranking"), W4 / asp("ICI_Fig4_b_heatmap"),
                W4 / asp("ICI_Fig4_c_grids"))
H4 = hA4 + hB4 + hC4 + 2 * g4
save([panel("ICI_Fig4_a_ranking", 0, 0, W4),
      panel("ICI_Fig4_b_heatmap", 0, hA4 + g4, W4),
      panel("ICI_Fig4_c_grids", 0, hA4 + g4 + hB4 + g4, W4)], W4, H4, OUT + "Fig4.svg")
