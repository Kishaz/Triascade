#!/usr/bin/env python3
"""Make the text of a figure SVG editable: one text object per label, not one per word.

Two things make svglite output awkward to edit by hand in Inkscape or Illustrator:
  1. ggtext / gridtext draw rich-text labels word by word, so every word is its own <text>;
  2. every <text> carries a fixed width (textLength), so an edited string is squeezed or
     stretched to the old width.
This rewrites the SVG in place. Words are joined with ordinary spaces into one <text> (a nested
<tspan> only where the style changes: bold, italic, colour, size), and all fixed widths are
removed so edited text reflows. Positions come from the rendered words, so nothing moves.

Usage:  python3 svg_join_text.py figure.svg [labels.json]

Without labels.json (any figure): consecutive words on one baseline that follow each other at
  word spacing are joined into one line of text. Plain labels are left as they are, minus the
  fixed width. Rotated labels are handled.
With labels.json (figures drawn with dlab(), which records its labels in drawing order as
  [{"txt": "<markup>", "hj": 0 | 0.5 | 1}, ...]): each label becomes ONE text object, its lines
  as <tspan sodipodi:role="line"> (Inkscape treats these as the lines of one text), anchored
  start / middle / end as it was drawn, with the drawn line height stated. The SVG must then
  contain no text other than those labels.
A file already processed is left alone (it carries a marker comment).
"""
import sys, re, json, html, math

MARK = "<!-- text joined by svg_join_text.py -->"
svg_path = sys.argv[1]; labels_path = sys.argv[2] if len(sys.argv) > 2 else None
s = open(svg_path, encoding="utf8").read()
if MARK in s:
    print(f"{svg_path}: already joined, left as it is"); sys.exit(0)

TEXT = re.compile(r"<text((?:\s+[\w:-]+='[^']*')*)\s*>(.*?)</text>", re.S)
ATTR = re.compile(r"([\w:-]+)='([^']*)'")
class W:                                             # one <text> as svglite wrote it
    def __init__(self, m):
        self.m = m; self.a = dict(ATTR.findall(m.group(1))); self.txt = m.group(2)
        self.len = float(self.a.get("textLength", "0px").replace("px", ""))
        self.style = style_of(self.a.get("style", "")); self.size = float(self.style.get("font-size", "10px").replace("px", ""))
        tr = self.a.get("transform")
        if tr:
            t = re.fullmatch(r"translate\(([-\d.]+),([-\d.]+)\) rotate\(([-\d.]+)\)", tr)
            assert t, f"unexpected transform: {tr}"
            self.x, self.y, self.ang = float(t.group(1)), float(t.group(2)), float(t.group(3))
        else:
            self.x, self.y, self.ang = float(self.a["x"]), float(self.a["y"]), 0.0
def style_of(st):
    d = {}
    for part in st.split(";"):
        if ":" in part:
            k, v = part.split(":", 1); d[k.strip()] = v.strip()
    return d
def css(d): return " ".join(f"{k}: {v};" for k, v in d.items())
def place(w):                                        # the position attributes of a text starting where w starts
    return (f"transform='translate({w.x:.2f},{w.y:.2f}) rotate({w.ang:g})'" if w.ang else f"x='{w.x:.2f}' y='{w.y:.2f}'")
def gap(a, b):                                       # distance from the end of a to the start of b, along and across the baseline
    t = math.radians(a.ang); dx, dy = math.cos(t), math.sin(t)
    ex, ey = a.x + dx * a.len, a.y + dy * a.len
    return (b.x - ex) * dx + (b.y - ey) * dy, abs(-(b.x - ex) * dy + (b.y - ey) * dx)
def runs_body(ws, spaces, base_css, base):           # words -> text with a tspan where the style differs from the base
    runs = []
    for w, sp in zip(ws, spaces):
        st = css(w.style)
        if runs and runs[-1][0] == st: runs[-1][1] += (" " if sp else "") + w.txt
        else: runs.append([st, w.txt, sp])
    body = ""
    for k, (st, txt, sp) in enumerate(runs):
        # the space between two runs is set in the main style where one of them has it
        if k > 0 and sp and st == base_css: txt = " " + txt
        if k < len(runs) - 1 and runs[k + 1][2] and runs[k + 1][0] != base_css: txt += " "
        if st == base_css: body += txt
        else:
            d = style_of(st); diff = {a: b for a, b in d.items() if base.get(a) != b}
            for a in base:                           # a property the base has and the run lacks: reset it
                if a not in d: diff[a] = {"font-weight": "normal", "font-style": "normal", "fill": "#000000"}.get(a, "inherit")
            body += f"<tspan style='{css(diff)}'>{txt}</tspan>"
    return body
def base_style(ws):                                  # the style most of the characters are set in
    weight = {}
    for w in ws: weight[css(w.style)] = weight.get(css(w.style), 0) + len(w.txt)
    b = max(weight, key=weight.get); return b, style_of(b)

words = [W(m) for m in TEXT.finditer(s)]
assert len(words) == s.count("<text"), "a text element is not in the form svglite writes"
new = {}; drop = set(); n_obj = 0                    # replacement for the first word of a group; the other words go

if labels_path:                                      # ---- one object per recorded label
    labels = json.load(open(labels_path, encoding="utf8"))
    def lines_of(txt):
        t = re.sub(r"<span[^>]*>|</span>", "", txt)
        return [html.unescape(ln.replace("**", "").replace("*", "")) for ln in re.split(r"<br\s*/?>", t)]
    it = iter(words)
    for lab in labels:
        lines = []
        for ln in lines_of(lab["txt"]):
            target = re.sub(r"\s+", "", ln); acc = ""; ws = []
            while acc != target:
                w = next(it); acc += html.unescape(w.txt); ws.append(w)
                assert target.startswith(acc), f"label text does not match the SVG: {ln!r} vs {acc!r}"
                assert not w.ang and "text-anchor" not in w.a
            lines.append(ws)
        allw = [w for ws in lines for w in ws]
        base_css, base = base_style(allw)
        anchor = {0: "start", 0.5: "middle", 1: "end"}[lab["hj"]]
        parts = []
        for ws in lines:
            left = ws[0].x; right = ws[-1].x + ws[-1].len
            x = {0: left, 0.5: (left + right) / 2, 1: right}[lab["hj"]]
            assert all(abs(w.y - ws[0].y) < 0.01 for w in ws), "words of one line are on different baselines"
            sp = [False] + [gap(a, b)[0] > 0.12 * max(a.size, b.size) for a, b in zip(ws, ws[1:])]
            parts.append((x, ws[0].y, runs_body(ws, sp, base_css, base)))
        st_text = f"{base_css} text-anchor: {anchor};"
        if len(parts) == 1:
            x, y, body = parts[0]
            el = f"<text xml:space='preserve' x='{x:.2f}' y='{y:.2f}' style='{st_text}'>{body}</text>"
        else:
            # Inkscape re-spaces the lines of a text object from its line-height, so state the one drawn
            pitch = (parts[-1][1] - parts[0][1]) / (len(parts) - 1)
            st_text += f" line-height: {pitch / float(base['font-size'].replace('px', '')):.4f};"
            inner = "".join(f"<tspan sodipodi:role='line' x='{x:.2f}' y='{y:.2f}'>{body}</tspan>" for x, y, body in parts)
            el = f"<text xml:space='preserve' x='{parts[0][0]:.2f}' y='{parts[0][1]:.2f}' style='{st_text}'>{inner}</text>"
        new[allw[0].m.start()] = el; n_obj += 1
        drop.update(w.m.start() for w in allw[1:])
    assert next(it, None) is None, "the SVG has text that is not in the label list"
else:                                                # ---- join the words of each line; keep plain labels
    groups = []
    for w in words:
        g = groups[-1] if groups else None
        ok = False
        if g and "text-anchor" not in w.a and "text-anchor" not in g[-1].a and w.ang == g[-1].ang and g[-1].len > 0:
            along, across = gap(g[-1], w); big = max(w.size, g[-1].size)
            ok = across < 0.03 * big and -0.06 * big <= along <= 0.5 * big
        if ok: g.append(w)
        else: groups.append([w])
    for g in groups:
        if len(g) == 1:
            w = g[0]; a = {k: v for k, v in w.a.items() if k not in ("textLength", "lengthAdjust")}
            el = "<text " + " ".join(f"{k}='{v}'" for k, v in a.items()) + f">{w.txt}</text>"
        else:
            base_css, base = base_style(g)
            sp = [False] + [gap(a, b)[0] > 0.12 * max(a.size, b.size) for a, b in zip(g, g[1:])]
            el = f"<text xml:space='preserve' {place(g[0])} style='{base_css}'>{runs_body(g, sp, base_css, base)}</text>"
            drop.update(w.m.start() for w in g[1:])
        new[g[0].m.start()] = el; n_obj += 1

out = []; pos = 0
for w in words:
    out.append(s[pos:w.m.start()]); pos = w.m.end()
    if w.m.start() in new: out.append(new[w.m.start()])
    elif s[pos:pos + 1] == "\n": pos += 1            # a dropped word leaves no blank line behind
out.append(s[pos:])
s2 = "".join(out)
if "sodipodi:role" in s2:
    s2 = s2.replace("<svg xmlns='http://www.w3.org/2000/svg'", "<svg xmlns='http://www.w3.org/2000/svg' xmlns:sodipodi='http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd'", 1)
i = s2.rindex("</svg>"); s2 = s2[:i] + MARK + "\n" + s2[i:]
assert "textLength" not in s2
open(svg_path, "w", encoding="utf8").write(s2)
print(f"{svg_path}: {len(words)} text elements -> {n_obj} editable text objects ({'by label' if labels_path else 'by line'})")
