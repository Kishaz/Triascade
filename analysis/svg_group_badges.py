#!/usr/bin/env python3
"""Group each numbered badge (its circle(s) + number text) into a <g> so a badge
can be selected and moved as ONE object in Inkscape/Illustrator.

svglite emits each badge as 1-2 <circle> lines immediately followed by the number
<text> line (text-anchor middle, ~6px, content = a 1-2 digit number). This wraps
those lines in <g class="badge" id="badge-N-k">. Rewrites in place; idempotent.

Usage: python3 svg_group_badges.py figure.svg
"""
import re, sys

TEXT_RE = re.compile(r"<text\b[^>]*>(.*?)</text>")
CIRCLE_RE = re.compile(r"^\s*<circle\b")


def badge_number(line):
    if "text-anchor='middle'" not in line and 'text-anchor="middle"' not in line:
        return None
    if "font-size: 6." not in line and "font-size:6." not in line:
        return None
    m = TEXT_RE.search(line)
    if not m:
        return None
    inner = re.sub(r"<[^>]+>", "", m.group(1)).strip()
    return inner if re.fullmatch(r"\d{1,2}", inner) else None


def main(path):
    txt = open(path, encoding="utf-8").read()
    if "class='badge'" in txt or 'class="badge"' in txt:
        print("badges already grouped; skipping")
        return
    out, grouped, k = [], 0, 0
    for line in txt.split("\n"):
        num = badge_number(line)
        if num is not None:
            circ = []
            while out and CIRCLE_RE.search(out[-1]):
                circ.insert(0, out.pop())
            k += 1
            out.append(f"<g class='badge' id='badge-{num}-{k}'>")
            out.extend(circ)
            out.append(line)
            out.append("</g>")
            grouped += 1
        else:
            out.append(line)
    open(path, "w", encoding="utf-8").write("\n".join(out))
    print(f"grouped {grouped} badges")


if __name__ == "__main__":
    main(sys.argv[1])
