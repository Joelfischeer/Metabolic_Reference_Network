"""
PNG versions of the workflow diagram and the methods table (same content and
design as workflow_diagram.tex / methods_table.tex, rendered with matplotlib
because no LaTeX compiler is needed for this).

Icons come from Font Awesome (fontawesomefree package), so run with:
    uv run --with fontawesomefree python methods_figures/make_pngs.py
"""
from pathlib import Path

import fontawesomefree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
FA = next((Path(fontawesomefree.__file__).parent / "static" / "fontawesomefree" / "webfonts").glob("fa-solid-900.ttf"))
FA_FP = FontProperties(fname=str(FA))
DPI = 300
plt.rcParams["font.family"] = "DejaVu Sans"

# ------------------------------------------------------------------ numbers ---
N = dict(organs=14, pairs=91, all="22,451", healthy="1,585", obese="3,594", both=200,
         organ_papers=700, years_pairs=20, years_organs=5)

COL = dict(db="#475569", filt="#2563A8", cond="#0F766E", key="#B45309", typ="#7C3AED",
           desc="#BE185D", organ="#15803D", viz="#1E293B", line="#64748B")
ICON = dict(database="", filter="", sliders="", lungs="", key="",
            tags="", comment="", robot="", diagram="")


def mix(hex_color, amount):
    """Blend a colour with white (amount = share of the colour)."""
    h = hex_color.lstrip("#")
    rgb = [int(h[i:i + 2], 16) for i in (0, 2, 4)]
    return tuple((amount * c + (1 - amount) * 255) / 255 for c in rgb)


def wrap(ax, text, fontsize, width_cm, weight="normal"):
    """Greedy word wrap using the real rendered text width (data units = cm)."""
    r = ax.figure.canvas.get_renderer()
    out = []
    for para in text.split("\n"):
        line = ""
        for w in para.split(" "):
            cand = (line + " " + w).strip()
            t = ax.text(0, 0, cand, fontsize=fontsize, fontweight=weight)
            wpx = t.get_window_extent(r).width
            t.remove()
            if wpx / ax.figure.dpi * 2.54 <= width_cm or not line:
                line = cand
            else:
                out.append(line)
                line = w
        out.append(line)
    return "\n".join(out)


# ============================================================ flow diagram ====

def workflow():
    fig = plt.figure(figsize=(28.7 / 2.54, 17.6 / 2.54), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-3.8, 24.9)
    ax.set_ylim(-14.4, 3.2)
    ax.axis("off")
    ax.set_aspect("equal")
    fig.canvas.draw()

    BOX_W, BOX_H = 6.1, 3.7

    def box(cx, cy, color, icon, title, body, w=BOX_W, h=BOX_H):
        ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                    boxstyle="round,pad=0,rounding_size=0.25",
                                    fc=mix(color, 0.07), ec=color, lw=1.6, zorder=3))
        ax.text(cx, cy + h / 2 - 0.85, wrap(ax, title, 10.5, w - 0.6, "bold"), ha="center", va="center",
                fontsize=10.5, fontweight="bold", color="#0f172a", zorder=5)
        ax.text(cx, cy + h / 2 - 1.5, wrap(ax, body, 8.3, w - 0.7), ha="center", va="top",
                fontsize=8.3, color="#1e293b", zorder=5, linespacing=1.35)
        x0, y0 = cx - w / 2 + 0.75, cy + h / 2
        ax.add_patch(Circle((x0, y0), 0.52, fc=color, ec="white", lw=2, zorder=6))
        ax.text(x0, y0, icon, fontproperties=FA_FP, fontsize=14, color="white",
                ha="center", va="center", zorder=7)

    def arrow(p0, p1):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=16, lw=2.0,
                                     color=COL["line"], shrinkA=0, shrinkB=0, zorder=2))

    def line(p0, p1):
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=COL["line"], lw=2.0, zorder=2, solid_capstyle="butt")

    X = [0, 6.9, 13.8, 20.7]
    Y = [0, -5.9, -11.9]

    # row 1
    box(X[0], Y[0], COL["db"], ICON["database"], "PubMed search",
        f"One query per organ pair, last {N['years_pairs']} years: organ MeSH terms and aliases AND "
        "22 metabolic AND 21 crosstalk terms (title/abstract)")
    box(X[1], Y[0], COL["filt"], ICON["filter"], "Papers per organ pair",
        "Local same-sentence co-location and multi-organ filter\n"
        f"n = {N['all']} unique papers\n({N['pairs']} pairs, {N['organs']} organs)")
    box(X[2], Y[0], COL["cond"], ICON["sliders"], "Condition subsets",
        f"A) all papers: {N['all']}\nB) healthy keywords: {N['healthy']}\nC) obese keywords: {N['obese']}\n"
        f"({N['both']} papers are in both B and C)")
    arrow((X[0] + BOX_W / 2, Y[0]), (X[1] - BOX_W / 2, Y[0]))
    arrow((X[1] + BOX_W / 2, Y[0]), (X[2] - BOX_W / 2, Y[0]))

    # row 2
    box(X[0], Y[1], COL["organ"], ICON["lungs"], "Organ literature",
        f"Per organ: MeSH/aliases AND metabolic terms, last {N['years_organs']} years\n"
        f"n = {N['organ_papers']} papers\n({N['organs']} organs × 50)")
    box(X[1], Y[1], COL["key"], ICON["key"], "Key player extraction",
        "Hormones, metabolites and proteins matched in title, abstract and MeSH terms; synonyms merged")
    box(X[2], Y[1], COL["typ"], ICON["tags"], "Connection types (LLM)",
        "llama3.2, 5 independent runs, majority vote: up to 3 of 11 connection types per pair")
    box(X[3], Y[1], COL["desc"], ICON["comment"], "Edge descriptions (LLM)",
        "25 random papers → 5 selected → cited summary of 3–4 sentences")

    # dashed container around the three per-pair analyses
    pad = 0.55
    left, right = X[1] - BOX_W / 2 - pad, X[3] + BOX_W / 2 + pad
    top, bottom = Y[1] + BOX_H / 2 + pad, Y[1] - BOX_H / 2 - pad
    ax.add_patch(FancyBboxPatch((left, bottom), right - left, top - bottom,
                                boxstyle="round,pad=0,rounding_size=0.3", fc="#00000008",
                                ec=COL["line"], lw=1.4, ls=(0, (5, 4)), zorder=1))
    ax.text(left + 0.15, bottom - 0.2, "For each organ–organ pair", ha="left", va="top",
            fontsize=9.5, style="italic", color=COL["line"])

    arrow((X[0], Y[0] - BOX_H / 2), (X[0], Y[1] + BOX_H / 2))
    bus_y = Y[0] - BOX_H / 2 - 1.0
    line((X[2], Y[0] - BOX_H / 2), (X[2], bus_y))
    line((X[1], bus_y), (X[3], bus_y))
    for x in (X[1], X[2], X[3]):
        arrow((x, bus_y), (x, Y[1] + BOX_H / 2))

    # row 3
    box(X[0], Y[2], COL["organ"], ICON["robot"], "Organ descriptions (LLM)",
        "6 papers → description of five sentences: substrates, energy pathways and hormonal regulation")
    box(X[2], Y[2], COL["viz"], ICON["diagram"], "Visualization",
        f"Interactive HTML dashboard (Cytoscape.js):\n{N['organs']} organs and {N['pairs']} connections, "
        "views for all papers, healthy and obese", w=9.4)
    arrow((X[0], Y[1] - BOX_H / 2), (X[0], Y[2] + BOX_H / 2))
    arrow((X[0] + BOX_W / 2, Y[2]), (X[2] - 9.4 / 2, Y[2]))
    arrow((X[2], bottom), (X[2], Y[2] + BOX_H / 2))

    out = HERE / "workflow_diagram.png"
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    print("[ok]", out.name)


# ================================================================== table =====

HEAD = ["Step", "Input", "Processing", "Use in the HTML dashboard"]
RATIO = [1.25, 1.7, 2.85, 1.85]
ROWS = [
    (
        '1  Literature search and filtering',
        '14 organs, 91 organ pairs; organ MeSH terms and aliases; 22 metabolic and 21 crosstalk keywords',
        'One PubMed query per pair (title/abstract), last 20 years, up to 20,000 hits, 0.4 s between requests; PMID, title, abstract, year and MeSH descriptors fetched via EFetch. Local rule-based filters keep a paper only if both organs and a crosstalk keyword occur in the same sentence (whole-word, negation-aware) and no third organ indicates a different relationship (multi-organ disambiguation). 187,757 hits retrieved (summed over pairs); 27,416 pair–paper assignments kept, i.e. 22,451 unique papers (87 of 91 pairs with at least one paper).',
        'Existence of a connection and its paper count; overlap-coefficient (OO) based edge weight; paper list with abstracts per connection; Literature Statistics tab (Sankey and bar plots).',
    ),
    (
        '2  Condition subselection',
        'Papers of step 1; 7 healthy and 11 obese keywords',
        'A paper is kept for a condition if its title or abstract contains at least one keyword of that condition (whole-word, negation-aware, anywhere in the text). The two sets overlap (200 papers). Healthy: 1,585 papers in 70 pairs; obese: 3,594 papers in 72 pairs. Computed on the fly while the HTML is built.',
        'Healthy/Obese switch: a connection is shown if it has at least one matching paper; paper counts, key players, OO weights and statistics are recomputed from the filtered papers.',
    ),
    (
        '3  Key player extraction',
        'Title, abstract and MeSH descriptors of the papers of each pair and condition',
        'Multi-pattern (Aho–Corasick) whole-word matching against curated vocabularies of hormones (171), metabolites (254) and proteins/transporters (45,138). Synonyms and paralog families are merged to one display name at match time; each term counts once per paper (count = number of papers).',
        'Ranked lists of hormones, metabolites and proteins with paper counts in the information panel of each connection; recalculated for the Healthy and Obese views.',
    ),
    (
        '4  Connection-type classification (LLM)',
        'First 100 papers of a pair (abstracts truncated to 400 characters); definitions of 11 connection types',
        'Local LLM (llama3.2 via Ollama, temperature 0.5) assigns 1–3 types per run; 5 independent runs; types named in the majority of runs (at least 3 of 5) are kept; if none reaches a majority the single most frequent type is kept.',
        'Connection-type label and colour of each edge and the type legend; Healthy/Obese views use the condition-specific classification where it exists and the all-papers classification otherwise.',
    ),
    (
        '5  Connection descriptions (LLM)',
        'Papers of a pair (condition-filtered for the condition files)',
        '25 papers drawn at random (all if fewer); the LLM selects the 5 most relevant; a second prompt writes a 3–4 sentence summary of the metabolic and hormonal basis of the axis that cites every selected paper.',
        'Cited summary with the five source papers in the information panel of each connection.',
    ),
    (
        '6  Organ descriptions (LLM)',
        '14 organs; organ MeSH term/aliases and the 22 metabolic keywords',
        'Separate PubMed search per organ (title/abstract, last 5 years, up to 50 papers; 700 in total); the first six papers (abstracts truncated to 400 characters) are given to the LLM, which writes exactly five sentences on consumed and produced substrates, energy pathways and hormonal regulation. Refusals are detected and re-prompted with a simplified prompt.',
        'Description and source papers in the information panel of each organ (node).',
    ),
    (
        '7  Assembly of the dashboard',
        'Caches of steps 1, 4, 5 and 6 (plus the condition-specific files); PET partial-correlation networks of the healthy and obese cohorts',
        'Undirected NetworkX graph (14 nodes, 91 edges) with all edge and node attributes, including the precomputed condition variants; export to a single HTML file rendered with Cytoscape.js (version 3.29.2, loaded from a CDN). The cohort partial-correlation networks (16 healthy and 9 obese edges) are added as comparison networks.',
        'Interactive network with All/Healthy/Obese views, comparison with the PET-based networks, and the Literature and Edge Statistics tabs.',
    ),
]
CAPTION = ("Pipeline steps and their contribution to the interactive HTML dashboard. Numbers refer to the "
           "current run (14 organs, 91 organ pairs); papers are unique papers unless stated otherwise.")


def table():
    W = 7.0                                     # inches: A4 portrait text width at 100%
    scale = W / sum(RATIO)
    colw = [r * scale for r in RATIO]
    fs, fs_head = 6.5, 7.2
    pad_x, pad_y, line_h = 0.05, 0.05, 0.113    # inches

    fig = plt.figure(figsize=(W, 12), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(-12, 0)
    ax.axis("off")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()

    def wrap_in(text, width_in, size, weight="normal", family="DejaVu Sans"):
        out = []
        for para in text.split("\n"):
            line = ""
            for w in para.split(" "):
                cand = (line + " " + w).strip()
                t = ax.text(0, 0, cand, fontsize=size, fontweight=weight, family=family)
                wpx = t.get_window_extent(r).width
                t.remove()
                if wpx / DPI <= width_in or not line:
                    line = cand
                else:
                    out.append(line)
                    line = w
            out.append(line)
        return out

    cap_lines = wrap_in(CAPTION, W - 0.3, 7.6, "bold")
    y = -0.2
    cap_h = len(cap_lines) * 0.14
    heights, wrapped = [], []
    for row in ROWS:
        cells = [wrap_in(c, colw[i] - 2 * pad_x, fs, "bold" if i == 0 else "normal")
                 for i, c in enumerate(row)]
        wrapped.append(cells)
        heights.append(max(len(c) for c in cells) * line_h + 2 * pad_y + 0.04)
    head_h = 0.26
    total = 0.2 + cap_h + 0.15 + head_h + sum(heights) + 0.2

    # caption
    ax.text(0.2, y, "\n".join(cap_lines), ha="left", va="top", fontsize=7.6, fontweight="bold", color="#0f172a",
            linespacing=1.3)
    y -= cap_h + 0.15
    ax.plot([0.1, W - 0.1], [y, y], color="#0f172a", lw=1.8)                    # top rule
    xs = [0.1]
    for w_ in colw:
        xs.append(xs[-1] + w_)
    xs = [x * (W - 0.2) / xs[-1] + 0.1 - 0.1 * (W - 0.2) / xs[-1] for x in xs]
    xs = [0.1 + (W - 0.2) * sum(RATIO[:i]) / sum(RATIO) for i in range(len(RATIO) + 1)]
    for i, h in enumerate(HEAD):
        ax.text(xs[i] + pad_x, y - head_h / 2, h, ha="left", va="center", fontsize=fs_head, fontweight="bold",
                color="#0f172a")
    y -= head_h
    ax.plot([0.1, W - 0.1], [y, y], color="#0f172a", lw=0.9)                    # mid rule
    for ri, (cells, h) in enumerate(zip(wrapped, heights)):
        if ri % 2 == 0:
            ax.add_patch(plt.Rectangle((0.1, y - h), W - 0.2, h, fc="#f1f5f9", ec="none", zorder=0))
        for ci, lines in enumerate(cells):
            ax.text(xs[ci] + pad_x, y - pad_y, "\n".join(lines), ha="left", va="top", fontsize=fs,
                    fontweight="bold" if ci == 0 else "normal", family="DejaVu Sans",
                    color="#0f172a" if ci == 0 else "#1e293b", linespacing=1.28)
        y -= h
    ax.plot([0.1, W - 0.1], [y, y], color="#0f172a", lw=1.8)                    # bottom rule

    used = -y + 0.15
    fig.set_size_inches(W, used)
    ax.set_ylim(-used, 0)
    out = HERE / "methods_table.png"
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    print("[ok]", out.name)


if __name__ == "__main__":
    workflow()
    table()
