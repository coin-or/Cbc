#!/usr/bin/env python3
"""Generate an HTML parameter reference page from 'mipster -dumpParameters' JSON.

Usage:
    python3 doc/generate_html_reference.py [path/to/mipster] > docs/parameters.html
"""

import subprocess, sys, os, json, re
from html import escape
from datetime import datetime


def load_params(cbc):
    out = subprocess.run([cbc, "-dumpParameters"],
                         capture_output=True, text=True, timeout=10).stdout
    start = out.index("["); end = out.rindex("]") + 1
    return json.loads(out[start:end])


def get_version(cbc):
    out = subprocess.run([cbc, "--version"],
                         capture_output=True, text=True, timeout=10).stdout
    m = re.search(r"(?:MIPster|CBC)\s+(\S+)", out, re.IGNORECASE)
    return m.group(1) if m else "devel"


def group_by_topic(params):
    topics = {}
    for p in params:
        t = p.get("topic", "")
        if not t or p.get("displayPriority", 0) == 0:
            continue
        topics.setdefault(t, []).append(p)
    return topics


ORDER = [
    "Stopping", "Cuts", "Heuristics",
    "MIP Preprocessing", "MIP Preprocessing \u2014 Bound Propagation", "LP Presolve",
    "Branching", "Tolerances", "Conflict Graph", "Strategy", "Solving",
    "Simplex", "Barrier", "Scaling", "Output", "I/O", "Parallelism",
]

CONSTRUCTIVE_HEURISTICS = {
    "feasibilityPump", "DivingCoefficient", "DivingFractional",
    "DivingGuided", "DivingLineSearch", "DivingPseudocost",
    "DivingVectorLength", "DivingSome", "roundingHeuristic",
    "greedyHeuristic", "naiveHeuristics", "pivotAndFix",
    "randomizedRounding", "Rens",
}
IMPROVEMENT_HEURISTICS = {
    "Rins", "VndVariableNeighborhoodSearch", "Dins",
    "proximitySearch", "dwHeuristic", "localTreeSearch",
}
IMPROVEMENT_2_HEURISTICS = {"combineSolutions", "combine2Solutions"}


def topic_anchor(topic):
    return re.sub(r"[^a-z0-9]+", "-", topic.lower()).strip("-")


def humanize_bound(v):
    s = str(v)
    if s == "2147483647":  return "INT_MAX"
    if s == "-2147483647": return "\u2212INT_MAX"
    if s in ("1e+60", "1e+20"):   return "\u221e"
    if s in ("-1e+60", "-1e+20"): return "\u2212\u221e"
    return s


BADGE_CLASS = {
    "integer": "badge-int", "double": "badge-dbl", "keyword": "badge-kw",
    "string": "badge-str", "action": "badge-act",
    "file": "badge-file", "directory": "badge-file",
}
TYPE_LABEL = {
    "integer": "Integer", "double": "Double", "keyword": "Keyword",
    "string": "String", "action": "Action", "file": "File", "directory": "Directory",
}


def param_card(p):
    name   = p["name"]
    ptype  = p.get("type", "")
    short  = p.get("shortHelp", "")
    long_h = p.get("longHelp", "")
    kws    = p.get("keywords", [])
    dv     = p.get("defaultValue", "")

    badge = ""
    if ptype in BADGE_CLASS:
        badge = (f'<span class="badge {BADGE_CLASS[ptype]}">'
                 f'{TYPE_LABEL.get(ptype, ptype)}</span>')

    meta = ""
    if dv not in ("", None):
        meta = f'<span class="pmeta">default:&nbsp;<code>{escape(str(dv))}</code></span>'

    body_parts = []

    if long_h and long_h != short:
        text = escape(long_h).replace("\\n", "<br>")
        body_parts.append(f'<p class="plong">{text}</p>')

    if kws:
        kw_html = " ".join(f'<code class="kw">{escape(k)}</code>' for k in kws)
        dv_note = (f'&ensp;<span class="pmeta">(default: <code>{escape(str(dv))}</code>)</span>'
                   if dv not in ("", None) else "")
        body_parts.append(f'<div class="pextra"><strong>Values:</strong> {kw_html}{dv_note}</div>')

    if ptype == "integer":
        lo = humanize_bound(p.get("lowerInt", ""))
        hi = humanize_bound(p.get("upperInt", ""))
        dv_note = (f'&ensp;<span class="pmeta">(default:&nbsp;{escape(str(dv))})</span>'
                   if dv not in ("", None) else "")
        body_parts.append(
            f'<div class="pextra"><strong>Range:</strong> '
            f'{escape(lo)}&thinsp;&ndash;&thinsp;{escape(hi)}{dv_note}</div>')
    elif ptype == "double":
        lo = humanize_bound(p.get("lowerDbl", ""))
        hi = humanize_bound(p.get("upperDbl", ""))
        dv_note = (f'&ensp;<span class="pmeta">(default:&nbsp;{escape(str(dv))})</span>'
                   if dv not in ("", None) else "")
        body_parts.append(
            f'<div class="pextra"><strong>Range:</strong> '
            f'{escape(lo)}&thinsp;&ndash;&thinsp;{escape(hi)}{dv_note}</div>')

    # data-search is used by JS: name + short help, lowercased for fast matching
    search_text = escape((name + " " + short).lower())
    data = f'data-search="{search_text}"'

    summary = (f'<span class="pname"><code>-{escape(name)}</code></span>'
               f'{badge}'
               f'<span class="pshort">{escape(short)}</span>'
               f'{meta}')

    if body_parts:
        body = "\n".join(body_parts)
        return (f'<details class="param" {data}>\n'
                f'  <summary>{summary}</summary>\n'
                f'  <div class="pbody">{body}</div>\n'
                f'</details>')
    else:
        return (f'<div class="param param-flat" {data}>\n'
                f'  {summary}\n'
                f'</div>')


def section_html(topic, params):
    anchor = topic_anchor(topic)
    cards = "\n".join(param_card(p) for p in params)
    return (f'<section id="{anchor}" class="topic-section">\n'
            f'  <h2>{escape(topic)}</h2>\n'
            f'{cards}\n'
            f'</section>\n')


def heuristics_section(params):
    constructive, improvement, improvement2, general = [], [], [], []
    for p in params:
        n = p["name"]
        if n in CONSTRUCTIVE_HEURISTICS:    constructive.append(p)
        elif n in IMPROVEMENT_HEURISTICS:   improvement.append(p)
        elif n in IMPROVEMENT_2_HEURISTICS: improvement2.append(p)
        else:                               general.append(p)

    anchor = topic_anchor("Heuristics")
    parts = [f'<section id="{anchor}" class="topic-section">\n  <h2>Heuristics</h2>\n']

    groups = [
        ("Constructive Heuristics", constructive,
         "Do <strong>not</strong> require an existing feasible solution — "
         "they construct one from scratch."),
        ("Improvement Heuristics", improvement,
         "Require <strong>at least one</strong> feasible solution; "
         "they improve upon the incumbent."),
        ("Improvement Heuristics (2+ solutions)", improvement2,
         "Require <strong>at least two</strong> feasible solutions; "
         "they combine or crossover multiple solutions."),
        ("General Heuristic Settings", general, ""),
    ]
    for title, group, intro in groups:
        if not group:
            continue
        parts.append(f'<div class="subsection">')
        parts.append(f'<h3>{title}</h3>')
        if intro:
            parts.append(f'<p class="sub-intro">{intro}</p>')
        for p in group:
            parts.append(param_card(p))
        parts.append('</div>')

    parts.append('</section>')
    return "\n".join(parts)


def ordered_topics(topics):
    seen = set()
    result = []
    for t in ORDER:
        if t in topics:
            result.append(t)
            seen.add(t)
    for t in topics:
        if t not in seen:
            result.append(t)
    return result


CSS = """\
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

    :root {
      --accent: #2563eb;
      --bg: #0f172a;
      --surface: #1e293b;
      --surface2: #253047;
      --text: #e2e8f0;
      --muted: #94a3b8;
      --border: #334155;
      --radius: 10px;
      --nav-h: 56px;
      --sb-w: 230px;
    }

    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      background: var(--bg);
      color: var(--text);
      line-height: 1.6;
    }
    a { color: var(--accent); text-decoration: none; }
    a:hover { text-decoration: underline; }
    code { font-family: "SF Mono", "Cascadia Code", Consolas, monospace; font-size: 0.875em; }

    /* ── Nav ── */
    nav.topnav {
      position: sticky; top: 0; z-index: 100;
      background: var(--surface);
      border-bottom: 1px solid var(--border);
      padding: 0 2rem;
      display: flex; align-items: center; justify-content: space-between;
      height: var(--nav-h);
    }
    .nav-brand { font-weight: 700; font-size: 1.1rem; color: var(--text); text-decoration: none; }
    .nav-links { display: flex; gap: 1.5rem; font-size: 0.9rem; }
    .nav-links a { color: var(--muted); }
    .nav-links a:hover { color: var(--text); text-decoration: none; }

    /* ── Two-column layout ── */
    .page-layout {
      display: flex;
      align-items: flex-start;
      max-width: 1260px;
      margin: 0 auto;
      padding: 0 1.5rem;
    }

    /* ── Sidebar ── */
    .sidebar {
      width: var(--sb-w);
      flex-shrink: 0;
      position: sticky;
      top: var(--nav-h);
      max-height: calc(100vh - var(--nav-h));
      overflow-y: auto;
      padding: 1.5rem 1rem 2rem 0;
      border-right: 1px solid var(--border);
    }
    .sidebar::-webkit-scrollbar { width: 4px; }
    .sidebar::-webkit-scrollbar-thumb { background: var(--border); border-radius: 2px; }

    /* search box + result count */
    .search-wrap {
      position: relative;
      margin-bottom: 0.35rem;
    }
    .search-wrap input {
      display: block; width: 100%;
      background: var(--surface2);
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 0.45rem 2rem 0.45rem 0.75rem;
      color: var(--text);
      font-size: 0.85rem;
      outline: none;
    }
    .search-wrap input:focus { border-color: var(--accent); }
    .search-clear {
      position: absolute; right: 0.5rem; top: 50%; transform: translateY(-50%);
      background: none; border: none; color: var(--muted);
      cursor: pointer; font-size: 1rem; line-height: 1; padding: 0;
      display: none;
    }
    .search-clear.visible { display: block; }
    .search-hint {
      font-size: 0.72rem; color: var(--muted);
      margin-bottom: 0.75rem; padding-left: 0.1rem;
    }
    .search-hint kbd {
      background: var(--surface2); border: 1px solid var(--border);
      border-radius: 3px; padding: 0 0.3rem; font-size: 0.7rem;
    }
    #match-count {
      font-size: 0.72rem; color: var(--muted);
      margin-bottom: 0.5rem; padding-left: 0.1rem; min-height: 1em;
    }

    /* sidebar nav links */
    .sb-link {
      display: flex; justify-content: space-between; align-items: center;
      padding: 0.28rem 0.5rem;
      border-radius: 6px;
      font-size: 0.82rem;
      color: var(--muted);
      transition: background .1s, color .1s, opacity .1s;
    }
    .sb-link:hover { background: var(--surface2); color: var(--text); text-decoration: none; }
    .sb-link.active { color: var(--text); background: var(--surface2); font-weight: 600; }
    .sb-link.dimmed { opacity: 0.3; pointer-events: none; }
    .sb-count {
      font-size: 0.7rem;
      background: var(--surface2);
      border-radius: 999px;
      padding: 0 0.45rem;
      color: var(--muted);
      min-width: 1.6rem;
      text-align: center;
      flex-shrink: 0;
    }

    /* ── Main content ── */
    .main-content {
      flex: 1;
      min-width: 0;
      padding: 2rem 0 4rem 2.5rem;
    }
    .page-header { margin-bottom: 2.5rem; }
    .page-header h1 { font-size: 2rem; font-weight: 800; margin-bottom: 0.4rem; }
    .page-header p { color: var(--muted); font-size: 0.9rem; }

    /* ── Topic sections ── */
    .topic-section {
      margin-bottom: 3rem;
      /* scroll-margin offsets the sticky nav when jumping to an anchor */
      scroll-margin-top: calc(var(--nav-h) + 0.75rem);
    }
    .topic-section h2 {
      font-size: 1.35rem; font-weight: 700;
      padding-bottom: 0.5rem;
      border-bottom: 1px solid var(--border);
      margin-bottom: 1rem;
    }
    .subsection { margin: 1.25rem 0; }
    .subsection h3 {
      font-size: 0.95rem; font-weight: 600; color: var(--muted);
      margin-bottom: 0.35rem; text-transform: uppercase; letter-spacing: 0.05em;
    }
    .sub-intro { font-size: 0.85rem; color: var(--muted); margin-bottom: 0.65rem; }

    /* ── Parameter cards ── */
    .param {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 8px;
      margin-bottom: 0.4rem;
      overflow: hidden;
      transition: border-color .12s;
    }
    .param:hover { border-color: #475569; }

    /* expandable params */
    details.param > summary {
      display: flex; align-items: center; flex-wrap: wrap; gap: 0.5rem;
      padding: 0.6rem 1rem;
      cursor: pointer;
      list-style: none;
      user-select: none;
    }
    details.param > summary::-webkit-details-marker { display: none; }
    details.param > summary::before {
      content: "\\25B6";
      font-size: 0.55rem;
      color: var(--muted);
      margin-right: 0.3rem;
      transition: transform .15s;
      flex-shrink: 0;
    }
    details.param[open] > summary::before { transform: rotate(90deg); }
    details.param[open] > summary { border-bottom: 1px solid var(--border); background: var(--surface2); }

    /* flat (no-expand) params */
    .param-flat {
      display: flex; align-items: center; flex-wrap: wrap; gap: 0.5rem;
      padding: 0.6rem 1rem 0.6rem 2rem;
    }

    /* param body (expanded content) */
    .pbody {
      padding: 0.8rem 1rem 0.9rem 2.35rem;
      font-size: 0.875rem;
    }
    .pbody .plong {
      color: var(--muted); margin-bottom: 0.5rem; line-height: 1.6;
    }
    .pbody .pextra { margin-top: 0.4rem; font-size: 0.83rem; color: var(--muted); }
    .pbody .pextra strong { color: var(--text); }
    .pbody .kw {
      background: var(--surface2); border-radius: 4px;
      padding: 0.1rem 0.35rem;
    }

    /* summary typography */
    .pname code { font-size: 0.9rem; font-weight: 600; color: var(--text); }
    .pshort { font-size: 0.84rem; color: var(--muted); flex: 1; min-width: 0; }
    .pmeta  { font-size: 0.76rem; color: var(--muted); margin-left: auto; white-space: nowrap; }
    .pmeta code { color: #93c5fd; }

    /* ── Type badges ── */
    .badge {
      display: inline-block; border-radius: 4px;
      padding: 0.08rem 0.4rem; font-size: 0.68rem;
      font-weight: 700; letter-spacing: 0.04em; flex-shrink: 0;
    }
    .badge-int  { background: #1e3a5f; color: #93c5fd; }
    .badge-dbl  { background: #2e1a5e; color: #c4b5fd; }
    .badge-kw   { background: #064e3b; color: #6ee7b7; }
    .badge-str  { background: #451a03; color: #fcd34d; }
    .badge-act  { background: #450a0a; color: #fca5a5; }
    .badge-file { background: #1e293b; color: #94a3b8; }

    /* ── No-results message ── */
    #no-results {
      display: none; text-align: center;
      color: var(--muted); padding: 4rem 2rem; font-size: 1rem;
    }

    /* ── Footer ── */
    footer {
      text-align: center;
      border-top: 1px solid var(--border);
      padding: 2rem; font-size: 0.85rem; color: var(--muted);
    }

    @media (max-width: 768px) {
      .sidebar { display: none; }
      .main-content { padding-left: 0; }
    }
"""

JS = """\
  const searchInput = document.getElementById('param-search');
  const clearBtn    = document.getElementById('search-clear');
  const matchCount  = document.getElementById('match-count');
  const noResults   = document.getElementById('no-results');
  const params      = Array.from(document.querySelectorAll('.param'));
  const sections    = Array.from(document.querySelectorAll('.topic-section'));
  const sbLinks     = Array.from(document.querySelectorAll('.sb-link'));

  function applyFilter(q) {
    let visible = 0;
    params.forEach(el => {
      const hit = !q || (el.dataset.search || '').includes(q);
      el.style.display = hit ? '' : 'none';
      if (hit) visible++;
    });

    // hide subsections that have no visible params
    document.querySelectorAll('.subsection').forEach(sub => {
      const anyVisible = sub.querySelectorAll('.param').some(
        el => el.style.display !== 'none'
      );
      sub.style.display = anyVisible ? '' : 'none';
    });

    // hide/show topic sections and dim sidebar links
    sections.forEach(sec => {
      const anyVisible = sec.querySelectorAll('.param').some(
        el => el.style.display !== 'none'
      );
      sec.style.display = anyVisible ? '' : 'none';
      const link = document.querySelector(`.sb-link[href="#${sec.id}"]`);
      if (link) link.classList.toggle('dimmed', !anyVisible);
    });

    // update counters
    if (q) {
      matchCount.textContent = `${visible} of ${params.length} parameters`;
    } else {
      matchCount.textContent = '';
    }
    noResults.style.display = (q && visible === 0) ? 'block' : 'none';
    clearBtn.classList.toggle('visible', q.length > 0);
  }

  searchInput.addEventListener('input', () => applyFilter(searchInput.value.trim().toLowerCase()));

  clearBtn.addEventListener('click', () => {
    searchInput.value = '';
    applyFilter('');
    searchInput.focus();
  });

  // Press "/" to focus the search box (unless already focused in an input)
  document.addEventListener('keydown', e => {
    if (e.key === '/' && document.activeElement !== searchInput &&
        document.activeElement.tagName !== 'INPUT') {
      e.preventDefault();
      searchInput.focus();
      searchInput.select();
    }
    if (e.key === 'Escape' && document.activeElement === searchInput) {
      searchInput.value = '';
      applyFilter('');
      searchInput.blur();
    }
  });

  // Highlight active sidebar link based on scroll position
  const observer = new IntersectionObserver(entries => {
    entries.forEach(e => {
      if (e.isIntersecting) {
        sbLinks.forEach(a => a.classList.remove('active'));
        const link = document.querySelector(`.sb-link[href="#${e.target.id}"]`);
        if (link) {
          link.classList.add('active');
          // scroll the link into view inside the sidebar if needed
          link.scrollIntoView({ block: 'nearest' });
        }
      }
    });
  }, { rootMargin: '-8% 0px -82% 0px' });
  sections.forEach(s => observer.observe(s));
"""


def emit(version, topics, out):
    date  = datetime.now().strftime("%B %Y")
    total = sum(len(v) for v in topics.values())
    order = ordered_topics(topics)

    sidebar_links = "\n".join(
        f'<a href="#{topic_anchor(t)}" class="sb-link">'
        f'{escape(t)}<span class="sb-count">{len(topics[t])}</span></a>'
        for t in order
    )

    sections = []
    for t in order:
        if t == "Heuristics":
            sections.append(heuristics_section(topics[t]))
        else:
            sections.append(section_html(t, topics[t]))
    content = "\n".join(sections)

    out.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>MIPster \u2014 Parameter Reference</title>
  <style>
{CSS}  </style>
</head>
<body>

<nav class="topnav">
  <a class="nav-brand" href="index.html">\u26a1 MIPster</a>
  <div class="nav-links">
    <a href="index.html#download">Download</a>
    <a href="index.html#features">Features</a>
    <a href="index.html#quickstart">Quick Start</a>
    <a href="parameters.html" style="color:var(--text)">Parameters</a>
    <a href="https://github.com/h-g-s/mipster">GitHub</a>
  </div>
</nav>

<div class="page-layout">

  <aside class="sidebar">
    <div class="search-wrap">
      <input type="search" id="param-search" autocomplete="off" spellcheck="false"
             placeholder="Search parameters\u2026" aria-label="Search parameters">
      <button id="search-clear" class="search-clear" aria-label="Clear search">&times;</button>
    </div>
    <div class="search-hint">Press <kbd>/</kbd> to focus &nbsp;&middot;&nbsp; <kbd>Esc</kbd> to clear</div>
    <div id="match-count"></div>
    <nav id="sidebar-nav">
{sidebar_links}
    </nav>
  </aside>

  <div class="main-content">
    <header class="page-header">
      <h1>Parameter Reference</h1>
      <p>MIPster {escape(version)} &mdash; {escape(date)} &mdash; {total} parameters across {len(order)} sections</p>
    </header>

{content}

    <div id="no-results">
      <p>No parameters match your search.<br>
      <small>Try a shorter term, or press <strong>Esc</strong> to clear.</small></p>
    </div>
  </div>

</div>

<footer>
  <p>MIPster is a fork of <a href="https://github.com/coin-or/Cbc">COIN-OR CBC</a>
  &middot; Licensed under <a href="https://github.com/h-g-s/mipster/blob/master/LICENSE">EPL-2.0</a>
  &middot; <a href="https://github.com/h-g-s/mipster">GitHub</a></p>
</footer>

<script>
{JS}</script>
</body>
</html>
""")


def main():
    cbc = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/prog/cbc/bin/mipster")
    params  = load_params(cbc)
    version = get_version(cbc)
    topics  = group_by_topic(params)
    total   = sum(len(v) for v in topics.values())
    print(f"Parsed {len(topics)} topics, {total} parameters", file=sys.stderr)
    emit(version, topics, sys.stdout)


if __name__ == "__main__":
    main()
