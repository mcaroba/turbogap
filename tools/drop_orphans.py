#!/usr/bin/env python3
"""Remove named entities from declaration statements in turbogap.f90.

Declarations here are shared: one statement declares many names across several
continuation lines, so a whole-line delete would take live names with it. This
rewrites each affected statement with the dropped entities removed, and deletes
the statement outright only when nothing is left.
"""
import re
import sys

PATH = 'src/turbogap.f90'
DROP = set(n.strip().lower() for n in open(sys.argv[1]).read().split(',') if n.strip())

src = open(PATH, encoding='utf-8', errors='replace').read().splitlines()
DECL = re.compile(r'^\s*(integer|real|double\s+precision|complex|logical|character|type\s*\(|class\s*\()', re.I)


def entities(rhs):
    out, depth, cur = [], 0, ''
    for ch in rhs:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(cur)
            cur = ''
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return [e.strip() for e in out if e.strip()]


def wrap(prefix, ents, indent):
    """Re-emit `prefix :: a, b, c` wrapped at a sensible width."""
    line = f'{indent}{prefix} :: '
    out, cur = [], line
    for n, e in enumerate(ents):
        piece = e + (', ' if n < len(ents) - 1 else '')
        if len(cur) + len(piece) > 96 and cur.strip() != (prefix + ' ::'):
            out.append(cur.rstrip() + ' &')
            cur = indent + '     ' + piece
        else:
            cur += piece
    out.append(cur.rstrip())
    return out


result, i, changed, removed_total, deleted = [], 0, 0, 0, 0
while i < len(src):
    raw = src[i]
    line = raw.split('!', 1)[0]
    if not DECL.match(line):
        result.append(raw)
        i += 1
        continue

    start = i
    chunk = line.rstrip()
    while chunk.endswith('&') and i + 1 < len(src):
        i += 1
        chunk = chunk[:-1] + src[i].split('!', 1)[0].strip().lstrip('&').rstrip()
    i += 1

    if '::' not in chunk:
        result.extend(src[start:i])
        continue

    prefix, rhs = chunk.split('::', 1)
    ents = entities(rhs)
    keep = []
    dropped_here = 0
    for e in ents:
        m = re.match(r'([a-zA-Z_][a-zA-Z0-9_]*)', e)
        if m and m.group(1).lower() in DROP:
            dropped_here += 1
        else:
            keep.append(e)

    if dropped_here == 0:
        result.extend(src[start:i])
        continue

    changed += 1
    removed_total += dropped_here
    indent = re.match(r'\s*', src[start]).group(0)
    if keep:
        result.extend(wrap(prefix.strip(), keep, indent))
    else:
        deleted += 1
    # keep any trailing comment from the first line, if it carried one
    if '!' in src[start] and keep:
        pass

open(PATH, 'w', encoding='utf-8').write('\n'.join(result) + '\n')
print(f'statements rewritten: {changed} (of which deleted entirely: {deleted})')
print(f'entities removed    : {removed_total} (requested {len(DROP)})')
assert removed_total == len(DROP), 'removed count does not match the request'
