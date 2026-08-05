#!/usr/bin/env python3
"""Find declarations in turbogap.f90 that nothing references any more.

A name counts as referenced if it appears anywhere outside its own declaration
statement -- including inside OTHER declaration statements, where it may be an
array bound or a kind. That is the N_CONTRIB trap: a parameter used only as a
dimension looks unused to a naive scan.
"""
import re
import sys

PATH = sys.argv[1] if len(sys.argv) > 1 else 'src/turbogap.f90'
WANT = set(n.strip().lower() for n in sys.argv[2].split(',')) if len(sys.argv) > 2 else None

src = open(PATH, encoding='utf-8', errors='replace').read().splitlines()
DECL = re.compile(r'^\s*(integer|real|double\s+precision|complex|logical|character|type\s*\(|class\s*\()', re.I)
ID = re.compile(r'\b([a-zA-Z_][a-zA-Z0-9_]*)\b')

# Join declaration statements, remembering their line span.
stmts, other = [], []
i = 0
while i < len(src):
    line = src[i].split('!', 1)[0]
    if DECL.match(line):
        start = i
        chunk = line.rstrip()
        while chunk.endswith('&') and i + 1 < len(src):
            i += 1
            chunk = chunk[:-1] + src[i].split('!', 1)[0].strip().lstrip('&').rstrip()
        stmts.append((start, i, chunk))
    else:
        other.append(line)
    i += 1


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


declared = {}   # name -> index into stmts
for idx, (a, b, chunk) in enumerate(stmts):
    if '::' not in chunk:
        continue
    for ent in entities(chunk.split('::', 1)[1]):
        m = re.match(r'([a-zA-Z_][a-zA-Z0-9_]*)', ent)
        if m:
            declared.setdefault(m.group(1).lower(), idx)

exec_text = '\n'.join(other).lower()
exec_names = set(x.lower() for x in ID.findall(exec_text))

orphans = []
for name, own in sorted(declared.items()):
    if WANT and name not in WANT:
        continue
    if name in exec_names:
        continue
    # referenced from a different declaration statement?
    used_in_other_decl = False
    for idx, (a, b, chunk) in enumerate(stmts):
        if idx == own:
            continue
        if re.search(r'\b' + re.escape(name) + r'\b', chunk, re.I):
            used_in_other_decl = True
            break
    if used_in_other_decl:
        continue
    orphans.append(name)

print(f'{len(orphans)} orphaned declaration names')
for n in orphans:
    a, b, chunk = stmts[declared[n]]
    print(f'  {n:<40} L{a+1}-{b+1}')
