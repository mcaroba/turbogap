#!/usr/bin/env python3
"""Compare two Fortran files procedure by procedure.

Reports, for each procedure present in both, its length on each side and the
real (whitespace-insensitive) diff between them -- so the merge surface can be
attributed to individual procedures rather than to the file as a whole.
"""
import re
import sys
import difflib

RE_START = re.compile(
    r'^\s*(?:pure\s+|elemental\s+|recursive\s+)*'
    r'(?:subroutine|(?:(?:real|integer|logical|complex|double\s+precision|character|type\s*\([^)]*\))\s*(?:\([^)]*\)\s*)?)?function)\s+'
    r'([a-zA-Z_]\w*)', re.I)
RE_END = re.compile(r'^\s*end\s+(subroutine|function)\b', re.I)


def procs(path):
    out, cur, buf = {}, None, []
    for raw in open(path, encoding='utf-8', errors='replace'):
        line = raw.rstrip('\n')
        code = line.split('!', 1)[0]
        if cur is None:
            m = RE_START.match(code)
            if m and not code.lstrip().lower().startswith('end'):
                cur, buf = m.group(1).lower(), [line]
        else:
            buf.append(line)
            if RE_END.match(code):
                out[cur] = buf
                cur, buf = None, []
    return out


a, b = sys.argv[1], sys.argv[2]
A, B = procs(a), procs(b)
shared = sorted(set(A) & set(B))
onlyA = sorted(set(A) - set(B))
onlyB = sorted(set(B) - set(A))


def norm(lines):
    out = []
    for l in lines:
        c = l.split('!', 1)[0]
        c = re.sub(r'\s+', '', c).lower()   # strip ALL whitespace, as diff -w does
        if c:
            out.append(c)
    return out


print(f'{"procedure":<44} {"A":>6} {"B":>6} {"realdiff":>9}')
print('-' * 70)
rows = []
for p in shared:
    d = list(difflib.unified_diff(norm(A[p]), norm(B[p]), n=0))
    changed = sum(1 for l in d
                  if (l.startswith('+') or l.startswith('-'))
                  and not l.startswith('+++') and not l.startswith('---'))
    rows.append((changed, p, len(A[p]), len(B[p])))
for changed, p, la, lb in sorted(rows, key=lambda r: -r[0]):
    flag = '  <-- identical' if changed == 0 else ''
    print(f'{p:<44} {la:>6} {lb:>6} {changed:>9}{flag}')

print(f'\nshared procedures: {len(shared)}   '
      f'identical ignoring whitespace/comments: {sum(1 for r in rows if r[0]==0)}')
print(f'total real diff across shared procedures: {sum(r[0] for r in rows)}')
print(f'\nA-only ({len(onlyA)}): {", ".join(onlyA)}')
print(f'\nB-only ({len(onlyB)}): {", ".join(onlyB)}')
print(f'B-only line count: {sum(len(B[p]) for p in onlyB)}')
