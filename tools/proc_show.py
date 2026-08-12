#!/usr/bin/env python3
"""Show the whitespace/comment-insensitive diff of one procedure between files."""
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


def norm(lines):
    out = []
    for l in lines:
        c = l.split('!', 1)[0]
        c = re.sub(r'\s+', '', c).lower()
        if c:
            out.append(c)
    return out


A = procs(sys.argv[1])
B = procs(sys.argv[2])
p = sys.argv[3].lower()
limit = int(sys.argv[4]) if len(sys.argv) > 4 else 120
d = list(difflib.unified_diff(norm(A[p]), norm(B[p]), 'master', 'gpu', n=1, lineterm=''))
print('\n'.join(d[:limit]))
print(f'\n[{len(d)} diff lines total]')
