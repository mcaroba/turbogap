#!/usr/bin/env python3
"""Put each keyword family's handlers in the same order on both branches.

The families already hold the same keywords; they were emitted in whatever order
the original chain happened to use, which differs between the trees, so `diff`
pairs unrelated handlers and reports divergence that is not there --
read_options_vdw read as 100 diverged lines while holding the same keywords.

Sorting alphabetically within each family is safe because, after the duplicates
below are dropped, every keyword appears in exactly one handler: nothing depends
on which branch of the chain is tested first.

Duplicate handlers are removed first. Each duplicate's copies were confirmed
character-identical, so the later one was unreachable and dropping it cannot
change behaviour.

Edits are collected against the original text and applied last-to-first, so no
offset is ever used after the string it indexes has been rewritten. Doing it the
other way round silently leaves some families unsorted.
"""
import re

PATH = 'src/read_files.f90'
src = open(PATH, encoding='utf-8').read()

SUB = re.compile(r'^( *)subroutine (read_options_\w+)\(', re.M)
HEAD = re.compile(r"^( *)(if|else if) \(keyword == ['\"]([a-z0-9_]+)['\"]")

edits = []
moved = dropped_total = 0

for m in SUB.finditer(src):
    name = m.group(2)
    body_start = src.index('\n', m.end()) + 1
    body_end = src.index(f'end subroutine {name}', body_start)
    lines = src[body_start:body_end].splitlines(keepends=True)

    first = IND = None
    for i, l in enumerate(lines):
        h = HEAD.match(l)
        if h:
            first, IND = i, h.group(1)
            break
    if first is None:
        continue
    close = next(i for i in range(first, len(lines))
                 if lines[i].rstrip('\n') == IND + 'else')

    handlers, cur = [], None
    for i in range(first, close):
        h = HEAD.match(lines[i])
        if h and h.group(1) == IND:
            if cur:
                handlers.append((cur[0], i, cur[1]))
            cur = (i, h.group(3))
    if cur:
        handlers.append((cur[0], close, cur[1]))

    seen, kept, dropped = set(), [], 0
    for s, e, kw in handlers:
        if kw in seen:
            dropped += 1
            continue
        seen.add(kw)
        kept.append((s, e, kw))
    dropped_total += dropped

    before = [k for _, _, k in kept]
    kept.sort(key=lambda h: h[2])
    if [k for _, _, k in kept] != before:
        moved += 1

    rebuilt = []
    for n, (s, e, _) in enumerate(kept):
        chunk = list(lines[s:e])
        h = chunk[0]
        if n == 0:
            h = re.sub(r'^( *)else if ', r'\1if ', h)
        elif not h.lstrip().startswith('else if'):
            h = re.sub(r'^( *)if ', r'\1else if ', h, count=1)
        chunk[0] = h
        rebuilt.extend(chunk)

    new_body = ''.join(lines[:first]) + ''.join(rebuilt) + ''.join(lines[close:])
    edits.append((body_start, body_end, new_body))

for a, b, text in reversed(edits):
    src = src[:a] + text + src[b:]

open(PATH, 'w', encoding='utf-8').write(src)
print(f'{len(edits)} families processed, {moved} reordered, '
      f'{dropped_total} duplicate handlers dropped')
