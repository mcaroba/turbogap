#!/usr/bin/env python3
"""Classify a driver block's identifiers into arguments / locals / free.

An identifier that is declared in the driver and used both inside and outside
the block must be an argument. Declared in the driver and used only inside ->
becomes a local of the extracted procedure. Not declared in the driver at all
-> comes from a USEd module and is free.
"""
import re
import sys

path, lo, hi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
src = open(path, encoding='utf-8', errors='replace').read().splitlines()

KW = set('''if then else elseif end do while call subroutine function return program module use
implicit none integer real double precision complex logical character dimension allocatable
pointer target intent in out inout parameter save data external interface contains type class
select case default exit cycle go to continue stop write read print open close allocate
deallocate nullify associated present size shape lbound ubound reshape trim adjustl adjustr len
index scan verify max min abs sqrt exp log sin cos tan sum product matmul dot_product transpose
all any count maxval minval maxloc minloc mod modulo nint floor ceiling int dble cmplx sign huge
tiny epsilon true false and or not eq ne lt le gt ge eqv neqv only public private optional result
recursive elemental pure where forall associate block procedure allocated null merge pack unpack
spread transfer achar iachar char ichar repeat kind c_ptr c_size_t c_int c_double c_char c_bool
c_loc c_null_ptr c_f_pointer dp sp'''.split())

ID = re.compile(r'\b([a-zA-Z_][a-zA-Z0-9_]*)\b')
DECL = re.compile(r'^\s*(integer|real|double\s+precision|complex|logical|character|type\s*\(|class\s*\(|procedure)', re.I)

# --- find the driver's declaration section: from `implicit none` to the first
#     executable statement is unreliable, so instead treat any line that looks
#     like a declaration (incl. its continuations) as declaring its names.
declared = set()
i = 0
while i < len(src):
    line = src[i].split('!', 1)[0]
    if DECL.match(line):
        chunk = line
        while chunk.rstrip().endswith('&') and i + 1 < len(src):
            i += 1
            chunk += src[i].split('!', 1)[0]
        rhs = chunk.split('::', 1)[1] if '::' in chunk else ''
        # strip initialisers and array specs
        rhs = re.sub(r'=[^,]*', '', rhs)
        rhs = re.sub(r'\([^)]*\)', '', rhs)
        for n in ID.finditer(rhs):
            declared.add(n.group(1).lower())
    i += 1


def ids(seq, skip_decl):
    """Count identifier uses. When skip_decl, skip declaration statements --
    including their continuation lines, which do not themselves start with a
    type keyword and would otherwise read as real uses."""
    out = {}
    in_decl_cont = False
    for raw in seq:
        line = raw.split('!', 1)[0]
        if line.lstrip().startswith('#'):
            continue
        is_decl_start = bool(DECL.match(line))
        if skip_decl and (is_decl_start or in_decl_cont):
            in_decl_cont = line.rstrip().endswith('&')
            continue
        in_decl_cont = False
        for m in ID.finditer(line):
            n = m.group(1).lower()
            if n in KW or n.isdigit():
                continue
            out[n] = out.get(n, 0) + 1
    return out


inb = ids(src[lo - 1:hi], True)
out = ids(src[:lo - 1] + src[hi:], True)

args, locals_, free = [], [], []
for n in sorted(inb):
    if n not in declared:
        free.append(n)
    elif out.get(n, 0) > 0:
        args.append(n)
    else:
        locals_.append(n)

print(f'block {lo}-{hi} ({hi-lo+1} lines)')
print(f'  ARGUMENTS (declared in driver, used outside too) : {len(args)}')
print(f'  LOCALS    (declared in driver, only used here)   : {len(locals_)}')
print(f'  FREE      (from USEd modules)                    : {len(free)}')
print('\n--- ARGUMENTS ---')
print('  ' + ', '.join(args))
print('\n--- LOCALS (move into the procedure) ---')
print('  ' + ', '.join(locals_))
print('\n--- FREE (sample) ---')
print('  ' + ', '.join(free[:60]))
