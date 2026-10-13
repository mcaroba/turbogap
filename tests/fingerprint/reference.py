"""The fingerprint written out with Python's unbounded integers, from the
REF/HASH lines fingerprintverify prints. Exit 1 on any mismatch."""
import sys

MASK = 0xFFFFFFFF
M1 = 1540483477
M2 = 668265261


def lane_step(h, w, m):
    k = w & MASK
    k = (k * m) & MASK
    k ^= k >> 24
    k = (k * m) & MASK
    return ((h * m) & MASK) ^ k


def lane_final(h, m):
    h ^= h >> 13
    h = (h * m) & MASK
    return h ^ (h >> 15)


def fingerprint(dims, words):
    h1 = (2166136261 ^ dims[0]) & MASK
    h2 = (3266489917 ^ dims[1]) & MASK
    h1 = lane_step(h1, dims[1], M1)
    h2 = lane_step(h2, dims[0], M2)
    for w in words:
        h1 = lane_step(h1, w, M1)
        h2 = lane_step(h2, w, M2)
    h = (lane_final(h1, M1) << 32) | lane_final(h2, M2)
    return h - (1 << 64) if h >= 1 << 63 else h


lines = sys.stdin.read().splitlines()
bad = 0
for ref, got in zip(lines[0::2], lines[1::2]):
    f = ref.split()
    dims = [int(f[2]), int(f[3])]
    words = [int(x) for x in f[4:]]
    want = fingerprint(dims, words)
    have = int(got.split()[1])
    print(f"{f[1]}: fortran {have} python {want} {'ok' if have == want else 'MISMATCH'}")
    bad += have != want
sys.exit(1 if bad else 0)
