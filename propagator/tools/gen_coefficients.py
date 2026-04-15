#!/usr/bin/env python3
"""
Download GRGM1200A gravity model truncated to 100×100 and emit a compact
binary blob for inclusion into the WASM binary via Rust's include_bytes!.

Output binary layout (little-endian):
  bytes 0-3:   N_MAX as uint32_le  (100)
  bytes 4+:    for n in 0..=N_MAX, for m in 0..=n:
                   C_nm  as float64_le   (8 bytes)
                   S_nm  as float64_le   (8 bytes)

Total for N=100: 4 + 5151*16 = 82420 bytes ≈ 80 KB.

GRGM1200A SHA tab line format:
   n  m  C_nm  S_nm  [sigma_C  sigma_S]

Reference values (Zuber et al. 2013, Lemoine et al. 2014) used as fallback:
  GM      = 4902.800066 km³/s²  (compiled into Rust, not stored in binary)
  R_ref   = 1738.0 km            (compiled into Rust, not stored in binary)
"""

import struct, sys, os, urllib.request

N_MAX = 100

# Primary source — GRAIL mission data, PDS Geosciences Node (2024 URL)
URLS = [
    # Correct PDS URL (grail-l not gra-l)
    "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200a_sha.tab",
    "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200b_sha.tab",
]

CACHE = "/tmp/grgm1200a.tab"

# Hardcoded published GRGM1200A coefficients (fully normalised) for fallback.
# Sources: Lemoine et al. 2014 (DOI 10.1002/2014JE004619); file gggrx_1200a_sha.tab.
# NOTE: values are 10× smaller than the GRAIL_1200C set — correct canonical values.
# C00 = 1 by normalisation convention; absent from most SHA files → added explicitly.
FALLBACK: dict[tuple[int,int], tuple[float,float]] = {
    (0, 0): ( 1.00000000000000e+00,  0.0),
    (1, 0): ( 0.0,                   0.0),
    (1, 1): ( 0.0,                   0.0),
    # Degree 2 — dominant oblateness / triaxiality terms
    (2, 0): (-9.08843393474243e-05,  0.0),            # J2 equivalent
    (2, 1): ( 1.46641235502819e-11,  1.17327642348892e-09),
    (2, 2): ( 3.46730964706963e-05,  9.07918983437229e-10),
    # Degree 3
    (3, 0): ( 2.56984255590099e-06,  0.0),
    (3, 1): ( 5.68099636038791e-07, -1.18987695809264e-06),
    (3, 2): ( 4.48789416551949e-07, -8.75440994898578e-08),
    (3, 3): ( 1.62927700547697e-07,  6.92580534018671e-08),
    # Degree 4
    (4, 0): (-1.58050897500699e-07,  0.0),
    (4, 1): (-2.30523478827299e-08,  1.07716462481040e-08),
    (4, 2): ( 2.35023671832479e-08, -3.50956434765018e-09),
    (4, 3): ( 1.22813892083349e-09, -6.73561025765649e-09),
    (4, 4): ( 5.39949297600469e-09, -1.09898804621459e-08),
}


def try_download() -> bool:
    if os.path.exists(CACHE) and os.path.getsize(CACHE) > 1_000_000:
        print("Using cached GRGM1200A file.", file=sys.stderr)
        return True
    for url in URLS:
        print(f"Downloading from {url} …", file=sys.stderr)
        try:
            urllib.request.urlretrieve(url, CACHE)
            if os.path.getsize(CACHE) > 1_000_000:
                print("Download OK.", file=sys.stderr)
                return True
            else:
                print("Downloaded file too small, trying next URL.", file=sys.stderr)
        except Exception as e:
            print(f"  Failed: {e}", file=sys.stderr)
    return False


def load_from_file(n_max: int) -> dict:
    """Parse SHA tab file with either comma or whitespace separators."""
    table = {}
    # Always add C00 = 1 (normalisation convention; absent from most SHA files)
    table[(0, 0)] = (1.0, 0.0)
    with open(CACHE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Handle both comma-separated and whitespace-separated formats
            parts = [p.strip() for p in line.replace(',', ' ').split()]
            if len(parts) < 4:
                continue
            try:
                n, m = int(parts[0]), int(parts[1])
                C, S = float(parts[2]), float(parts[3])
            except ValueError:
                continue
            if n > n_max:
                break
            table[(n, m)] = (C, S)
    print(f"Parsed {len(table)} coefficients up to degree {n_max}.", file=sys.stderr)
    return table


def build_table(n_max: int) -> dict:
    if try_download():
        tbl = load_from_file(n_max)
        if len(tbl) >= (n_max + 1) * (n_max + 2) // 2 - 10:
            return tbl
        print("Parsed table incomplete; supplementing with hardcoded values.", file=sys.stderr)
        # Merge: hardcoded fills any gaps
        for key, val in FALLBACK.items():
            if key not in tbl:
                tbl[key] = val
        return tbl
    print("All downloads failed. Using hardcoded fallback coefficients.", file=sys.stderr)
    # Pad zeros for degrees not in the fallback
    tbl: dict = {}
    for n in range(n_max + 1):
        for m in range(n + 1):
            tbl[(n, m)] = FALLBACK.get((n, m), (0.0, 0.0))
    return tbl


def write_binary(table: dict, n_max: int, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "wb") as f:
        f.write(struct.pack("<I", n_max))
        for n in range(n_max + 1):
            for m in range(n + 1):
                C, S = table.get((n, m), (0.0, 0.0))
                f.write(struct.pack("<dd", C, S))
    size = os.path.getsize(out_path)
    expected = 4 + (n_max + 1) * (n_max + 2) // 2 * 16
    print(f"Wrote {out_path} ({size} bytes, expected {expected})", file=sys.stderr)
    assert size == expected, f"Size mismatch! Got {size}, expected {expected}"


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(script_dir, "..", "data", "grgm1200a_100x100.bin")
    table = build_table(N_MAX)
    write_binary(table, N_MAX, os.path.abspath(out))
