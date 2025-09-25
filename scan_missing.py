#!/usr/bin/env python3
"""
Find & insert missing HCAL channels using im.goodChannel(), run sanity checks,
print a summary, and report duplicates.

- Prints ONLY missing lines as: "col1 col2 col3 col4 1.00000 0.00000"
  (sorted by col1, col2, col3, col4).
- Optionally writes a verified "filled" file that contains original + inserted lines
  in sorted order, but ONLY if all sanity checks pass.
- Provides a compact summary (totals by subdetector, triplets, quads, φ multiplicities).
- Reports duplicate text lines and duplicate quadruplets with counts.

Usage:
  python scan_missing.py input.txt > missing.txt
  python scan_missing.py input.txt --sanity-only filled.txt
  python scan_missing.py input.txt --write-filled filled.txt
"""

import argparse
import sys
from collections import Counter
from im import goodChannel
import re
# Subdetector mapping and expected totals
SUBDET_MAP = {1: "HB", 2: "HE", 4: "HF"}
VALID_SUBDET_IDS = set(SUBDET_MAP.keys())

EXPECTED_TOTAL_LINES = 17424
EXPECTED_UNIQUE_LINES = 17424
EXPECTED_TRIPLETS = 324
EXPECTED_QUADS = 17424

def parse_args():
    p = argparse.ArgumentParser(
        description="Insert missing HCAL channels using im.goodChannel(), sanity-check, summarize, and report duplicates."
    )
    p.add_argument("input", help="Input text file (whitespace-separated; at least 4 cols).")
    p.add_argument("--strict", action="store_true", help="Error on malformed lines (default: skip).")
    p.add_argument("--write-filled", metavar="OUTPUT.txt",
                   help="Write verified file (original + inserted), sorted by (col1,col2,col3,col4).")
    p.add_argument("--sanity-only", metavar="FILE.txt",
                   help="Run sanity checks + summary + duplicate report on an existing file and exit.")
    return p.parse_args()

def read_key_line_pairs(path, strict=False):
    """
    Read lines with ≥4 columns (skip blank/comment lines).
    Cleans extraneous whitespace/tabs and skips empty lines.
    Returns:
      keys: set[(c1, ieta, iphi, depth)]
      lines_by_key: dict key -> original cleaned line text (first occurrence kept)
      raw_lines: list[str] of all cleaned non-comment, non-blank lines
      key_counter: Counter of quadruplet keys (duplicates show count > 1)
      line_counter: Counter of exact text lines (duplicates show count > 1)
    """
    keys = set()
    lines_by_key = {}
    raw_lines = []
    key_counter = Counter()
    line_counter = Counter()
    with open(path, "r", encoding="latin-1") as f:
        for ln, raw in enumerate(f, 1):
            # normalize whitespace (convert tabs -> spaces, collapse multiple spaces)
            s = raw.strip()
            s = re.sub(r"\s+", " ", s)
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 4:
                if strict:
                    raise ValueError(f"Line {ln}: expected ≥4 columns, got {len(parts)} -> {raw!r}")
                else:
                    continue
            try:
                c1 = int(parts[0]); ieta = int(parts[1]); iphi = int(parts[2]); depth = int(parts[3])
            except ValueError:
                if strict:
                    raise
                else:
                    continue
            key = (c1, ieta, iphi, depth)
            raw_lines.append(s)
            line_counter[s] += 1
            key_counter[key] += 1
            keys.add(key)
            lines_by_key.setdefault(key, s)
    return keys, lines_by_key, raw_lines, key_counter, line_counter

def build_expected_keys():
    """
    Build the full set of expected (c1, ieta, iphi, depth) keys from geometry via im.goodChannel().
    We call goodChannel(subdet, ieta, iphi, depth, mod=1).
    """
    expected = set()
    for c1, subdet in SUBDET_MAP.items():
        # ieta in [-41..-1, 1..41]
        for ieta in list(range(-41, 0)) + list(range(1, 42)):
            # iphi in [1..72]
            for iphi in range(1, 73):
                # depth in [1..7] (goodChannel will reject invalid ones)
                for depth in range(1, 8):
                    if goodChannel(subdet, ieta, iphi, depth, 1):
                        expected.add((c1, ieta, iphi, depth))
    return expected

def expected_phi_set_for_triplet(c1, ieta, depth):
    """
    Return the set of iphi (1..72) expected by geometry for a given (det,ieta,depth),
    computed via im.goodChannel(subdet, ieta, iphi, depth, 1).
    """
    subdet = SUBDET_MAP.get(c1)
    return {iphi for iphi in range(1, 73) if goodChannel(subdet, ieta, iphi, depth, 1)}

def sanity_checks_from_keys(all_lines_count, keys_set, lines_by_key,
                            key_counter=None, line_counter=None):
    """
    Run sanity checks on a *candidate* dataset represented by keys_set.
    If counters are provided, also report duplicates.
    Returns (ok:bool, messages:list[str])
    """
    msgs = []
    ok = True

    uniq_lines_count = len(lines_by_key)
    unique_quads = len(keys_set)
    triplets = {(c1, ieta, depth) for (c1, ieta, iphi, depth) in keys_set}
    unique_trips = len(triplets)

    # 1) total lines
    if all_lines_count != EXPECTED_TOTAL_LINES:
        ok = False
        msgs.append(f"[FAIL] Total lines = {all_lines_count}, expected {EXPECTED_TOTAL_LINES}")
    else:
        msgs.append(f"[OK]   Total lines = {all_lines_count}")

    # 2) unique lines
    if uniq_lines_count != EXPECTED_UNIQUE_LINES:
        ok = False
        msgs.append(f"[FAIL] Unique lines = {uniq_lines_count}, expected {EXPECTED_UNIQUE_LINES}")
    else:
        msgs.append(f"[OK]   Unique lines = {uniq_lines_count}")

    # 3) unique triplets
    if unique_trips != EXPECTED_TRIPLETS:
        ok = False
        msgs.append(f"[FAIL] Unique (det,ieta,depth) triplets = {unique_trips}, expected {EXPECTED_TRIPLETS}")
    else:
        msgs.append(f"[OK]   Unique (det,ieta,depth) triplets = {unique_trips}")

    # 4) unique quadruplets
    if unique_quads != EXPECTED_QUADS:
        ok = False
        msgs.append(f"[FAIL] Unique (det,ieta,depth,iphi) quadruplets = {unique_quads}, expected {EXPECTED_QUADS}")
    else:
        msgs.append(f"[OK]   Unique (det,ieta,depth,iphi) quadruplets = {unique_quads}")

    # 5) per-triplet φ multiplicity: 72/36/18 and matches geometry
    per_trip_fail_examples = []
    for (c1, ieta, depth) in sorted(triplets):
        have_phis = {iphi for (d, e, iphi, dep) in keys_set if (d, e, dep) == (c1, ieta, depth)}
        expected_phis = expected_phi_set_for_triplet(c1, ieta, depth)
        n = len(have_phis)
        if n not in (72, 36, 18):
            ok = False
            per_trip_fail_examples.append(
                f"(det={c1}, ieta={ieta}, depth={depth}): have {n} φs (allowed: 72/36/18)"
            )
        if have_phis != expected_phis:
            ok = False
            missing = sorted(expected_phis - have_phis)
            extras  = sorted(have_phis - expected_phis)
            per_trip_fail_examples.append(
                f"(det={c1}, ieta={ieta}, depth={depth}): φ mismatch; "
                f"missing={missing[:6]}{'...' if len(missing)>6 else ''}, "
                f"extras={extras[:6]}{'...' if len(extras)>6 else ''}"
            )
            if len(per_trip_fail_examples) >= 5:
                break
    if per_trip_fail_examples:
        msgs.append("[FAIL] Per-triplet φ counts/geometry mismatches:")
        msgs.extend(f"       - {m}" for m in per_trip_fail_examples)
    else:
        msgs.append("[OK]   Per-triplet φ multiplicities match geometry (all are 72/36/18 as expected)")

    # 6) duplicate reporting (if counters given)
    if line_counter is not None:
        dup_text = [(line, cnt) for line, cnt in line_counter.items() if cnt > 1]
        dup_text.sort(key=lambda x: -x[1])
        if dup_text:
            msgs.append(f"[WARN] Duplicate TEXT lines: {len(dup_text)} unique texts duplicated")
            for line, cnt in dup_text[:20]:
                msgs.append(f"       x{cnt} : {line}")
            if len(dup_text) > 20:
                msgs.append(f"       ... and {len(dup_text)-20} more")
        else:
            msgs.append("[OK]   No duplicate TEXT lines")

    if key_counter is not None:
        dup_keys = [(key, cnt) for key, cnt in key_counter.items() if cnt > 1]
        dup_keys.sort(key=lambda x: -x[1])
        if dup_keys:
            msgs.append(f"[WARN] Duplicate QUADRUPLETS: {len(dup_keys)} duplicated channel keys")
            for (c1, ieta, iphi, depth), cnt in dup_keys[:20]:
                msgs.append(f"       x{cnt} : ({c1} {ieta} {iphi} {depth})")
            if len(dup_keys) > 20:
                msgs.append(f"       ... and {len(dup_keys)-20} more")
        else:
            msgs.append("[OK]   No duplicate QUADRUPLETS")

    return ok, msgs

def print_summary(keys_set, stream=sys.stderr):
    """Compact summary: totals per subdetector & φ multiplicity histogram."""
    triplets = {(c1, ieta, depth) for (c1, ieta, iphi, depth) in keys_set}
    quads_by_subdet = Counter()
    trips_by_subdet = Counter()
    for (c1, ieta, iphi, depth) in keys_set:
        quads_by_subdet[SUBDET_MAP[c1]] += 1
    for (c1, ieta, depth) in triplets:
        trips_by_subdet[SUBDET_MAP[c1]] += 1

    mult_hist = Counter()
    for (c1, ieta, depth) in triplets:
        nphi = sum(1 for (d, e, iphi, dep) in keys_set if (d, e, dep) == (c1, ieta, depth))
        mult_hist[nphi] += 1

    print("\n=== SUMMARY ===", file=stream)
    print(f"Total quadruplets: {len(keys_set)}", file=stream)
    print(f"Total triplets:    {len(triplets)}", file=stream)
    for sd in ("HB", "HE", "HF"):
        print(f"  {sd}: triplets={trips_by_subdet[sd]:4d}  quads={quads_by_subdet[sd]:5d}", file=stream)
    print("φ multiplicity (triplet counts):", file=stream)
    for k in sorted(mult_hist.keys(), reverse=True):  # usually 72,36,18
        print(f"  {k:>2} φ : {mult_hist[k]} triplets", file=stream)
    print("===============\n", file=stream)

def write_filled_file(output_path, all_keys_sorted, lines_by_key):
    with open(output_path, "w", encoding="latin-1") as out:
        for key in all_keys_sorted:
            if key in lines_by_key:
                out.write(lines_by_key[key] + "\n")
            else:
                c1, ieta, iphi, depth = key
                out.write(f"{c1} {ieta} {iphi} {depth} 1.00000 0.00000\n")

def main():
    args = parse_args()

    # Mode 1: sanity check an existing file and exit (with summary + duplicate report)
    if args.sanity_only:
        keys, lines_by_key, raw_lines, key_counter, line_counter = read_key_line_pairs(
            args.sanity_only, strict=args.strict
        )
        ok, msgs = sanity_checks_from_keys(
            all_lines_count=len(raw_lines),
            keys_set=keys,
            lines_by_key=lines_by_key,
            key_counter=key_counter,
            line_counter=line_counter,
        )
        print("\n".join(msgs))
        print_summary(keys, stream=sys.stdout)
        sys.exit(0 if ok else 1)

    # Normal flow: read input, compute expected, identify missing, print missing
    present_keys, present_lines_by_key, present_raw_lines, key_counter, line_counter = read_key_line_pairs(
        args.input, strict=args.strict
    )
    expected_keys = build_expected_keys()

    missing = sorted(expected_keys - present_keys, key=lambda t: (t[0], t[1], t[2], t[3]))

    # Always print missing lines to stdout
    for c1, ieta, iphi, depth in missing:
        print(f"{c1} {ieta} {iphi} {depth} 1.00000 0.00000")

    # If writing a filled file, assemble candidate and VERIFY before writing (also show summary + dups)
    if args.write_filled:
        candidate_keys = expected_keys  # filled file contains full expected geometry

        # Build candidate "lines_by_key": present as-is + placeholders for missing
        candidate_lines_by_key = dict(present_lines_by_key)
        for key in (candidate_keys - present_keys):
            c1, ieta, iphi, depth = key
            candidate_lines_by_key[key] = f"{c1} {ieta} {iphi} {depth} 1.00000 0.00000"

        # For the candidate dataset, counts: every key will appear exactly once
        candidate_total_lines = len(candidate_keys)

        ok, msgs = sanity_checks_from_keys(
            all_lines_count=candidate_total_lines,
            keys_set=candidate_keys,
            lines_by_key=candidate_lines_by_key,
            # constructed candidate has no duplicates by design; omit counters here
        )
        print("\n".join(msgs), file=sys.stderr)
        print_summary(candidate_keys, stream=sys.stderr)

        if not ok:
            sys.stderr.write("[ABORT] Sanity checks failed. Not writing --write-filled output.\n")
            sys.exit(1)

        all_keys_sorted = sorted(candidate_keys, key=lambda t: (t[0], t[1], t[2], t[3]))
        write_filled_file(args.write_filled, all_keys_sorted, present_lines_by_key)
        print(f"[OK]   Wrote verified file: {args.write_filled}", file=sys.stderr)

if __name__ == "__main__":
    main()
