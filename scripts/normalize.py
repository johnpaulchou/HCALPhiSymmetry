#!/usr/bin/env python3

import argparse
import common
import itertools
import sys
from collections import defaultdict
from statistics import median



def read_file(filename):
    data = {}

    with open(filename, "r") as f:
        for line_number, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines or comments
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) != 6:
                print(f"Warning: Skipping malformed line {line_number}: {line}")
                continue

            try:
                subdet = int(parts[0])
                ieta   = int(parts[1])
                iphi   = int(parts[2])
                depth  = int(parts[3])
                corr   = float(parts[4])
                err    = float(parts[5])
            except ValueError as e:
                print(f"Warning: Type conversion error on line {line_number}: {e}")
                continue

            key = (subdet, ieta, iphi, depth)

            if key in data:
                print(f"Warning: Duplicate entry for key {key} on line {line_number}")

            data[key] = (corr, err)

    return data


def compute_medians(data):
    grouped = defaultdict(list)

    # Group corr values by (subdet, ieta, depth)
    for (subdet, ieta, iphi, depth), (corr, err) in data.items():
        group_key = (subdet, ieta, depth)
        grouped[group_key].append(corr)

    # Compute median per group
    medians = {}
    for group_key, corr_values in grouped.items():
        medians[group_key] = median(corr_values)

    return medians


def main():

    # setup parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputfilename", help="name of file containing the corrections")
    parser.add_argument("outputfilename", help="name of the output file that contains missing channels and is normalized")
    args=parser.parse_args()

    # read the data
    data = read_file(args.inputfilename)
    print(f"Read {len(data)} valid entries")

    # loop over all possible keys and add the missing ones to the data (with corr=1.0 and err=0.0)
    for subdetindex,subdetnum in enumerate(common.subdetnums):
        for depth in range(1, common.ndepths[subdetindex]+1):
            ietas1 = range(common.minabsietas[subdetindex],common.maxabsietas[subdetindex]+1)
            ietas2 = range(-common.maxabsietas[subdetindex],-common.minabsietas[subdetindex]+1)
            ietas = itertools.chain(ietas1, ietas2)
            for ieta in ietas:
                for iphi in range(common.miniphi, common.maxiphi+1):

                    val=common.goodChannel(subdetnum, ieta, iphi, depth, 0)
                    key=(subdetnum, ieta, iphi, depth)

                    if key not in data and val:
                        print("Missing key="+str(key))
                        data[key]=(1.0,0.0)

    # check that the length is correct now
    if len(data) != 17424:
        print(f"Warning: Expected 17424 entries, found {len(data)}")

    # compute the medians
    medians = compute_medians(data)
    print(f"Computed medians for {len(medians)} (subdet, ieta, depth) groups")
    for (subdet, ieta, depth), med in sorted(medians.items()):
        print(f"{subdet:3d} {ieta:4d} {depth:2d} {med:12.6f}")

    # Reprint normalized data
    with open(args.outputfilename, "w") as out:
        for (subdet, ieta, iphi, depth), (corr, err) in data.items():
            group_key = (subdet, ieta, depth)
            med = medians[group_key]

            if med == 0:
                print(f"Warning: Zero median for group {group_key}")
                continue

            corr_new = corr / med
            err_new  = err  / med

            out.write(f"{subdet} {ieta} {iphi} {depth} "
                      f"{corr_new:.6f} {err_new:.6f}\n")


    
if __name__ == "__main__":
    main()
