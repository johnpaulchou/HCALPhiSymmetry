#!/usr/bin/env python3

import argparse
import array
import json
import math

import ROOT


def load_hf_corrections(filename):
    corrections = {}

    with open(filename) as corr_file:
        for line_number, line in enumerate(corr_file, start=1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue

            fields = line.split()
            if len(fields) != 6:
                raise ValueError(
                    "{}:{}: expected 6 fields, found {}".format(
                        filename, line_number, len(fields)
                    )
                )

            try:
                subdet = int(fields[0])
            except ValueError as error:
                raise ValueError(
                    "{}:{}: invalid subdetector '{}'".format(
                        filename, line_number, fields[0]
                    )
                ) from error

            if subdet != 4:
                continue

            try:
                key = (int(fields[1]), int(fields[2]), int(fields[3]))
                corr = float(fields[4])
                float(fields[5])  # correrr is not needed for the MET calculation
            except ValueError as error:
                raise ValueError(
                    "{}:{}: invalid HF correction row".format(
                        filename, line_number
                    )
                ) from error

            if not math.isfinite(corr):
                raise ValueError(
                    "{}:{}: correction is not finite".format(filename, line_number)
                )
            if key in corrections:
                raise ValueError(
                    "{}:{}: duplicate HF correction for {}".format(
                        filename, line_number, key
                    )
                )
            corrections[key] = corr

    return corrections


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputfns", nargs="+", help="name(s) of .root file(s) to be processed")
    parser.add_argument("--corrfile", help="file with the HF response corrections")
    parser.add_argument("-o", "--outputfn", default="hfmet.root", help="name of the .root file that stores the histograms")
    parser.add_argument("--hf-geometry-json", default="hf_geometry.json", dest="hfgeofile", help='.json file produced by the "cmsRun printHF.py" command')
    parser.add_argument("--numevents", type=int, default=-1, help="number of events to process")
    
    args = parser.parse_args()

    # Store the cell center and its direction from (0,0,0). The direction from
    # the PV must be recalculated event by event.
    with open(args.hfgeofile) as geometry_file:
        hfdata = json.load(geometry_file)

    hf_geometry = {}
    for cell in hfdata["hf_cells"]:
        key = (cell["ieta"], cell["iphi"], cell["depth"])
        x = float(cell["x_cm"])
        y = float(cell["y_cm"])
        z = float(cell["z_cm"])
        distance = math.sqrt(x*x + y*y + z*z)
        if distance == 0.0:
            raise ValueError("HF cell {} is located at the origin".format(key))
        hf_geometry[key] = (x, y, z, x/distance, y/distance)

    hf_corrections = None
    if args.corrfile is not None:
        hf_corrections = load_hf_corrections(args.corrfile)
        print("Loaded HF corrections:", len(hf_corrections))

    chain = ROOT.TChain("phisym/phisymtree")
    for filename in args.inputfns:
        chain.Add(filename)

    total_entries = chain.GetEntries()
    nentries = (
        total_entries
        if args.numevents < 0
        else min(args.numevents, total_entries)
    )
    print("Total entries to loop over:", nentries)

    event = array.array("L", [0])
    nHits = array.array("i", [0])
    max_nhits = 20000
    hitEnergy = array.array("f", [0]*max_nhits)
    ieta = array.array("i", [0]*max_nhits)
    iphi = array.array("i", [0]*max_nhits)
    depth = array.array("i", [0]*max_nhits)
    subdet = array.array("i", [0]*max_nhits)
    pvx = array.array("f", [0])
    pvy = array.array("f", [0])
    pvz = array.array("f", [0])

    main_branches = [
        "event",
        "nHits",
        "hitEnergy",
        "ieta",
        "iphi",
        "depth",
        "subdet",
        "pvx",
        "pvy",
        "pvz",
    ]
    chain.SetBranchStatus("*", 0)
    for branch_name in main_branches:
        chain.SetBranchStatus(branch_name, 1)

    chain.SetBranchAddress("event", event)
    chain.SetBranchAddress("nHits", nHits)
    chain.SetBranchAddress("hitEnergy", hitEnergy)
    chain.SetBranchAddress("ieta", ieta)
    chain.SetBranchAddress("iphi", iphi)
    chain.SetBranchAddress("depth", depth)
    chain.SetBranchAddress("subdet", subdet)
    chain.SetBranchAddress("pvx", pvx)
    chain.SetBranchAddress("pvy", pvy)
    chain.SetBranchAddress("pvz", pvz)

    hMETwrtPV = ROOT.TH1D("hMETwrtPV", "HF MET w.r.t. PV;HF MET [GeV];Events", 100, 0, 100)
    hMETXwrtPV = ROOT.TH1D("hMETXwrtPV", "HF MET x-component w.r.t. PV;HF MET_{x} [GeV];Events",100,-100,100)
    hMETYwrtPV = ROOT.TH1D("hMETYwrtPV", "HF MET y-component w.r.t. PV;HF MET_{y} [GeV];Events", 100, -100, 100)
    hMETwrt000 = ROOT.TH1D("hMETwrt000", "HF MET w.r.t. (0,0,0);HF MET [GeV];Events", 100, 0, 100)
    hMETXwrt000 = ROOT.TH1D("hMETXwrt000", "HF MET x-component w.r.t. (0,0,0);HF MET_{x} [GeV];Events", 100, -100, 100)
    hMETYwrt000 = ROOT.TH1D("hMETYwrt000", "HF MET y-component w.r.t. (0,0,0);HF MET_{y} [GeV];Events", 100, -100, 100)

    processed_events = 0

    for ievt in range(nentries):
        if ievt % 10000 == 0:
            print("Event:", ievt)
        chain.GetEntry(ievt)

        # phiSymTree stores -999 in all three coordinates when no quality PV
        # is found. Skip those events so both MET definitions use one sample.
        if pvx[0] < -900.0 or pvy[0] < -900.0 or pvz[0] < -900.0:
            continue

        sum_px_pv = 0.0
        sum_py_pv = 0.0
        sum_px_000 = 0.0
        sum_py_000 = 0.0

        for ihit in range(nHits[0]):
            if subdet[ihit] != 4:
                continue

            cell_key = (ieta[ihit], iphi[ihit], depth[ihit])
            geometry = hf_geometry.get(cell_key)
            if geometry is None:
                raise KeyError(
                    "No HF geometry for cell {} in event {}".format(
                        cell_key, event[0]
                    )
                )

            x, y, z, unit_x_000, unit_y_000 = geometry
            energy = hitEnergy[ihit]
            if hf_corrections is not None:
                corr = hf_corrections.get(cell_key)
                if corr is None:
                    raise KeyError(
                        "No HF correction for cell {} in event {}".format(
                            cell_key, event[0]
                        )
                    )
                energy *= corr

            # Treat the calorimeter deposit as massless and point it from the
            # selected reference vertex toward the HF cell center.
            sum_px_000 += energy*unit_x_000
            sum_py_000 += energy*unit_y_000

            dx = x - pvx[0]
            dy = y - pvy[0]
            dz = z - pvz[0]
            distance_from_pv = math.sqrt(dx*dx + dy*dy + dz*dz)
            sum_px_pv += energy*dx/distance_from_pv
            sum_py_pv += energy*dy/distance_from_pv

        # Missing transverse momentum is the negative vector sum of the
        # reconstructed transverse momenta.
        met_x_pv = -sum_px_pv
        met_y_pv = -sum_py_pv
        met_x_000 = -sum_px_000
        met_y_000 = -sum_py_000

        hMETXwrtPV.Fill(met_x_pv)
        hMETYwrtPV.Fill(met_y_pv)
        hMETwrtPV.Fill(math.hypot(met_x_pv, met_y_pv))
        hMETXwrt000.Fill(met_x_000)
        hMETYwrt000.Fill(met_y_000)
        hMETwrt000.Fill(math.hypot(met_x_000, met_y_000))
        processed_events += 1

    output_file = ROOT.TFile(args.outputfn, "RECREATE")
    if not output_file or output_file.IsZombie():
        raise OSError("Could not create output file {}".format(args.outputfn))

    for hist in (
        hMETwrtPV,
        hMETXwrtPV,
        hMETYwrtPV,
        hMETwrt000,
        hMETXwrt000,
        hMETYwrt000,
    ):
        hist.Write()
    output_file.Close()

    print("Events filled:", processed_events)
    print("Wrote:", args.outputfn)


if __name__ == "__main__":
    main()
