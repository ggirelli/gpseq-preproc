#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from glob import glob
import logging
import os
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import regex as re  # type: ignore
from rich.logging import RichHandler  # type: ignore
import sys
from typing import Dict

version = "0.0.1"

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)

parser = argparse.ArgumentParser(
    description="Assemble summary table after running rereprpr.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument("root", type=str, help="Path to root folder.")

parser.add_argument(
    "--version",
    action="version",
    version=f"{sys.argv[0]} v{version}",
)

args = parser.parse_args()

logging.info(f"Looking into '{args.root}'...")

patterns: Dict[str, str] = {}
patterns[
    "SNPsplit_unassignable"
] = "^([0-9]+) reads were unassignable \(([0-9]+\.[0-9]+%)\)"
patterns[
    "SNPsplit_genome1"
] = "^([0-9]+) reads were specific for genome 1 \(([0-9]+\.[0-9]+%)\)"
patterns[
    "SNPsplit_genome2"
] = "^([0-9]+) reads were specific for genome 2 \(([0-9]+\.[0-9]+%)\)"
patterns["SNPsplit_weird"] = (
    "^([0-9]+) reads did not contain one of the expected"
    + " bases at known SNP positions \(([0-9]+\.[0-9]+%)\)"
)
patterns[
    "SNPsplit_conflict"
] = "^([0-9]+) contained conflicting allele-specific SNPs \(([0-9]+\.[0-9]+%)\)"

dataframe = pd.DataFrame()

fastq_dir_path = os.path.join(args.root, "fastq")
assert os.path.isdir(fastq_dir_path)

fastq_list = glob(os.path.join(fastq_dir_path, "*"))
library_id_list = [os.path.basename(x).split("_")[0] for x in fastq_list]
dataframe["prep_run"] = np.repeat(os.path.basename(args.root), len(library_id_list))

logging.info(f"Found {len(library_id_list)} library IDs: {library_id_list}")
dataframe["library_id"] = library_id_list

logging.info("Reading filtered alignment counts...")
for library_iid in range(len(library_id_list)):
    library_id = library_id_list[library_iid]
    log_path = os.path.join(args.root, "mapping", f"{library_id}.clean_count.txt")
    assert os.path.isfile(log_path)
    with open(log_path) as LH:
        dataframe.loc[library_iid, "mapped"] = int(LH.readlines()[0].strip())

logging.info("Reading unassignable counts...")
for library_iid in range(len(library_id_list)):
    library_id = library_id_list[library_iid]
    log_path = os.path.join(
        args.root, "mapping", f"{library_id}.clean.SNPsplit_report.txt"
    )
    assert os.path.isfile(log_path)
    with open(log_path) as LH:
        matched = False
        for line in LH:
            match = re.match(patterns["SNPsplit_unassignable"], line)
            if match is None:
                continue
            matched = True
            dataframe.loc[library_iid, "unassigned"] = int(match.groups()[0])
            dataframe.loc[library_iid, "unassigned%"] = match.groups()[1]
        assert matched, f"missing unassignable count line [{library_id}]"

logging.info("Reading genome 1 counts...")
for library_iid in range(len(library_id_list)):
    library_id = library_id_list[library_iid]
    log_path = os.path.join(
        args.root, "mapping", f"{library_id}.clean.SNPsplit_report.txt"
    )
    assert os.path.isfile(log_path)
    with open(log_path) as LH:
        matched = False
        for line in LH:
            match = re.match(patterns["SNPsplit_genome1"], line)
            if match is None:
                continue
            matched = True
            dataframe.loc[library_iid, "genome1"] = int(match.groups()[0])
            dataframe.loc[library_iid, "genome1%"] = match.groups()[1]
        assert matched, f"missing genome 1 count line [{library_id}]"

logging.info("Reading genome 2 counts...")
for library_iid in range(len(library_id_list)):
    library_id = library_id_list[library_iid]
    log_path = os.path.join(
        args.root, "mapping", f"{library_id}.clean.SNPsplit_report.txt"
    )
    assert os.path.isfile(log_path)
    with open(log_path) as LH:
        matched = False
        for line in LH:
            match = re.match(patterns["SNPsplit_genome2"], line)
            if match is None:
                continue
            matched = True
            dataframe.loc[library_iid, "genome2"] = int(match.groups()[0])
            dataframe.loc[library_iid, "genome2%"] = match.groups()[1]
        assert matched, f"missing genome 2 count line [{library_id}]"

logging.info("Reading weird counts...")
for library_iid in range(len(library_id_list)):
    library_id = library_id_list[library_iid]
    log_path = os.path.join(
        args.root, "mapping", f"{library_id}.clean.SNPsplit_report.txt"
    )
    assert os.path.isfile(log_path)
    with open(log_path) as LH:
        matched = False
        for line in LH:
            match = re.match(patterns["SNPsplit_weird"], line)
            if match is None:
                continue
            matched = True
            dataframe.loc[library_iid, "weird"] = int(match.groups()[0])
            dataframe.loc[library_iid, "weird%"] = match.groups()[1]
        assert matched, f"missing weird count line [{library_id}]"

logging.info("Reading conflicting counts...")
for library_iid in range(len(library_id_list)):
    library_id = library_id_list[library_iid]
    log_path = os.path.join(
        args.root, "mapping", f"{library_id}.clean.SNPsplit_report.txt"
    )
    assert os.path.isfile(log_path)
    with open(log_path) as LH:
        matched = False
        for line in LH:
            match = re.match(patterns["SNPsplit_conflict"], line)
            if match is None:
                continue
            matched = True
            dataframe.loc[library_iid, "conflict"] = int(match.groups()[0])
            dataframe.loc[library_iid, "conflict%"] = match.groups()[1]
        assert matched, f"missing conflicting count line [{library_id}]"


dataframe.sort_values("library_id").to_csv(
    os.path.join(args.root, "summary_table-allelic.tsv"), sep="\t", index=False
)
