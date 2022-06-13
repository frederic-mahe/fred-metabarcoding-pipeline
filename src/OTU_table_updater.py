#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
   update an OTU table with new taxonomic assignments
"""

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2022/01/31"
__version__ = "$Revision: 2.0"

import sys
import argparse


# *************************************************************************** #
#                                                                             #
#                                  Functions                                  #
#                                                                             #
# *************************************************************************** #

if __name__ == '__main__':
    """
    Parse arguments from command line.
    """
    parser = argparse.ArgumentParser(
        description="update an OTU table with new taxonomic assignments.")

    parser.add_argument("--old_otu_table",
                        dest="old_otu_table",
                        required=True,
                        help="old OTU table")

    parser.add_argument("--new_taxonomy",
                        dest="new_taxonomy_file",
                        required=True,
                        help="new taxonomic assignments")

    parser.add_argument("--new_otu_table",
                        dest="new_otu_table",
                        required=True,
                        help="new OTU table")

    ARGS = parser.parse_args()


def parse_taxonomy(new_taxonomy_file):
    """
    Map amplicons and taxonomic assignments.
    """
    separator = "\t"
    amplicons = [dict() for i in range(0, 256)]

    with open(new_taxonomy_file, "r") as new_taxonomy_data:
        print("PROGRESS: parsing taxonomy", file=sys.stderr)
        for line in new_taxonomy_data:
            amplicon, abundance, identity, taxonomy, references = \
                line.strip().split(separator)
            index = int(amplicon[0:2], 16)
            amplicons[index][amplicon] = (identity, taxonomy, references)

    return amplicons


def update_otu_table(old_otu_table, amplicons, new_otu_table):
    """
    Update taxonomy, identity and references.

    Column numbers are:
    1	OTU
    2	total
    3	cloud
    4	amplicon
    5	length
    6	abundance
    7	chimera
    8	spread
    9	quality
    10	sequence
    11	identity
    12	taxonomy
    13	references
    """
    separator = "\t"
    is_first_line = True
    with open(old_otu_table, "r") as old_otu_data:
        with open(new_otu_table, "w") as new_otu_file:
            print("PROGRESS: parsing and updating old OTU table",
                  file=sys.stderr)
            for line in old_otu_data:
                line = line.strip().split(separator)

                # header line is printed as-is
                if is_first_line:
                    is_first_line = False
                    print("\t".join(line), file=new_otu_file)
                    continue

                # update identity, taxonomy and references
                amplicon = line[3]
                index = int(amplicon[0:2], 16)
                line[10], line[11], line[12] = amplicons[index][amplicon]
                print("\t".join(line), file=new_otu_file)


def main():
    """
    update an OTU table with new taxonomic assignments.
    """
    # capture arguments
    old_otu_table = ARGS.old_otu_table
    new_taxonomy_file = ARGS.new_taxonomy_file
    new_otu_table = ARGS.new_otu_table

    # Parse the new taxomomic results
    amplicons = parse_taxonomy(new_taxonomy_file)

    # Sort by decreasing abundance (and alphabetical name? no)
    # not possible as-of-now, old table must be memoized first...

    # Parse the old OTU table and write a new one
    update_otu_table(old_otu_table, amplicons, new_otu_table)


# *************************************************************************** #
#                                                                             #
#                                     Body                                    #
#                                                                             #
# *************************************************************************** #

if __name__ == '__main__':

    main()

sys.exit(0)
