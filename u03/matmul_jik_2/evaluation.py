#!/usr/bin/python
# author: David Scherfgen
# This script expects the Dinero output files to be in subdirectories.
# The first hierarchy denotes the problem size (like 16, 128 and 1024).
# In each subdirectory, there must be another subdirectory named like
# the aspect that is being varied (like l1assoc, l2size, linesize, ...).
# In these subdirectories, the script looks for a text file named after
# the concrete value that was used for the selected aspect.
#
# For example, the Dinero results for a run with
# - problem size = 1024
# - L2 associativity = 2m
#
# the results should be in:
# ./1024/l2assoc/2m.txt
#
# The output is a CSV file that can be opened with Excel or Open Office.

import re
import sys

# through command line parameters (see below)
#problems = ("100", "1000")

aspects = {
           "l1dsize"   : ("L1 Size", ("1k", "2k", "4k", "8k", "16k", "64k", "128k", "256k", "512k")),
           "l2usize"   : ("L2 Size", ("64k", "128k", "256k", "512k", "1m", "2m", "4m", "8m", "16m")),
           "l3usize"   : ("L3 Size", ("512k", "1m", "2m", "4m", "8m", "16m", "32m", "64m", "128m")),
           "l1dassoc"  : ("L1 Associativity", ("1", "2", "3", "4", "5", "6", "7", "8")),
           "l2uassoc"  : ("L2 Associativity", ("1", "2", "3", "4", "8", "16", "32", "64", "128")),
           "l3uassoc"  : ("L3 Associativity", ("1", "2", "3", "4", "8", "16", "32", "64", "128")),
           "lsize"     : ("Line Size", ("8", "16", "32", "64", "128", "256", "512", "1024"))
           }

# costs for misses (in cycles)
costL1=2
costL2=6
costL3=29
costMM=98

# matches lines like:
# Demand Fetches		     4333426
demand_fetches_re = re.compile(r"Demand Fetches(\s*)(?P<value>\d+)")

# matches lines like:
# Demand Misses		      744276
demand_misses_re = re.compile(r"Demand Misses(\s*)(?P<value>\d+)")

if len(sys.argv) < 2:
    print "error: no program arguments given\nusage: python evaluation.py arg_1 ... arg_n\n\twhere arg_i is the i-th program argument\n\te.g.python evaluation.py 100 200 300"
    sys.exit(1)

problems = sys.argv[1:]

for aspect in aspects.items():
    aspect_id = aspect[0]
    aspect_name = aspect[1][0]
    aspect_options = aspect[1][1]
    
    f = open(aspect_id + ".csv", "w")
    f.write(aspect_name + "\n")
    f.write("problem;" + aspect_id + ";l1access;l1miss;l2access;l2miss;l3access;l3miss;linesize;cycles;relcycles\n")
    
    for problem in problems:
        ref_cycles = 0
        for aspect_option in aspect_options:
            f.write(problem + ";" + aspect_option + ";")
            g = open(problem + "/" + aspect_id + "/" + aspect_option + ".txt", "r")
            values = []
            for line in g:
                result = re.search(demand_fetches_re, line)
                value = 0
                if result:
                    value = int(result.group("value"))
                    f.write(str(value) + ";")
                else:
                    result = re.search(demand_misses_re, line)
                    if result:
                        value = int(result.group("value"))
                        f.write(str(value) + ";")
                if value != 0: values.append(value)
            g.close()

            # do cost calculation
            line_size = 64
            if aspect_id == "lsize": line_size = int(aspect_option)
            f.write(str(line_size) + ";")
            linesfactor = (line_size // 64 - 1)
            line_penalty = [costL2,costL3,costMM]
            if(linesfactor < 0): linesfactor = 0
            # values[0], values[1]: demand fetches / demand misses L1
            # values[2], values[3]: demand fetches / demand misses L2
            # values[4], values[5]: demand fetches / demand misses L3
            if len(values) > 0:
                cycles  = (values[0] - values[1]) * costL1 # L1 hits
                cycles += (values[2] - values[3]) * costL2 + (values[1] * linesfactor * line_penalty[0]) # L2 hits
                cycles += (values[4] - values[5]) * costL3 + (values[3] * linesfactor * line_penalty[1]) # L3 hits
                cycles += values[5]               * costMM + (values[5] * linesfactor * line_penalty[2]) # main memory
                if ref_cycles == 0: ref_cycles = cycles
                f.write(str(cycles) + ";")
                rel_cycles = float(cycles) / ref_cycles
                f.write(str(rel_cycles).replace(".", ","))
            f.write("\n")
    f.close()
