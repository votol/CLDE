#!/usr/bin/env python3

import argparse
import os
import errno

parser = argparse.ArgumentParser(description='Generation parameters')
parser.add_argument('-i','--input', required=True, help='input cl file')
parser.add_argument('-o','--output', required=True, help='output header file')
args = parser.parse_args()

try:
    with open(args.input, "r") as infile:
        content = infile.read()
    if not os.path.exists(os.path.dirname(args.output)):
        try:
            os.makedirs(os.path.dirname(args.output))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    content = content.replace("\n","\\n\\\n")
    content = content.replace('"', '\\"')
    with open(args.output,"w") as outfile:
        outfile.write("#pragma once\n")
        variable_name = os.path.basename(args.output)
        variable_name = variable_name[:variable_name.find('.')]
        outfile.write("const char* ")
        outfile.write(variable_name + "_kernel")
        outfile.write(' = "')
        outfile.write(content)
        outfile.write('";\n')
except RuntimeError as e:
    print(str(e))
    if os.path.isfile(args.output):
        os.remove(args.output)
    exit(-1)
