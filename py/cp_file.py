#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import shutil

def readfile(filename, output):
    file = open(filename, 'r')
    for line in file.readlines():
        tmp = re.split('\s+', line.strip())
        shutil.copy(tmp[0], output)


def main():
    parser = argparse.ArgumentParser(description='csd solvent')
    parser.add_argument('filename', 
                        help='log file')
    parser.add_argument('output', 
                        help='output dir') 

    args = parser.parse_args()
    filename = args.filename
    output = args.output

    readfile(filename, output)

if __name__ == "__main__":
    main()