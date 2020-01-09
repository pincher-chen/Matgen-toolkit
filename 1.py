#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import shutil

def readfile(filename):
    file = open(filename, 'r')

    new = open("new.txt", 'w')
    for line in file.readlines():
        t1 = re.split(':', line.strip())
        t2 = re.split('\s+', t1[1])
        print(t2)
        if len(t2) > 1:
            new.write(line)
    new.close()
    file.close()


def main():
    parser = argparse.ArgumentParser(description='csd solvent')
    parser.add_argument('filename', 
                        help='log file')

    args = parser.parse_args()
    filename = args.filename

    readfile(filename)

if __name__ == "__main__":
    main()