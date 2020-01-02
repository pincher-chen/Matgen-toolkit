#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import os
import shutil


def read(input_path, output_path):
    files = []
    base_path = input_path
    
    for file in os.listdir(input_path):  
        file_path = os.path.join(input_path, file)  
        if os.path.splitext(file_path)[1]=='.cif':  
            files.append(file)

    # print(files)
    error_dir = output_path + "/error_file"

    tMap = {"unknown": []}
    for file in files:
        sType = parse_cif(os.path.join(base_path, file))

        if sType not in tMap and sType is not None:
            tMap[sType] = []
        
        if sType is None:
            tMap["unknown"].append(file)
        else:
            tMap[sType].append(file)
    return tMap

def write(output_path, tMap):
    file = open(output_path + "/" + "icsd_structure_type.txt", 'w')    
    for item in tMap:
        file.write(item + " : ")
        file.write(str(tMap[item]) + "\n")
    file.close()


def parse_cif(filename):
    print(filename)

    file = open(filename, 'r')
    # lines = file.read().decode('utf8').split('\n')
    
    for line in file.readlines():
        if "_chemical_name_structure_type" in line:
            sType = re.split('\s+', line.strip())[1]
            return sType
    return None

def main():
    parser = argparse.ArgumentParser(description='fingerprint log analysis')
    parser.add_argument('dir', 
                        help='icsd file path')
    parser.add_argument('output',
                        help="output path")
    args = parser.parse_args()

    input_path = args.dir
    output_path = args.output

    tMap = read(input_path, output_path)
    write(output_path, tMap)


if __name__ == "__main__":
    main()
