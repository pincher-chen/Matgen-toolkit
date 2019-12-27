#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re


def getFather(relation, x):
    if relation[x] == x:
        return x
    return getFather(relation, relation[x])

def readfile(filename, threshold):
    items = set()

    eles = []
    relation = {}

    file = open(filename, 'r')
    for line in file.readlines():
        tmp = re.split('\s+', line.strip())
        eles.append(tmp)
        relation[tmp[0]] = tmp[0]
        relation[tmp[1]] = tmp[1]
    file.close()
    
    for ele in eles:
        if float(ele[2]) < threshold:
            x = getFather(relation, ele[0])
            y = getFather(relation, ele[1])
            if x == y:
                continue
            else:
                relation[y] = x

    sim = {}
    print("Total file:", len(relation))
    for item in relation:
        # print(item, relation[item])
        if relation[item] == item:
            items.add(item)
            sim[item] = set()
    
    for item in relation:
        sim[getFather(relation, item)].add(item)

    # print(sim)

    return items, sim

def writeFile(output, items, threshold):
    file = open(output + "/" + "analysis_" + str(threshold) + ".txt", 'w')
    for item in items:
        file.write(item + "\n")
    file.close()

    # file = open(output + "/" + "similar_file_" + str(threshold) + ".txt", 'w')

def writeSim(output, sim, threshold):
    file = open(output + "/" + "similarity_" + str(threshold) + ".txt", 'w')
    for item in sim:
        file.write(str(sim[item]) + "\n")
    file.close()

def main():
    parser = argparse.ArgumentParser(description='fingerprint log analysis')
    parser.add_argument('filename', 
                        help='log file')
    parser.add_argument('threshold', 
                        help='partitioning threshold in similarity calculation',
                        type=float)
    args = parser.parse_args()

    filename = args.filename
    threshold = args.threshold

    items, sim = readfile(filename, threshold)
    # print(items)
    writeFile("./", items, threshold)

    writeSim("./", sim, threshold)

    print("Result: ", len(items))

if __name__ == "__main__":
    main()
