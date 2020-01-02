#!/bin/env python
import pymongo
import chardet
import argparse
import re
import os
import shutil
import pathlib

files = []
sim = {}
sim_file = ""


def read_dir(path):
    for item in os.listdir(path):
        cur = os.path.join(path, item)
        if os.path.isdir(cur):
            read_dir(cur)
        else:
            if os.path.splitext(cur)[1] == '.cif':
                files.append(cur)
            else:
                global sim_file
                sim_file = cur


def get_clean(str):
    s = str.strip()

    if len(s) == 0:
        return s
    
    if s[0] == '"' or s[0] == '\'':
        s = s[1:]
    if s[-1] == '"' or s[-1] == '\'':
        s = s[:-1]
    return s

def nameModified(name):
    t = name.split('-')
    return t[0] + ".cif"

def data2db(db, input_path, output_path):

    print("Get All File")
    read_dir(input_path)

    print("Get Similar File Log")
    if sim_file == "":
        print("Similar file log Not found!")
        return

    sim_in = open(sim_file, 'r')
    for line in sim_in.readlines():
        items = re.split(':', line)
        sfile = re.split('\s+', items[1].strip())
        sfile = [nameModified(x) for x in sfile]
        sim[items[0].strip()] = " ".join(sfile)
    # print(sim)

    error_file = output_path + "/mongodb_icsd_error_log.txt"
    err_out = open(error_file, 'w')
    
    print("ICSD Data:")
    for file in files:
        # data = parse_cif(file)
        # print(data)
        print(file)
        try:
            data = parse_cif(file)
            # print(data)
            db.insert_one(data)
        except BaseException as e:
            # print(file, e)
            err_out.write(file + ":" + str(e) + "\n")
        finally:
            pass

    err_out.close()

def parse_cif(filename):
    # print(os.path.split(filename))
    paths = pathlib.Path(filename).parts
    data = {
        "icsd_code": "",
        "name_systematic": "",
        "formula_structural": "",
        "formula_sum": "",
        "name_structure_type": "",

        "element_type": "",
        "space_group": "",
        "similar_file": "",
        "cif_data": ""
    }


    f = open(filename, 'rb')
    encoding = chardet.detect(f.read()).get('encoding')
    f.close()

    file = open(filename, mode='r', encoding=encoding)

    buff = ""
    tag = ""

    for line in file.readlines():
        buff += line
        if len(tag) > 0:
            data[tag] = get_clean(line)
            tag = ""
        if "code_ICSD" in line:
            if len(list(filter(None, re.split('\s+', line)))) > 1:
                data["icsd_code"] = get_clean(" ".join(re.split('\s+', line)[1:]))
            else:
                tag="icsd_code"
        if "name_systematic" in line:
            if len(list(filter(None, re.split('\s+', line)))) > 1:
                data["name_systematic"] = get_clean(" ".join(re.split('\s+', line)[1:]))
            else:
                tag += "name_systematic"
        if "formula_structural" in line:
            if len(list(filter(None, re.split('\s+', line)))) > 1:
                data["formula_structural"] = get_clean(" ".join(re.split('\s+', line)[1:]))
            else:
                tag="formula_structural"
        if "formula_sum" in line:
            if len(list(filter(None, re.split('\s+', line)))) > 1:
                data["formula_sum"] = get_clean(" ".join(re.split('\s+', line)[1:]))
            else:
                tag = "formula_sum"
        if "_name_structure_type" in line:
            if len(list(filter(None, re.split('\s+', line)))) > 1:
                data["name_structure_type"] = get_clean(" ".join(re.split('\s+', line)[1:]))
            else:
                tag = "name_structure_type"
    file.close()

    if paths[-1] in sim:
        data["similar_file"] = sim[paths[-1]]
    else:
        data["similar_file"] = ""

    data["element_type"] = paths[-3]
    data["space_group"] = paths[-4]
    data["cif_data"] = buff
    return data


def main():
    # mongodb
    mongodb_addr = "12.11.70.12:10101"
    db_name = "ucsd"
    col_name = "icsd"

    myclient = pymongo.MongoClient(mongodb_addr)
    mydb = myclient[db_name]
    mycol = mydb[col_name]

    parser = argparse.ArgumentParser(description='MongoDB - icsd')
    parser.add_argument('dir',
                        help='icsd file path')
    parser.add_argument('output',
                        help="error file output path")
    args = parser.parse_args()

    input_path = args.dir
    output_path = args.output

    data2db(mycol, input_path, output_path)


if __name__ == "__main__":
    main()
