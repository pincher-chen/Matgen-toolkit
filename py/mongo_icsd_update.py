#!/bin/env python
import pymongo
import chardet
import argparse
import re
import os
import shutil
import pathlib

def nameModified(name):
    t = name.split('-')
    return t[0] + ".cif"

def data2db(db, output_path):
    error_file = output_path + "/mongodb_icsd_error_log.txt"
    err_out = open(error_file, 'w')

    print("ICSD Data:")
    files = db.find()
    for file in files:
        print(file["similar_file"])
        try:
            if file["similar_file"] == []:
                db.update({'_id': file["_id"]}, {'$set': {"similar_file": ""}})
            else:
                sfile = [nameModified(x) for x in file["similar_file"].split(" ")]
                fstr = " ".join(sfile)
                db.update({'_id': file["_id"]}, {'$set': {"similar_file": fstr}})
        except BaseException as e:
            # print(file, e)
            err_out.write(file["icsd_code"] + ":" + str(e) + "\n")
        finally:
            pass
    err_out.close()


def main():
    # mongodb
    mongodb_addr = "12.11.70.12:10101"
    db_name = "ucsd"
    col_name = "icsd"

    myclient = pymongo.MongoClient(mongodb_addr)
    mydb = myclient[db_name]
    mycol = mydb[col_name]

    parser = argparse.ArgumentParser(description='MongoDB - icsd')
    parser.add_argument('output', help="error file output path")
    args = parser.parse_args()

    output_path = args.output
    data2db(mycol, output_path)


if __name__ == "__main__":
    main()
