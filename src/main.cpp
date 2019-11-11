#include <iostream>
#include "iohelper.h"
#include "cmdline.h"

using namespace std;

int main(int argc, char *argv[]) {
    // 命令行参数获取
    cmdline::parser parser;
    parser.add<string>("cif_in", 'i', "input MOF cif file", true, "");
    parser.add<string>("output_path", 'o', "output filepath", true, "");
    parser.add<double>("skin_distance", 'd', "the skin distance(coefficient) you want to use", false, 0.25);
    parser.add("solvent", 's', "output the solvent was found");
    parser.add("model", 'm', "remove solvent molecules anyway");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cerr << parser.usage();
        return 0;
    }


    CIF cif = CIF(parser.get<string>("cif_in"));

    try {
        cif.split_cif();
    } 
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }


    return 0;
}