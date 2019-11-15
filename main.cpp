#include <iostream>
#include "include/cif.h"
#include "include/cmdline.h"

using namespace std;

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("cif_in", 'i', "input MOF cif file", true, "");
    parser.add<string>("output_path", 'o', "output filepath", false);
    parser.add<double>("skin_distance", 'd', "the skin distance(coefficient) you want to use", false, 0.25);
    parser.add("solvent", 's', "output the solvent was found");
    parser.add("force", 'f', "remove solvent molecules anyway");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cerr << parser.usage();
        return 0;
    }

    double skin_distance = 0.25;
    double coefficient = 0;
    if(parser.exist("skin_distance")) {
        double distance = parser.get<double>("distance");
        if(distance < 1.0) {
            skin_distance = distance;
            coefficient = -1.0;
        }
        else {
            coefficient = distance;
            skin_distance = -1.0;
        }
    }

    // judge remove solvent molecules anyway
    bool isForce = false;
    if(parser.exist("force")) {
        isForce = true;
    }

    CIF cif = CIF(parser.get<string>("cif_in"));

    try {
        cif.parse_file();
        cif.get_known_res();
        cif.build_base_cell(skin_distance, coefficient);
        cif.find_solvent();
        cif.export_modify_result(parser.get<string>("output_path"), isForce);
    } 
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }

    return 0;
}