#include <iostream>
#include "include/cifhelper.h"
#include "include/cmdline.h"

using namespace std;

int main(int argc, char *argv[]) {
    // 命令行参数获取
    cmdline::parser parser;
    parser.add<string>("cif_in", 'i', "input MOF cif file", true, "");
    parser.add<string>("output_path", 'o', "output filepath", false);
    parser.add<double>("skin_distance", 'd', "the skin distance(coefficient) you want to use", false, 0.25);
    parser.add("solvent", 's', "output the solvent was found");
    parser.add("model", 'm', "remove solvent molecules anyway");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cerr << parser.usage();
        return 0;
    }

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

    CIF cif = CIF(parser.get<string>("cif_in"));

    try {
        cif.split_cif();
        cif.seek_crystal_info();
        cif.cal_crystal_info();
        cif.deal_site_loop();
        cif.get_known_solvent();
        cif.get_atom_radius();
        cif.get_atom_coordinates();
        cif.calc_bond_distance();
        cif.judge_if_have_bond();
        cif.connect_network();
        cif.find_solvent();
        cif.export_modify_result(parser.get<string>("output_path"));
    } 
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }


    return 0;
}