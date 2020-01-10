#include <iostream>
#include "../include/cif.h"
#include "../include/cmdline.h"

using namespace std;

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("cif_in", 'i', "input MOF cif file", true, "");
    parser.add<string>("output_path", 'o', "output filepath", false);
    parser.add("force", 'f', "remove solvent molecules anyway");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    string output = "";
    if(parser.exist("output_path")) {
        output = parser.get<string>("output_path");
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
        cif.build_base_cell();
        cif.find_solvent();        
        if(!output.empty()) {
            cif.export_modify_result(output, isForce);
        }
    } 
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }

    return 0;
}