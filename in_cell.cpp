#include <iostream>
#include <iomanip>
#include <fstream>
#include "include/cmdline.h"
#include "include/cif.h"
#include "include/cif_func.h"
#include "include/func.h"

using namespace std;

void export_in_cell_result(string input, string output, CIF &cif, map<string, set<vector<double>>> &all_atoms);

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_path", 'i', "input MOF cif file", true, "");
    parser.add<string>("output_path", 'o', "output file path", true, "");

    parser.parse_check(argc, argv);


    if(argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }
    
    string input = parser.get<string>("input_path");
    string output = parser.get<string>("output_path");

    
    try {
        get_res();

        CIF cif = CIF(input);

        cif.parse_file();

        map<string, set<vector<double>>> all_atoms = in_cell(cif);
        export_in_cell_result(input, output, cif, all_atoms);
    }
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }
}


void export_in_cell_result(string input, string output, CIF &cif, map<string, set<vector<double>>> &all_atoms) {
    // export
    cout << "Exporting in-cell result..." << endl;
    if(!is_folder_exist(output)) {
        make_dir(output);
    }
    
    string name = output + "/" + del_split(get_filename(input), '.')[0] + "_in_cell.cif";
    ofstream out(name);

    map<string, vector<string>> loop_dict = cif.get_loop_dict();

    if(loop_dict.count("chaos_loop")) {
        for(auto str : loop_dict["chaos_loop"]) {
            if(str.find("END") != string::npos) {
                continue;
            }
            if(str.find("_symmetry_") != string::npos || str.find("space_group") != string::npos) {
                continue;
            }
            out << str << endl;
        }
        out << "_symmetry_space_group_name_H-M    \"P 1\"" << endl;
        out << "_symmetry_int_tables_number 1" << endl;
    }

    if(loop_dict.count("sym_loop")) {
        out << "loop_" << endl;
        out << "_symmetry_equiv_pos_as_xyz" << endl;
        out << "'x, y, z'" << endl;
        for(auto str : loop_dict["sym_loop"]) {
            if(str.find("x") != string::npos || str.find("y") != string::npos || str.find("z") != string::npos) {
                continue;
            }
            if(str.find("END") != string::npos) {
                continue;
            }
            out << str << endl;
        }
    }

    if(loop_dict.count("type_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["type_loop"]) {
            if(str.find("END") != string::npos) {
                continue;
            }
            out << str << endl;
        }
    }

    if(loop_dict.count("bond_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["bond_loop"]) {
            if(str.find("END") != string::npos) {
                continue;
            }
            out << str << endl;
        }
    }

    if(loop_dict.count("unknow_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["unknow_loop"]) {
            if(str.find("END") != string::npos) {
                continue;
            }
            out << str << endl;
        }
    }

    if(loop_dict.count("aniso_site_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["aniso_site_loop"]) {
            if(str.find("END") != string::npos) {
                continue;
            }
            out << str << endl;
        }
    }


    out << "loop_" << endl;
    out << "_atom_site_label" << endl << "_atom_site_type_symbol" << endl \
        << "_atom_site_fract_x" << endl << "_atom_site_fract_y" << endl << "_atom_site_fract_z" << endl \
        << "_atom_site_B_iso_or_equiv" << endl;
    // map<string, set>
    for(auto i = all_atoms.begin(); i != all_atoms.end(); i++) {
        string species = i->first;
        int cnt = 1;
        // set
        for(auto j = (i->second).begin(); j != (i->second).end(); j++) {
            // vector x-y-z
            out << setw(10) << left << species + to_string(cnt++) << setw(10) << left << species << setw(15) << left << (*j)[0] << setw(15) << left << (*j)[1] << setw(15) << left << (*j)[2] << setw(15) << left << "1.000" << endl;
        }
    }

    out << endl << endl << "#END"; 
    out.close();

    cout << "Export file " << name << " successfully！" << endl;
}