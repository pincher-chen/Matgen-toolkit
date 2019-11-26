#include <iostream>
#include <iomanip>
#include <fstream>
#include "include/cmdline.h"
#include "include/cif.h"
#include "include/func.h"

#define eps 1e-6

using namespace std;

double modify_num(double num);

void get_symm_info(CIF &cif, vector<vector<double>> &trans_arr, vector<vector<vector<double>>> &symm_arr);

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
        CIF cif = CIF(input);

        cif.parse_file();

        // n - array 1*3
        vector<vector<double>> trans_arr;
        // n - array 3*3
        vector<vector<vector<double>>> symm_arr;
        
        get_symm_info(cif, trans_arr, symm_arr);

        map<string, vector<double>> atom_cd = cif.get_atom_cd();
        
        // species - atom coordinates
        map<string, set<vector<double>>> all_atoms;
        
        for(auto iter = atom_cd.begin(); iter != atom_cd.end(); iter++) {
            string species = cif.get_species(iter->first);
                        
            int size = symm_arr.size();
            for(int i = 0; i < size; i++) {
                double x = modify_num(symm_arr[i][0][0] * iter->second[0] + trans_arr[i][0]);
                double y = modify_num(symm_arr[i][1][1] * iter->second[1] + trans_arr[i][1]);
                double z = modify_num(symm_arr[i][2][2] * iter->second[2] + trans_arr[i][2]);
                all_atoms[species].insert(vector<double> {x, y, z});
            }

        }

        export_in_cell_result(input, output, cif, all_atoms);
    }
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }
}

void get_symm_info(CIF &cif, vector<vector<double>> &trans_arr, vector<vector<vector<double>>> &symm_arr) {
    vector<vector<string>> symm = cif.get_symm();
    string name_HM = cif.get_name_HM();
    string name_Hall = cif.get_name_Hall();

    if(symm.empty()) {
        if(name_Hall.empty() && name_HM.empty()) {
            cout << "P1 symmetry is assumed!" << endl;
            symm.push_back(vector<string>{"x", "y", "z"});
        }
        else if(!name_Hall.empty()) {
            symm = SymOpsHall[name_Hall];
        }
        else if(!name_HM.empty()) {
            symm = SymOpsHall[HM2Hall[name_HM]];
        }
    }


    for(auto item : symm) {
        vector<vector<double>> s(3, vector<double>(3, 0.0));
        vector<double> t(3, 0.0);
        vector<string> ex{"-x", "-y", "-z"};

        // x - y - z
        for(int i = 0; i < 3; i++) {
            if(item[i].find(ex[i]) != string::npos) {
                s[i][i] = -1.0;
            }
            else {
                s[i][i] = 1.0;
            }

            t[i] = frac2double(item[i]);
        }

        symm_arr.push_back(s);
        trans_arr.push_back(t);
    }
}

void export_in_cell_result(string input, string output, CIF &cif, map<string, set<vector<double>>> &all_atoms) {
    // export
    cout << "Exporting in cell result..." << endl;
    if(!is_folder_exist(output)) {
        make_dir(output);
    }
    
    string name = output + "/" + del_split(get_filename(input), '.')[0] + "_in_cell.cif";
    ofstream out(name);

    map<string, vector<string>> loop_dict = cif.get_loop_dict();

    if(loop_dict.count("chaos_loop")) {
        for(auto str : loop_dict["chaos_loop"]) {
            out << str << endl;
        }
    }

    if(loop_dict.count("sym_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["sym_loop"]) {
            out << str << endl;
        }
    }

    if(loop_dict.count("type_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["type_loop"]) {
            out << str << endl;
        }
    }

    if(loop_dict.count("bond_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["bond_loop"]) {
            out << str << endl;
        }
    }

    if(loop_dict.count("unknow_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["unknow_loop"]) {
            out << str << endl;
        }
    }

    if(loop_dict.count("aniso_site_loop")) {
        out << "loop_" << endl;
        for(auto str : loop_dict["aniso_site_loop"]) {
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
            cout.setf(ios::left);
            out << setw(10) << left << species + to_string(cnt++) << setw(10) << left << species << setw(15) << left << (*j)[0] << setw(15) << left << (*j)[1] << setw(15) << left << (*j)[2] << setw(15) << left << "1.000" << endl;
        }
    }

    out << endl << endl << "#END"; 
    out.close();

    cout << "Export file " << name << " successfully！" << endl;
}

double modify_num(double num) {
    if(num - 1.0 > eps) {
        return modify_num(num - 1.0);
    }
    else if(num < 0) {
        return modify_num(num + 1.0);
    }

    return num;
}