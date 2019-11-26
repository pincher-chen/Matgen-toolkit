#include <iostream>
#include <iomanip>
#include <fstream>
#include "include/cmdline.h"
#include "include/cif.h"
#include "include/func.h"

using namespace std;

double modify_num(double num);

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

        vector<vector<string>> symm = cif.get_symm();

        map<string, vector<double>> atom_cd = cif.get_atom_cd();
        
        // n - array 1*3
        vector<vector<double>> trans_arr;
        // n - array 3*3
        vector<vector<vector<double>>> symm_arr;
        
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

        // species - atom coordinates
        map<string, set<vector<double>>> all_atoms;

        for(auto iter = atom_cd.begin(); iter != atom_cd.end(); iter++) {
            string species = cif.get_species(iter->first);
            
            all_atoms[species].insert(iter->second);
            
            int size = symm_arr.size();
            for(int i = 0; i < size; i++) {
                double x = modify_num(symm_arr[i][0][0] * iter->second[0] + trans_arr[i][0]);
                double y = modify_num(symm_arr[i][1][1] * iter->second[1] + trans_arr[i][1]);
                double z = modify_num(symm_arr[i][2][2] * iter->second[2] + trans_arr[i][2]);
                // if(x < 0 || x > 1 || y < 0 || y > 1 || z < 0 || z > 1) {
                //     cerr << x << " " << y << " " << z << endl;
                //     printVec(iter->second);
                //     cerr << endl;

                //     cerr << i << endl;
                //     for(auto num : trans_arr[i]) {
                //         cerr << num << " " ;
                //     }
                //     cerr << endl;
                // }
                all_atoms[species].insert(vector<double> {x, y, z});
            }

        }

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
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }
}

double modify_num(double num) {
    if(num >= 0.0 && num <= 1.0) {
        return num;
    }
    else if(num > 1.0) {
        return num - 1.0;
    }
    else if(num < 0) {
        return num + 1.0;
    }
}