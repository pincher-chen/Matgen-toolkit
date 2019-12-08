#ifndef CIF_FUNC_H
#define CIF_FUNC_H

#include "conf.h"
#include "cif.h"

#define eps 1e-6


void get_res() noexcept(false) {
    cout << "Getting some known resources..." << endl;
    get_atom_radius();

    get_known_solvent();

    get_elements();

    get_HM2Hall();

    get_Hall2Number();

    get_Number2Hall();

    get_Hall2HM();

    get_Rhomb2HexHall();

    get_AP2Number();

    get_Number2AP();

    get_SymOpsHall();
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

// used in in-cell for modify x/y/z
double modify_num(double num) {
    if(num - 1.0 > eps) {
        return modify_num(num - 1.0);
    }
    else if(num < 0) {
        return modify_num(num + 1.0);
    }

    return num;
}

// return in-cell atoms
map<string, set<vector<double>>> in_cell(CIF &cif) {
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

    return all_atoms;
}

#endif