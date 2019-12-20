#ifndef CIF_FUNC_H
#define CIF_FUNC_H


#include "conf.h"
#include "exception.h"
#include "conf.h"
#include "cif.h"

#define EPS 1e-6


void get_res(bool log = true) noexcept(false) {
    if(log) {
        cout << "Getting some known resources..." << endl;
    }
    
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
            int index = i;

            if(item[i].find("x") != string::npos) {
                index = 0;
            }
            else if(item[i].find("y") != string::npos) {
                index = 1;
            }
            else {
                index = 2;
            }

            if(item[i].find(ex[index]) != string::npos) {
                s[index][index] = -1.0;
            }
            else {
                s[index][index] = 1.0;
            }

            t[index] = frac2double(item[i]);
        }

        symm_arr.push_back(s);
        trans_arr.push_back(t);
    }

    for(int i = 0; i < trans_arr.size(); i++) {
        for(int j = 0; j < symm_arr[i].size(); j++) {
            printVec(symm_arr[i][j]);
            cout << " ";
        }
        cout << endl;

        printVec(trans_arr[i]);
        cout << endl;
    }
}

// used in in-cell for modify x/y/z
double modify_num(double num) {
    if(num - 1.0 > EPS) {
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


// judge whether contains specific atoms
// input likes A&B or A
bool is_exist_atoms(CIF &cif, string &s) noexcept(false) {
    vector<string> atoms = del_split(s, '&');
    set<string> all_atoms = cif.get_atoms();

    for(auto &atom : atoms) {
        atom = element_format(atom);
        if(!all_atoms.count(atom)) {
            return false;
        }
    }

    return true;
}


// judge whether contains specific bonds
// input likes A-B&C-D or A-B
bool is_exist_bonds(CIF &cif, string &s) noexcept(false) {
    vector<string> bonds = del_split(s, '&');
    vector<Connected> connect;
    for(auto &bond : bonds) {
        vector<string> atoms = del_split(bond, '-');
        connect.push_back(Connected(element_format(atoms[0]), element_format(atoms[1])));
    }

    set<Connected> all_bonds = cif.get_bonds();
    for(auto &cn : connect) {
        if(!all_bonds.count(cn)) {
            return false;
        }
    }

    return true;
}



#endif