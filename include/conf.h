#ifndef CONF_H
#define CONF_H

#include <string>
#include <fstream>
#include <iostream>
#include <regex>
#include <algorithm> 
#include <vector>
#include <map>
#include <set>

#include "exception.h"
#include "func.h"

using std::string;
using std::ifstream;
using std::ios;
using std::vector;
using std::map;
using std::pair;
using std::set;
using std::to_string;
using std::regex;
using std::regex_search;

/* model */

struct Solvent {
    string chem_formula;
    vector<string> elements;
    
    Solvent(string _chem_formula, vector<string> _elements) : \
        chem_formula(_chem_formula), elements(_elements) {}

    Solvent(const Solvent& other) {
        chem_formula = other.chem_formula;
        elements = other.elements;
    }

    friend bool operator<(const Solvent& a, const Solvent& b) {
        return a.chem_formula < b.chem_formula;
    }
};

struct Radius {
    string species;
    string radius;
    string O_metal;

    Radius(string _species, string _radius, string _O_metal) : \
        species(_species), radius(_radius), O_metal(_O_metal) {}

    Radius(const Radius &radius) {
        this->species = radius.species;
        this->radius = radius.radius;
        this->O_metal = radius.O_metal;
    }
};

struct Connected {
    string atom_a;
    string atom_b;

    Connected(string _a, string _b) : atom_a(_a), atom_b(_b){}

    bool operator==(const Connected &other) {
        return  this->atom_a == other.atom_a && this->atom_b == other.atom_b;
    }

    friend bool operator<(const Connected& a, const Connected& b) {
        return a.atom_a < b.atom_a || a.atom_b < b.atom_b;
    }
};

struct Cell {
    int atom_num;

    // n - array 1*3
    vector<vector<double>> trans_arr;

    // n - array 3*3
    vector<vector<vector<double>>> symm_arr;

    // 3 * 3
    double lattice[3][3];

    // atom coordinate, n * 3
    vector<vector<double>> position;

    // atom type
    vector<int> types;
};


/* data */

// chem_formula - Solvent
map<string ,Solvent> solvent_dict;

// species - Radius
map<string, Radius> radius_dict;

// elements, like [H, He, ...]
vector<string> elements;

// HM2Hall
map<string, string> HM2Hall;

// Hall2Number
map<string, int> Hall2Number;

// Number2Hall
map<int, string> Number2Hall;

// Hall2HM
map<string, string> Hall2HM;

// Rhomb2HexHall
map<string, string> Rhomb2HexHall;

// AP2Number
map<int, int> AP2Number;

// Number2AP
map<int, int> Number2AP;

// SymOpsHall
map<string, vector<vector<string>>> SymOpsHall;

/* IO Reader */

// convert element like fe -> Fe, c -> C
string element_format(string &ele) noexcept(false) {
    if(ele.size() == 1) {
        ele[0] = toupper(ele[0]);
    }
    else if(ele.size() == 2) {
        ele[0] = toupper(ele[0]);
        ele[1] = tolower(ele[1]);
    }
    else {
        throw Exception(ele + " is illegal");
    }
    return ele;
}

// get atom radius
void get_atom_radius() noexcept(false) {
    if(!radius_dict.empty()) {
        return;
    }

    ifstream in("./conf/radius.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("radius file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        in.close();
        
        vector<string> all_radius = re_split(str, "\\n+");
        for(auto rad : all_radius) {
            vector<string> rad_info = re_split(rad, "\\s+");
            Radius radius = Radius(rad_info[0], rad_info[1], rad_info.back());
            radius_dict.insert(pair<string, Radius>(rad_info[0], radius));
        }
    }
}

// get known solvent
void get_known_solvent() noexcept(false) {
    if(!solvent_dict.empty()) {
        return;
    }

    ifstream in("./conf/solvent.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("known solvent file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> all_solvent = re_split(str, "\\n+");
        for(auto sol : all_solvent) {
            vector<string> ele = del_split(sol, ' ');
            string formula = ele[0];
            ele.erase(ele.begin());

            Solvent s = Solvent(formula, ele);
            solvent_dict.insert(pair<string, Solvent>(formula, s));
        }
    }
}

// get elements
void get_elements() noexcept(false) {
    if(!elements.empty()) {
        return;
    }

    ifstream in("./conf/element.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("element file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        elements = re_split(str, "\\n+");
        for(auto &ele : elements) {
            ele = element_format(ele);
        }
    }
}

// get HM2Hall
void get_HM2Hall() noexcept(false) {
    if(!HM2Hall.empty()) {
        return;
    }

    ifstream in("./conf/HM2Hall.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("HM2Hall file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");
            HM2Hall.insert(pair<string, string>(info[0], info[1]));
        }

    }
}

// get Hall2Number
void get_Hall2Number() noexcept(false) {
    if(!Hall2Number.empty()) {
        return;
    }

    ifstream in("./conf/Hall2Number.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("Hall2Number file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");
            Hall2Number.insert(pair<string, int>(info[0], stoi(info[1])));
        }

    }
}

// get Number2Hall
void get_Number2Hall() noexcept(false) {
    if(!Number2Hall.empty()) {
        return;
    }

    ifstream in("./conf/Number2Hall.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("Number2Hall file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");
            Number2Hall.insert(pair<int, string>(stoi(info[0]), info[1]));
        }

    }
}

// get Hall2HM
void get_Hall2HM() noexcept(false) {
    if(!Hall2HM.empty()) {
        return;
    }

    ifstream in("./conf/Hall2HM.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("Hall2HM file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");
            Hall2HM.insert(pair<string, string>(info[0], info[1]));
        }
    }
}

// get Rhomb2HexHall
void get_Rhomb2HexHall() {
    if(!Rhomb2HexHall.empty()) {
        return;
    }

    ifstream in("./conf/Rhomb2HexHall.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("Rhomb2HexHall file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");
            Rhomb2HexHall.insert(pair<string, string>(info[0], info[1]));
        }
    }
}

// get AP2Number
void get_AP2Number() {
    if(!AP2Number.empty()) {
        return;
    }

    ifstream in("./conf/AP2Number.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("AP2Number file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");
            AP2Number.insert(pair<int, int>(stoi(info[0]), stoi(info[1])));
        }
    }
}

// get Number2AP
void get_Number2AP() {
    if(!Number2AP.empty()) {
        return;
    }

    ifstream in("./conf/Number2AP.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("Number2AP file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");
            Number2AP.insert(pair<int, int>(stoi(info[0]), stoi(info[1])));
        }
    }
}

// get SymOpsHall
void get_SymOpsHall() {
    if(!SymOpsHall.empty()) {
        return;
    }

    ifstream in("./conf/SymOpsHall.txt", ios::in);
    if(!in.is_open()) {
        throw Exception("SymOpsHall file does not exist or fails to open!");
    }
    else {
        string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        vector<string> lines = re_split(str, "\\n+");
        for(auto line : lines) {
            if(regex_search(line, regex("^\\s*#"))) {
                continue;
            }

            vector<string> info = re_split(line, "\\s+:\\s+");

            string key = info[0];

            vector<string> symm = del_split(info[1], ' ');
            vector<vector<string>> values;
            for(auto item : symm) {
                item = item.substr(1, item.size()-2);
                vector<string> tmp = del_split(item, ',');
                values.push_back(tmp);
            }
            SymOpsHall[key] = values;
        }
    }   
}

#endif