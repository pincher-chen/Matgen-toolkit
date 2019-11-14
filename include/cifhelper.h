#ifndef IO_HELPER_H
#define IO_HELPER_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include <cmath>
#include <algorithm> 

#include "exception.h"
#include "func.h"


using std::string;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::vector;
using std::map;
using std::pair;

#define PI acos(-1);

double skin_distance = 0.25;
double coefficient = 0;

struct Solvent {
    string chem_formula;
    vector<string> elements;
    
    Solvent(string _chem_formula, vector<string> _elements) : \
        chem_formula(_chem_formula), elements(_elements) {}
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

struct Distance {
    string atom_a;
    string atom_b;
    double distance;

    Distance(string _a, string _b, double _dis) : \
        atom_a(_a), atom_b(_b), distance(_dis) {}
};


class CIF
{

public:
    static vector<Solvent> solvent_list;
    static map<string, Radius> radius_dict;
    
    CIF(string filename) { 
        this->filename = filename;
    }

    // get known solvent
    void get_known_solvent() throw(Exception) {
        if(!solvent_list.empty()) {
            return;
        }

        ifstream in("./conf/solvent.txt", ios::in);
        if(in == nullptr) {
            throw Exception("known solvent file does not exist or fails to open!");
        }
        else {
            string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            vector<string> all_solvent = re_split(str, "\\n+");
            for(auto sol : all_solvent) {
                vector<string> ele = del_split(sol, ' ');
                string formula = ele[0];
                ele.erase(ele.begin());
                solvent_list.push_back(Solvent(formula, ele));
            }
        }
        
        // for(auto item : solvent_list) {
        //     cerr << item.chem_formula << " ";
        //     for(auto it : item.elements) {
        //         cerr << it << " ";
        //     }
        //     cerr << endl;
        // }
    }

    // get atom radius
    void get_atom_radius() throw(Exception) {
        if(!radius_dict.empty()) {
            return;
        }

        ifstream in("./conf/radius.txt", ios::in);
        if(in == nullptr) {
            throw Exception("radius file does not exist or fails to open!");
        }
        else {
            string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

            
            vector<string> all_radius = re_split(str, "\\n+");
            for(auto rad : all_radius) {
                vector<string> rad_info = re_split(rad, "\\s+");
                Radius radius = Radius(rad_info[0], rad_info[1], rad_info.back());
                radius_dict.insert(pair<string, Radius>(rad_info[0], radius));
            }
        }
    }

    // split the information in the cif file and save in map
    void split_cif() throw(Exception) {
        ifstream in(this->filename, ios::in);
        if(in == nullptr) {
            throw Exception(this->filename + " does not exist or fails to open!");
        }
        else {
            string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            this->cif_buf = str;
            // cerr << str << endl;

            vector<string> cif_loop_blocks  = re_split(str, "\\s+?loop_\\n");

            for (string content : cif_loop_blocks) {
                if (content.find("###########") != string::npos ||
                    std::regex_search(content, std::regex("^data_.*?")) ||
                    std::regex_search(content, std::regex("^\\s*_symmetry_space_group\\w+?"))) {

                    vector<string> chaos_loop = del_split(content, '\n');
                    loop_dict["chaos_loop"] = chaos_loop;

                    // for (auto item : loop_dict["chaos_loop"]) {
                    //     cerr << item << endl;
                    // }
                }
                else if (std::regex_search(content, std::regex("^\\s*_symmetry_equiv_pos\\w+?"))) {
                    vector<string> sym_loop = del_split(content, '\n');
                    loop_dict["sym_loop"] = sym_loop;
                }
                else if (std::regex_search(content, std::regex("^\\s*_atom_site\\w+?"))) {
                    vector<string> site_loop = del_split(content, '\n');
                    loop_dict["site_loop"] = site_loop;
                }
                else if (std::regex_search(content, std::regex("^\\s*_atom_type\\w+?"))) {
                    vector<string> type_loop = del_split(content, '\n');
                    loop_dict["type_loop"] = type_loop;
                }
                else if (std::regex_search(content, std::regex("^\\s*_geom_bond\\w+?"))) {
                    vector<string> bond_loop = del_split(content, '\n');
                    loop_dict["bond_loop"] = bond_loop;
                }
                else {
                    vector<string> unknow_loop = del_split(content, '\n');
                    loop_dict["unknow_loop"] = unknow_loop;
                }
            }
        }
    }

    // get crystal cell parameters from cif file
    void seek_crystal_info() throw(Exception) {
        string str = this->cif_buf;
        std::regex rgx("(_cell_[l|a].*?)\\s+(.*)");
        
        try {
            for (std::sregex_token_iterator it(str.begin(), str.end(), rgx), end; it != end; it++) {
                string cur = it->str();
                // cerr << cur << endl;
                vector<string> params = del_split(cur, ' ');
                string key = params[0];
                double value = get_num(params[1]);

                if(fabs(value + 1.0) < 0.00001) {
                    throw Exception("");
                }

                // cerr << key << " " << value << endl;
                cell_params[key] = value;
            }
        }
        catch(...) {
            throw Exception("cell parameter extract failed");
        }
    }

    // calculate lattice contstant and angle
    void cal_crystal_info() throw(Exception) {
        try {
            double a = this->cell_params["_cell_length_a"];
            double b = this->cell_params["_cell_length_b"];
            double c = this->cell_params["_cell_length_c"];
            
            double ap = this->cell_params["_cell_angle_alpha"] / 180 * PI;
            double bt = this->cell_params["_cell_angle_beta"] / 180 * PI;
            double ga = this->cell_params["_cell_angle_gamma"] / 180 * PI;

            // cerr << a << " " << b << " " << c << " " << ap << " " << bt << " " << ga << endl;

            double bc2 = pow(b, 2) + pow(c, 2) - 2 * b * c * cos(ap);
            
            double h1 = a;
            double h2 = b * cos(ga);
            double h3 = b * sin(ga);
            double h4 = c * cos(bt);
            double h5 = (pow(h2 - h4, 2) + pow(h3, 2) + pow(c, 2) - pow(h4, 2) - bc2) / (2 * h3);
            double h6 = sqrt(pow(c, 2) - pow(h4, 2) - pow(h5, 2));
            
            this->cell.push_back({h1, 0.0, 0.0});
            this->cell.push_back({h2, h3, 0.0});
            this->cell.push_back({h4, h5, h6});

            this->angel.push_back(ap);
            this->angel.push_back(bt);
            this->angel.push_back(ga);
        }
        catch(...) {
            throw Exception("cell parameter extract failed");
        }
    }
    
    // get the site information of the atom in the loop
    void deal_site_loop() throw(Exception) {
        try {
            vector<string> site_loop = this->loop_dict["site_loop"];
            vector<string> data_value;

            for(auto data : site_loop) {
                if(std::regex_match(data, std::regex("\\s*_\\w+.*?"))) {
                    this->site_label.push_back(trim(data));
                }
                else {
                    data_value.push_back(data);
                }
            }

            for(auto &item : data_value) {
                if(item == "#END") {
                    break;
                }
                this->site_value.push_back(del_split(item, ' '));
            }

            // for(auto item : this->site_label) {
            //     cerr << item << " ";
            // }
            // cerr << endl;

            // for(auto item : this->site_value) {
            //     for(auto v : item) {
            //         cerr << v << " ";
            //     }
            //     cerr << endl;
            // }
        }
        catch(...) {
            throw Exception("failed to divide data");
        }
    }

    // get the coordinates information of the atom
    void get_atom_coordinates() throw(Exception) {
        try {
            int x_index = std::distance(std::begin(this->site_label), \
                                        std::find(this->site_label.begin(), this->site_label.end(), "_atom_site_fract_x"));
            int y_index = std::distance(std::begin(this->site_label), \
                                        std::find(this->site_label.begin(), this->site_label.end(), "_atom_site_fract_y"));
            int z_index = std::distance(std::begin(this->site_label), \
                                        std::find(this->site_label.begin(), this->site_label.end(), "_atom_site_fract_z"));

            int species_index = std::distance(std::begin(this->site_label), \
                                            std::find(this->site_label.begin(), this->site_label.end(), "_atom_site_type_symbol"));
            int atom_site_index = std::distance(std::begin(this->site_label), \
                                                std::find(this->site_label.begin(), this->site_label.end(), "_atom_site_label"));
            
            // cerr << x_index << " " << y_index << " " << z_index << " " << species_index << " " << atom_site_index << endl;
            
            for(auto atom_info : this->site_value) {
                double x = get_num(atom_info[x_index]);
                double y = get_num(atom_info[y_index]);
                double z = get_num(atom_info[z_index]);

                string species = atom_info[species_index];
                string atom_site = atom_info[atom_site_index];

                string key = atom_site + "-" + species;
                this->atom_cd[key] = vector<double>{x, y, z};
            }

            // printMap(this->atom_cd);

        }
        catch(std::exception &e) {
            throw Exception(string(e.what()));
        }
    }

    // convert fractional coordinates to real coordinates and calculate bond lengths between different atoms
    vector<double> get_cart(vector<double> &atom_cd) {
        vector<double> lx_vec = this->cell[0];
        vector<double> ly_vec = this->cell[1];
        vector<double> lz_vec = this->cell[2];

        double frac_a = atom_cd[0];
        double frac_b = atom_cd[1];
        double frac_c = atom_cd[2];

        double x_cart = lx_vec[0] * frac_a + ly_vec[0] * frac_b + lz_vec[0] * frac_c;
        double y_cart = lx_vec[1] * frac_a + ly_vec[1] * frac_b + lz_vec[1] * frac_c;
        double z_cart = lx_vec[2] * frac_a + ly_vec[2] * frac_b + lz_vec[2] * frac_c;

        return vector<double>{x_cart, y_cart, z_cart};
    }

    double clac_distance(vector<double> &a_cd, vector<double> &b_cd) {
        // coordinate transformation
        vector<double> a_cart = get_cart(a_cd);
        vector<double> b_cart = get_cart(b_cd);

        double distance = sqrt(pow(fabs(a_cart[0] - b_cart[0]), 2) + \
                                pow(fabs(a_cart[1] - b_cart[1]), 2) + \
                                pow(fabs(a_cart[2] - b_cart[2]), 2));
        return distance;
    }

    // calculate bond distance
    void calc_bond_distance() {
        for(auto a = this->atom_cd.begin(); a != this->atom_cd.end(); a++) {
            for(auto b = this->atom_cd.begin(); b != this->atom_cd.end(); b++) {
                if(a == b) {
                    continue;
                }
                vector<double> a_cd = a->second;
                vector<double> b_cd = b->second;            

                double distance = clac_distance(a_cd, b_cd);

                string name = a->first + "-" + b->first;

                this->atom_dist[name] = distance;
            }
        }

        // printMap(this->atom_dist);
    }

    void judge_if_have_bond() throw(Exception) {
        for(auto iter = this->atom_dist.begin(); iter != this->atom_dist.end(); iter++) {
            vector<string> values = del_split(iter->first, '-');
            string atom_a = values[1];
            string atom_b = values[3];
            double distance = iter->second;

            auto a = radius_dict.find(atom_a), b = radius_dict.find(atom_b);

            if(a == radius_dict.end()) {
                throw Exception(atom_a + " 's radius not found");
            }

            if(b == radius_dict.end()) {
                throw Exception(atom_b + " 's radius not found");
            }

            double radius_a = get_num(a->second.radius);
            double radius_b = get_num(b->second.radius);

            if(skin_distance > 0) {
                double bond_ab = distance - radius_a - radius_b;
                if(bond_ab < skin_distance) {
                    this->have_bond_list.push_back(iter->first);
                }
            }
            else {
                double bond_ab = (radius_a + radius_b) * coefficient;
                if(distance < bond_ab) {
                    this->have_bond_list.push_back(iter->first);
                }
            }
        }
        // printVec(this->have_bond_list);
        cerr << "The number of bonded atom pairs is " << this->have_bond_list.size() << endl;
    }



private:
    string filename;
    string cif_buf;
   
    map<string, vector<string>> loop_dict;

    // label - atom coordinates
    map<string, vector<double>> atom_cd;

    // atom&atom - distance
    map<string, double> atom_dist;
    vector<string> have_bond_list;

    map<string, double> cell_params;
    vector<vector<double>> cell;
    vector<double> angel;

    vector<string> site_label;
    vector<vector<string>> site_value;
};

vector<Solvent> CIF::solvent_list;
map<string, Radius> CIF::radius_dict;

#endif