#ifndef IO_HELPER_H
#define IO_HELPER_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include "exception.h"

using std::string;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::vector;
using std::map;


class CIF {

public:
    CIF(string filename) {
        this->filename = filename;
    }

    // divide strings according to delimiters
    static vector<string> del_split(const string &str, char delimiter) {
        vector<string> result;
        string cur = str;
        while (!cur.empty()) {
            int index = cur.find_first_of(delimiter);
            if (index == -1) {
                result.push_back(cur);
                cur.clear();
            }
            else {
                result.push_back(cur.substr(0, index));
                cur = cur.substr(index + 1, cur.size() - index - 1);
            }
        }
        return result;
    }

    // divide strings according to regular matches
    static vector<string> re_split(const string &str, string pattern) {
        string cur = str;
        std::regex rgx(pattern);
        vector<string> result;
        for (std::sregex_token_iterator iter(cur.begin(), cur.end(), rgx, -1), end; iter != end; ++iter) {
            // cout << iter->str() << endl;
            result.push_back(iter->str());
        }
        return result;
    }

    // split the information in the cif file and save in map
    void split_cif() throw(Exception) {
        ifstream in(this->filename, ios::in);
        if(in == nullptr) {
            throw Exception(this->filename + " does not exist or fails to open!");
        }
        else {
            string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            // cerr << str << endl;

            vector<string> cif_loop_blocks  = re_split(str, "\\s+?loop_\\n");

            for(string content : cif_loop_blocks) {
                if(content.find("###########") != string::npos || \
                    std::regex_match(content, std::regex("data_.*?")) || \
                    std::regex_match(content, std::regex("\\s*_symmetry_space_group\\w+?"))) {
                    
                    vector<string> chaos_loop = del_split(content, '\n');
                    loop_dict["chaos_loop"] = chaos_loop;

                    // for (auto item : loop_dict["chaos_loop"]) {
                    //     cerr << item << endl;
                    // }
                }
                else if(std::regex_match(content, std::regex("\\s*_symmetry_equiv_pos\\w+?"))) {
                    vector<string> sym_loop = del_split(content, '\n');
                    loop_dict["sym_loop"] = sym_loop;
                }
                else if(std::regex_match(content, std::regex("\\s*_atom_site\\w+?"))) {
                    vector<string> site_loop = del_split(content, '\n');
                    loop_dict["site_loop"] = site_loop;
                }
                else if(std::regex_match(content, std::regex("\\s*_atom_type\\w+?"))) {
                    vector<string> type_loop = del_split(content, '\n');
                    loop_dict["type_loop"] = type_loop;
                }
                else if (std::regex_match(content, std::regex("\\s*_geom_bond\\w+?"))) {
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
private:
    string filename;
    map<string, vector<string>> loop_dict;
};

#endif