#include <iostream>
#include <cctype>
#include "../include/cif.h"
#include "../include/cmdline.h"

using namespace std;



int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_dir", 'i', "csd folder location", true, "");
    parser.add<string>("output_dir", 'o', "classification result export location", true);
    parser.add<string>("specific_metal", 'm', "only remove specified metal elements(the input form likes Fe-Cu-Zn), default is all metal elements", false, "");
    parser.add<double>("skin_distance", 'd', "the skin distance(coefficient) you want to use", false, 0.25);

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    double skin_distance = 0.25;
    double coefficient = 0;
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

    set<string> sp_metal;
    if(parser.exist("specific_metal")) {
        string str = parser.get<string>("specific_metal");
        vector<string> sp = del_split(str, '-');
        for(auto &item : sp) {
            if(item.size() == 1) {
                item[0] = toupper(item[0]);
            }
            else if(item.size() == 2) {
                item[0] = toupper(item[0]);
                item[1] = tolower(item[1]);
            }
            else {
                throw Exception("The input metal " + item + " does not exist");
            }

            sp_metal.insert(item);
        }
    }

    string input = parser.get<string>("input_dir");
    string output = parser.get<string>("output_dir");

    if(!is_folder_exist(input)) {
        cout << "CSD data not found!" << endl;
    }

    vector<string> files = get_all_files(input, "cif");
    try {
        for(auto item : files) {
            cout << "---------------- " << item << " ----------------" << endl;
            CIF cif = CIF(item);
            cif.parse_file();
            cif.get_known_res();
            
            bool satisfy = true;

            // judge whether contain metal
            set<string> atoms = cif.get_atoms();
            for(auto iter = atoms.begin(); iter != atoms.end(); iter++) {
                if(cif.get_atom_state(*iter) == "1") {
                    if(sp_metal.empty() || sp_metal.count(*iter)) {
                        satisfy = false;
                    }
                    break;
                }
            }
            if(!satisfy) {
                cout << item << " contains metal!" << endl << endl;
                continue;
            }


            // judge whether contain disorder-molecule
            cif.build_base_cell(skin_distance, coefficient);
            satisfy = cif.judge_if_have_disorder();
            if(!satisfy) {
                cout << item << " contains disorder-molecule!" << endl << endl;
                continue;
            }


            // judge wheather contain known solvent
            cif.find_solvent(false);
            vector<string> solvent_chk_list = cif.get_solvent_chk_list();
            for(auto &sol : solvent_chk_list) {
                if(solvent_dict.count(sol) > 0) {
                    satisfy = false;
                    break;
                }
            }
            if(!satisfy) {
                cout << item << " contains known solvent!" << endl << endl;
                continue;
            }

            // copy
            if(solvent_chk_list.size() == 1) {
                cp_file(item, output + "/csd_normal/");
            }
            else {
                cp_file(item, output + "/csd_warning/");
            }
            cout << endl;
        }

    }
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }
    return 0;
}