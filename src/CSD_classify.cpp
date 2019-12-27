#include <iostream>
#include <cctype>
#include "../include/cif.h"
#include "../include/cif_func.h"
#include "../include/cmdline.h"

using namespace std;

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_dir", 'i', "csd folder location", true, "");
    parser.add<string>("output_dir", 'o', "classification result export location", true);
    parser.add<double>("skin_distance", 'd', "the skin distance(coefficient) you want to use", false, 0.25);
    parser.add<string>("remove", 'r', "only remove the cif which contains special elements or special bonds(the input form likes special meatal/special bonds(Fe|Cu/Fe-O|C-O&C-H) or only input one of them, please use '/' as separators for elements and bonds)", false, "");
    parser.add<string>("keep", 'k', "only keep the cif which contains special elements and special bond(the input form likes special meatal/special bonds(Fe|Cu/Fe-O|C-O&C-H) or only input one of them, please use '/' as separators for elements and bonds", false, "");
    parser.add("log", 'l', "print the detail log, no log by default");
    parser.add("unique", 'u', "remove duplicate files");

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

    bool log = false;
    if(parser.exist("log")) {
        log = true;
    }

    bool unique = false;
    if(parser.exist("unique")) {
        unique = true;
    }

    bool keep = false, rm = false;
    string content;
    if(parser.exist("keep")) {
        keep = true;
        content = parser.get<string>("keep");
    }

    if(parser.exist("remove")) {
        rm = true;
        content = parser.get<string>("remove");
    }

    if(rm && keep) {
        cerr << "Please only choose remove or keep" << endl;
    }
    
    vector<string> sp_elements, sp_bonds;
    // deal with special content
    if(!content.empty()) {
        vector<string> strs = del_split(content, '/');
        if(strs.size() == 1) {
            if(strs[0].find("-") == string::npos) {
                // sp elements
                sp_elements = del_split(strs[0], '|');
            }
            else {
                // sp bonds
                sp_bonds = del_split(strs[0], '|');
            }
        }
        else if(strs.size() == 2) {
            sp_elements = del_split(strs[0], '|');
            sp_bonds = del_split(strs[1], '|');
        }
        else {
            cerr << "Please check your input form" << endl;
        }
    }


    string input = parser.get<string>("input_dir");
    string output = parser.get<string>("output_dir");

    if(!is_folder_exist(input)) {
        cout << "CSD data not found!" << endl;
    }

    vector<string> files = get_all_files(input, "cif");
    set<string> files_set;
    for(auto item : files) {
        try {
            if(log) {
                cout << "---------------- " << item << " ----------------" << endl;
            }

            if(unique && files_set.count(get_pname(get_filename(item)))) {
                if(log) {
                    cout << item << " file  duplicates!" << endl << endl;
                }
                continue;
            }

            CIF cif = CIF(item, log);
            cif.parse_file();
            cif.get_known_res();
            
            bool satisfy = true;

            // judge whether contain special metal                
            if(!rm && !keep) {
                set<string> atoms = cif.get_atoms();
                for(auto iter = atoms.begin(); iter != atoms.end(); iter++) {
                    if(cif.get_atom_state(*iter) == "1") {
                        satisfy = false;
                        break;
                    }
                }
                if(!satisfy) {
                    if(log) {
                        cout << item << " contains metal!" << endl << endl;
                    }
                    continue;
                }
            } 
            else if(rm) {
                for(auto &ele : sp_elements) {
                    if(is_exist_atoms(cif, ele)) {
                        satisfy = false;
                        break;
                    }
                }
                if(!satisfy) {
                    if(log) {
                        cout << item << " contains special elements!" << endl << endl;
                    }
                    continue;
                }
            }
            else {
                satisfy = false;
                for(auto &ele : sp_elements) {
                    if(is_exist_atoms(cif, ele)) {
                        satisfy = true;
                        break;
                    }
                }
                // if(sp_elements.empty()) {
                //     satisfy = true;
                // }
                if(!satisfy) {
                    if(log) {
                        cout << item << " doesn't contain special elements!" << endl << endl;
                    }
                    continue;
                }
            }
            
            // judge whether contain disorder-molecule
            cif.build_base_cell(skin_distance, coefficient);
            satisfy = cif.judge_if_have_disorder();
            if(!satisfy) {
                if(log) {
                    cout << item << " contains disorder-molecule!" << endl << endl;
                }
                continue;
            }

            // judge whether contain bonds
            if(!rm && !keep) {

            } 
            else if(rm) {
                for(auto &ele : sp_bonds) {
                    if(is_exist_bonds(cif, ele)) {
                        satisfy = false;
                        break;
                    }
                }
                if(!satisfy) {
                    if(log) {
                        cout << item << " contains special bonds!" << endl << endl;
                    }
                    continue;
                }
            }
            else {
                satisfy = false;
                for(auto &ele : sp_bonds) {
                    if(is_exist_bonds(cif, ele)) {
                        satisfy = true;
                        break;
                    }
                }
                if(!satisfy) {
                    if(log) {
                        cout << item << " doesn't contain special bonds!" << endl << endl;
                    }
                    continue;
                }
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
                if(log) {
                    cp_file(item, output + "/csd_solvent/");
                    cout << item << " contains known solvent!" << endl << endl;
                }
                continue;
            }

            // copy
            if(solvent_chk_list.size() == 1) {
                cp_file(item, output + "/csd_normal/");
            }
            else {
                cp_file(item, output + "/csd_warning/");
            }

            if(unique) {
                // cout << get_pname(get_filename(item)) << endl;
                files_set.insert(get_pname(get_filename(item)));
            }

            if(log) {
                cout << endl;
            }
        }
        catch(Exception err) {
            if(log) {
                cerr << err.msg << endl;
            }
            return -1;
        }
    }

    return 0;
}