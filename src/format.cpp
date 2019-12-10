#include <iostream>
#include <iomanip>
#include <fstream>
#include "../include/cif.h"
#include "../include/cif_func.h"
#include "../include/cmdline.h"

using namespace std;

const string conf_dir = "../conf/";

/* export */
void export_gaussion_format_result(string input, string output, string mode, string coord_type, CIF &cif);

void export_vasp_format_result(string input, string output, string mode, string coord_type, CIF &cif);

int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input", 'i', "input file", true, "");

    parser.add<string>("output", 'o', "output path of the conversion result ", true);

    // mode - in-cell/asymmetric
    parser.add<string>("mode", 'm', "the mode of the format conversion(in-cell/asymmetric)", false, "asymmetric");

    // format - vasp_file/opt-freq
    parser.add<string>("type", 't', "the type of the result format(gjf/vasp), convert format to gaussion format or vasp format", true, "");

    // coord_type - cart/frac
    parser.add<string>("coord_type", 'c', "the type of the coordinate(fract/cart), the coordinates of the atom in the conversion result are fractional coordinates or cartesian coordinates", false, "fract");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }


    string mode = "asymmetric";
    if(parser.exist("mode")) {
        mode = parser.get<string>("mode");
    }

    string coord_type = "fract";
    if(parser.exist("coord_type")) {
        coord_type = parser.get<string>("coord_type");
    }

    string input = parser.get<string>("input");
    string output = parser.get<string>("output");
    string type = parser.get<string>("type");

    if(!is_folder_exist(input)) {
        cout << "file not found!" << endl;
        return -1;
    }


    try {
        get_res();

        CIF cif = CIF(input);
        cif.parse_file();

        if(type == "gjf") {
            export_gaussion_format_result(input, output, mode, coord_type, cif);
        }
        else if(type == "vasp") {
            export_vasp_format_result(input, output, mode, coord_type, cif);
        }
        else {
            cerr << "Not support" << endl;
        }
    }
    catch(Exception err) {
        cout << err.msg << endl;
        return -1;
    }
    return 0;
}

void export_gaussion_format_result(string input, string output, string mode, string coord_type, CIF &cif) {
    // export
    cout << "Exporting format result..." << endl;
    if(!is_folder_exist(output)) {
        make_dir(output);
    }

    string basic = conf_dir + "opt-freq.conf";
    ifstream in(basic, ios::in);
    string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();

    string name = output + "/" + del_split(get_filename(input), '.')[0] + ".gjf";
    ofstream out(name);
    out << str << endl;

    map<string, set<vector<double>>> atom_cd = cif.get_atom_cd_set();
    if(mode == "in-cell") {
        atom_cd = in_cell(cif);
    }

    // map<string, set>
    for(auto i = atom_cd.begin(); i != atom_cd.end(); i++) {
        string species = i->first;
        int cnt = 1;
        // set
        for(auto j = (i->second).begin(); j != (i->second).end(); j++) {
            // vector x-y-z

            vector<double> cd = *j;
            if(coord_type == "cart") {
                cd = cif.frac2cart(cd);
            }
            // out.setf(ios::showpoint);
            out.flags(ios::fixed);
            out.precision(8);
            out << setw(8) << left << " " + species << setw(20) << right << cd[0] << setw(15) << right << cd[1] << setw(15) << right << cd[2] << endl;
        }
    }
    out << endl;
    out.close();

    cout << "Export file " << name << " successfully!" << endl;
}

void export_vasp_format_result(string input, string output, string mode, string coord_type, CIF &cif) {
    // export
    cout << "Exporting format result..." << endl;
    if(!is_folder_exist(output)) {
        make_dir(output);
    }

    string name = output + "/" + del_split(get_filename(input), '.')[0] + ".vasp";
    ofstream out(name);
    out << del_split(get_filename(input), '.')[0] << endl;
    out << "1.0" << endl;

    out.flags(ios::fixed);
    out.precision(10);
    
    // box size
    vector<vector<double>> lattice = cif.get_lattice();
    // for(int i = 0; i < 3; i++) {
    //     for(int j = 0; j < 3; j++) {
    //         cout << lattice[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    
    for(int i = 0; i < 3; i++) {
        out << setw(20) << right << lattice[0][i] << setw(20) << right << lattice[1][i] << setw(20) << right << lattice[2][i] << endl;
    }

    map<string, set<vector<double>>> atom_cd = cif.get_atom_cd_set();
    if(mode == "in-cell") {
        atom_cd = in_cell(cif);
    }

    // atom type
    for(auto i = atom_cd.begin(); i != atom_cd.end(); i++) {
        out << setw(10) << right << i->first;        
    }
    out << endl;

    // atom number
    for(auto i = atom_cd.begin(); i != atom_cd.end(); i++) {
        out << setw(10) << right << i->second.size();        
    }
    out << endl;

    if(coord_type == "cart") {
        out << "Cartesian" << endl;
    }
    else {
        out << "Direct" << endl;
    }
    
    // map<string, set>
    for(auto i = atom_cd.begin(); i != atom_cd.end(); i++) {
        string species = i->first;
        int cnt = 1;
        // set
        for(auto j = (i->second).begin(); j != (i->second).end(); j++) {
            // vector x-y-z

            vector<double> cd = *j;
            if(coord_type == "cart") {
                cd = cif.frac2cart(cd);
            }
            // out.setf(ios::showpoint);
            out.flags(ios::fixed);
            out.precision(8);
            // out << setw(8) << left << " " + species + to_string(cnt++) << setw(20) << right << cd[0] << setw(15) << right << cd[1] << setw(15) << right << cd[2] << endl;
            out << setw(20) << right << cd[0] << setw(15) << right << cd[1] << setw(15) << right << cd[2] << endl;        
        }
    }
    cout << "Export file " << name << " successfully!" << endl;
}