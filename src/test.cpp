#include <iostream>
#include <cctype>
#include "../include/cif.h"
#include "../include/cif_func.h"
#include "../include/cif_sim.h"
#include "../include/conf.h"
#include "../include/cmdline.h"
#include "../include/fp/fplib.h"


void get_cell(CIF &cif, Cell *cell);

vector<double> test_get_fp_nonperiodic(CIF &cif, Cell *cell);

void test_get_fp_periodic(CIF &a, Cell *ac, CIF &b, Cell *bc);


int main(int argc, char *argv[]) {
    // cmd
    cmdline::parser parser;
    parser.add<string>("input_a", 'a', "icsd location", true, "");
    parser.add<string>("input_b", 'b', "icsd location", true, "");
    parser.add<string>("output_dir", 'o', "classification result export location", true);
    parser.add("log", 'l', "print the detail log, no log by default");

    parser.parse_check(argc, argv);
    
    if (argc == 1 || parser.exist("help")) {
        cout << parser.usage();
        return 0;
    }

    string input_a = parser.get<string>("input_a");
    string input_b = parser.get<string>("input_b");    
    string output = parser.get<string>("output_dir");
    
    get_res();
    CIF a = CIF(input_a);
    a.parse_file();
    
    CIF b = CIF(input_b);
    b.parse_file();
    double d =  get_fp_similarity(a, b);
    cout << d << endl;
    // Cell ac, bc;
    // get_cell(a, &ac);
    // get_cell(b, &bc);

    // test_get_fp_periodic(a, &ac, b, &bc);
}

void get_cell(CIF &cif, Cell *cell) {
    vector<vector<double>> tmp = cif.get_lattice();
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            cell->lattice[j][i] = tmp[i][j];
        }
    }

    map<string, vector<double>> atom_cd = cif.get_atom_cd();
    for(auto iter = atom_cd.begin(); iter != atom_cd.end(); iter++) {
        string specie = cif.get_species(iter->first);
        cell->position.push_back(iter->second);
        auto atom = find(elements.begin(), elements.end(), specie);
        if(atom == elements.end()) {
            throw Exception(specie + " 's information not found");
        }
        cell->types.push_back(distance(elements.begin(), atom));
    }

    cell->atom_num = atom_cd.size();
}

vector<double> test_get_fp_nonperiodic(CIF &cif, Cell *cell) {
    double position[cell->atom_num][3];
    for(int i = 0; i < cell->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position[i][j] = cell->position[i][j];
        }
    }

    int types[cell->atom_num];
    for(int i = 0; i < cell->types.size(); i++) {
        types[i] = cell->types[i];
    }

    int lseg, nid;
    double *fp;

    lseg = 4;
    nid = cell->atom_num * lseg;
    fp = (double *) malloc(sizeof(double) * nid);

    int znucl[elements.size()];
    for(int i = 0; i < elements.size(); i++) {
        znucl[i] = i + 2;
    }

    get_fp_nonperiodic(nid, cell->atom_num, elements.size(), types, position, znucl, fp);
    
    vector<double> vec;
    for(int i = 0; i < nid; i++) {
        vec.push_back(fp[i]);
    }
    free(fp);

    printVec(vec);
    cerr << endl;


    return vec;
}


void test_get_fp_periodic(CIF &a, Cell *ac, CIF &b, Cell *bc) {
    double cutoff = 1.0;
    int lmax = 0, lseg, l;
    int natx = 100;

    double position1[ac->atom_num][3];
    for(int i = 0; i < ac->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position1[i][j] = ac->position[i][j];
        }
    }

    double position2[bc->atom_num][3];
    for(int i = 0; i < bc->position.size(); i++) {
        for(int j = 0; j < 3; j++) {
            position2[i][j] = bc->position[i][j];
        }
    }

    int types1[ac->atom_num];
    for(int i = 0; i < ac->types.size(); i++) {
        types1[i] = ac->types[i];
    }

    int types2[bc->atom_num];
    for(int i = 0; i < bc->types.size(); i++) {
        types2[i] = bc->types[i];
    }

    int znucl[elements.size()];
    for(int i = 0; i < elements.size(); i++) {
        znucl[i] = i + 2;
    }

    int nat = ac->atom_num;
    int ntyp = elements.size();

    double **sfp1, **lfp1, **sfp2, **lfp2;
    if (lmax == 0) {
        lseg = 1; l = 1;
    } else if (lmax == 1) {
        lseg = 4; l = 2;
    } else {
        cerr << "ORBITAL ERROR."; 
        return;
    }

    sfp1 = (double **) malloc(sizeof(double)*nat);
    lfp1 = (double **) malloc(sizeof(double)*nat);

    for(int i = 0; i < nat; i++ ) {
        sfp1[i] = (double *) malloc(sizeof(double) * l * (ntyp + 1));
        lfp1[i] = (double *) malloc(sizeof(double) * (natx * lseg));
    }

    sfp2 = (double **) malloc(sizeof(double)*nat);
    lfp2 = (double **) malloc(sizeof(double)*nat);

    for(int i = 0; i < nat; i++ ) {
        sfp2[i] = (double *) malloc(sizeof(double) * l * (ntyp + 1));
        lfp2[i] = (double *) malloc(sizeof(double) * (natx * lseg));
    }

    get_fp_periodic(lmax, nat, ntyp, types1, ac->lattice, position1, znucl, natx,  cutoff, sfp1, lfp1);
    get_fp_periodic(lmax, nat, ntyp, types2, bc->lattice, position2, znucl, natx,  cutoff, sfp2, lfp2);

    // for(int i = 0; i < nat; i++) {
    //     for(int j = 0; j < natx * lseg; j++) {
    //         cout << lfp1[i][j] << " ";
    //     }
    //     cout << endl;
    //     for(int j = 0; j < natx * lseg; j++) {
    //         cout << lfp2[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    auto s = a.get_atom_cd_set();
    int f[s.size()];
    double res_long = get_fpdistance_periodic(nat, ntyp, types1, natx * lseg, lfp1, lfp2, f);
    double res_short = get_fpdistance_periodic(nat, ntyp, types1, l * (ntyp + 1), sfp1, sfp2, f);
    cout << res_long << " " << res_short << endl;

    for(int i = 0; i < nat; i++) {
        free(sfp1[i]);
        free(lfp1[i]);
        free(sfp2[i]);
        free(lfp2[i]);
        sfp1[i] = NULL;
        lfp1[i] = NULL;
        sfp2[i] = NULL;
        lfp2[i] = NULL;
    }
    free(lfp1);
    free(sfp1);
    free(lfp2);
    free(sfp2);
}