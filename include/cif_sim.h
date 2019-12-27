#ifndef CIF_SIMILARITY_H
#define CIF_SIMILARITY_H

#include "conf.h"
#include "exception.h"
#include "conf.h"
#include "cif.h"
#include "./fp/fplib.h"

const int lmax = 0;
const double cutoff = 5.0;

int natx;

double get_fp_similarity(CIF &a, CIF &b);

void get_fp_periodic(CIF &a);

int get_natx(int nat, double lattice[3][3], double cutoff);

double get_cutoff();

/* base on fringerprint */
double get_fp_similarity(CIF &a, CIF &b) {
    int lseg, l, nat = 0;
    int ntype = elements.size();;

    int znucl[elements.size()];
    int i, j;

    for(i = 0; i < elements.size(); i++) {
        znucl[i] = i + 2;
    }

    if (lmax == 0) {
        lseg = 1; l = 1;
    } else if (lmax == 1) {
        lseg = 4; l = 2;
    } else {
        cerr << "ORBITAL ERROR."; 
        return -1;
    }

    // cif a    
    double a_lattice[3][3];
    vector<vector<double>> a_tmp = a.get_lattice();
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            a_lattice[j][i] = a_tmp[i][j]; 
        }
    }

    // async
    map<string, vector<double>> a_atom_cd = a.get_atom_cd();
    nat = a_atom_cd.size();
    double a_position[a_atom_cd.size()][3];
    int types[a_atom_cd.size()];
    i = 0;
    for(auto iter = a_atom_cd.begin(); iter != a_atom_cd.end(); iter++, i++) {
        string specie = a.get_species(iter->first);
        auto atom = find(elements.begin(), elements.end(), specie);
        if(atom == elements.end()) {
            throw Exception(specie + " 's information not found");
        }
        types[i] = distance(elements.begin(), atom);

        for(j = 0; j < 3; j++) {
            a_position[i][j] = iter->second[j];
        }
    }

    // // in-cell
    // map<string, set<vector<double>>> a_atom_cd = in_cell(a);
    // for(auto it = a_atom_cd.begin(); it != a_atom_cd.end(); it++) {
    //     nat += it->second.size();
    // }

    // natx = get_natx(nat, a_lattice, cutoff);

    // double a_position[nat][3];
    // int types[nat];
    // i = 0;
    // for(auto it = a_atom_cd.begin(); it != a_atom_cd.end(); it++) {
    //     string specie = it->first;
        
    //     auto atom = find(elements.begin(), elements.end(), specie);
    //     if(atom == elements.end()) {
    //         throw Exception(specie + " 's information not found");
    //     }
    //     int type = distance(elements.begin(), atom);

    //     for(auto ait = it->second.begin(); ait != it->second.end(); ait++, i++) {
    //         for(j = 0; j < 3; j++) {
    //             a_position[i][j] = (*ait)[j];
    //         }
    //         types[i] = type;
    //     }
    // }


    // cif b
    double b_lattice[3][3];
    vector<vector<double>> b_tmp = b.get_lattice();
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            b_lattice[j][i] = b_tmp[i][j]; 
        }
    }

    // async
    map<string, vector<double>> b_atom_cd = b.get_atom_cd();
    double b_position[b_atom_cd.size()][3];
    i = 0;
    for(auto iter = b_atom_cd.begin(); iter != b_atom_cd.end(); iter++, i++) {
        for(j = 0; j < 3; j++) {
            b_position[i][j] = iter->second[j];
        }        
    }

    // // in-cell
    // map<string, set<vector<double>>> b_atom_cd = in_cell(b);

    // double b_position[nat][3];
    // int types[nat];
    // i = 0;
    // for(auto it = b_atom_cd.begin(); it != b_atom_cd.end(); it++) {
    //     string specie = it->first;
        
    //     auto atom = find(elements.begin(), elements.end(), specie);
    //     if(atom == elements.end()) {
    //         throw Exception(specie + " 's information not found");
    //     }
    //     int type = distance(elements.begin(), atom);

    //     for(auto ait = it->second.begin(); ait != it->second.end(); ait++, i++) {
    //         for(j = 0; j < 3; j++) {
    //             b_position[i][j] = (*ait)[j];
    //         }
    //         types[i] = type;
    //     }
    // }

    natx = get_natx(nat, a_lattice, cutoff);

    double **a_sfp, **a_lfp, **b_sfp, **b_lfp;
    
    a_sfp = (double **) malloc(sizeof(double)*nat);
    a_lfp = (double **) malloc(sizeof(double)*nat);

    b_sfp = (double **) malloc(sizeof(double)*nat);
    b_lfp = (double **) malloc(sizeof(double)*nat);

    for(int i = 0; i < nat; i++ ) {
        a_sfp[i] = (double *) malloc(sizeof(double) * l * (ntype + 1));
        a_lfp[i] = (double *) malloc(sizeof(double) * (natx * lseg));
        
        b_sfp[i] = (double *) malloc(sizeof(double) * l * (ntype + 1));
        b_lfp[i] = (double *) malloc(sizeof(double) * (natx * lseg));
    }

    get_fp_periodic(lmax, nat, ntype, types, a_lattice, a_position, znucl, natx, cutoff, a_sfp, a_lfp);
    get_fp_periodic(lmax, nat, ntype, types, b_lattice, b_position, znucl, natx, cutoff, b_sfp, b_lfp);

    auto s = a.get_atom_cd_set();
    int f[s.size()];
    double dis = get_fpdistance_periodic(nat, ntype, types, natx * lseg, a_lfp, b_lfp, f);

    for(i = 0; i < nat; i++) {
        free(a_sfp[i]);
        free(a_lfp[i]);
        free(b_sfp[i]);
        free(b_lfp[i]);
        a_sfp[i] = NULL;
        a_lfp[i] = NULL;
        b_sfp[i] = NULL;
        b_lfp[i] = NULL;
    }
    free(a_sfp);
    free(a_lfp);
    free(b_sfp);
    free(b_lfp);

    return dis;
}

void get_fp_periodic(CIF &a) {
    int lseg, l, nat = 0;
    int ntype = elements.size();;

    int znucl[elements.size()];
    int i, j;

    for(i = 0; i < elements.size(); i++) {
        znucl[i] = i + 2;
    }

    if (lmax == 0) {
        lseg = 1; l = 1;
    } else if (lmax == 1) {
        lseg = 4; l = 2;
    } 

    // cif a    
    double a_lattice[3][3];
    vector<vector<double>> a_tmp = a.get_lattice();
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            a_lattice[j][i] = a_tmp[i][j]; 
        }
    }

    // in-cell
    // map<string, set<vector<double>>> a_atom_cd = in_cell(a);
    // for(auto it = a_atom_cd.begin(); it != a_atom_cd.end(); it++) {
    //     nat += it->second.size();
    // }

    // double a_position[nat][3];
    // int types[nat];
    // i = 0;
    // for(auto it = a_atom_cd.begin(); it != a_atom_cd.end(); it++) {
    //     string specie = it->first;
        
    //     auto atom = find(elements.begin(), elements.end(), specie);
    //     if(atom == elements.end()) {
    //         throw Exception(specie + " 's information not found");
    //     }
    //     int type = distance(elements.begin(), atom);

    //     for(auto ait = it->second.begin(); ait != it->second.end(); ait++, i++) {
    //         for(j = 0; j < 3; j++) {
    //             a_position[i][j] = (*ait)[j];
    //         }
    //         types[i] = type;
    //     }
    // }

    // async
    map<string, vector<double>> a_atom_cd = a.get_atom_cd();
    nat = a_atom_cd.size();
    natx = get_natx(nat, a_lattice, cutoff);
    double a_position[a_atom_cd.size()][3];
    int types[a_atom_cd.size()];
    i = 0;
    for(auto iter = a_atom_cd.begin(); iter != a_atom_cd.end(); iter++, i++) {
        string specie = a.get_species(iter->first);

        for(j = 0; j < 3; j++) {
            a_position[i][j] = iter->second[j];
        }
        
        auto atom = find(elements.begin(), elements.end(), specie);
        if(atom == elements.end()) {
            throw Exception(specie + " 's information not found");
        }
        types[i] = distance(elements.begin(), atom);
    }

    double **a_sfp, **a_lfp;
    
    a_sfp = (double **) malloc(sizeof(double)*nat);
    a_lfp = (double **) malloc(sizeof(double)*nat);

    for(int i = 0; i < nat; i++ ) {
        a_sfp[i] = (double *) malloc(sizeof(double) * l * (ntype + 1));
        a_lfp[i] = (double *) malloc(sizeof(double) * (natx * lseg));        
    }

    get_fp_periodic(lmax, nat, ntype, types, a_lattice, a_position, znucl, natx, cutoff, a_sfp, a_lfp, true);

    for(i = 0; i < nat; i++) {
        free(a_sfp[i]);
        free(a_lfp[i]);
        a_sfp[i] = NULL;
        a_lfp[i] = NULL;
    }
    free(a_sfp);
    free(a_lfp);
}

int get_natx(int nat, double lattice[3][3], double cutoff) {
    int ixyz = get_ixyz(lattice, cutoff);
    return pow(ixyz * 2, 3) * nat;
}

double get_cutoff() {
    return cutoff;
}

#endif