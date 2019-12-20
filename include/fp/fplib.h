#ifndef FPLIB_H
#define FPLIB_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "rcov.h"
#include "fingerprint.h"
#include "apc.h"
#include <stdio.h>


void get_fp_nonperiodic(int nid, int nat, int ntyp, int types[], double rxyz[][3], int znucl[], double fp[]);

void get_fp_periodic(int lmax, int nat, int ntyp, int types[], double lat[3][3],
                     double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp, double **lfp, bool print);

void get_fp_periodic_short(int lmax, int nat, int ntyp, int types[], double lat[3][3],
                           double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp);

void get_fp_periodic_long(int lmax, int nat, int ntyp, int types[], double lat[3][3],
                          double rxyz[][3], int znucl[], int natx, double cutoff, double **lfp);

double get_fpdistance_periodic(int nat, int ntyp, int types[], int fp_len,
                               double **fp1, double **fp2, int f[]);

void get_fp_nonperiodic(int nid, int nat, int ntyp, int types[], double rxyz[][3], int znucl[], double fp[])
{

    int i, j;
    int lwork, lda, info;
    //double om[nid][nid];
    double **om;
    double amp[nat];
    double a[nid * nid];
    double w[nid];
    double wkopt;
    double *work;
    double rcov[nat];

    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov(znucl[types[i] - 1]);

    for (i = 0; i < nat; i++)
        amp[i] = 1.0;

    nid = 4 * nat;

    if ((om = (double **)malloc(sizeof(double) * nid)) == NULL)
    {
        // fprintf(stderr, "Memory could not be allocated.");
        // exit(1);
        throw Exception("Memory could not be allocated.");
    }
    for (i = 0; i < nid; i++)
    {
        if ((om[i] = (double *)malloc(sizeof(double) * nid)) == NULL)
        {
            // fprintf(stderr, "Memory could not be allocated.");
            // exit(1);
            throw Exception("Memory could not be allocated.");
        }
    }
    creat_om(4, nat, rxyz, rcov, amp, om);

    for (i = 0; i < nid; i++)
        for (j = 0; j < nid; j++)
            a[i * nid + j] = om[j][i];

    lda = nid;
    lwork = -1;
    dsyev("V", "U", &nid, a, &lda, w, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double *)malloc(lwork * sizeof(double));
    dsyev("V", "U", &nid, a, &lda, w, work, &lwork, &info);
    if (info > 0)
    {
        // fprintf(stderr, "Error: dsyev 1");
        // exit(1);
        throw Exception("Error: dsyev 1");
    }
    if (w[0] < -1E-12)
    {
        printf("w[0] = %g\n", w[0]);
        // fprintf(stderr, "Error: Negative w");
        // exit(1);
        throw Exception("Error: Negative w");
    }

    for (i = 0; i < nid; i++)
        fp[i] = w[nid - 1 - i];

    free(work);
    for (i = 0; i < nid; i++)
        free(om[i]);
    free(om);
}

void get_fp_periodic(int lmax, int nat, int ntyp, int types[], double lat[3][3],
                     double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp, double **lfp, bool print = false)
{
    int i, ixyz, flag = 0;
    double rcov[nat];
    int lseg, l;

    if (lmax == 0)
    {
        lseg = 1;
        l = 1;
    }
    else if (lmax == 1)
    {
        lseg = 4;
        l = 2;
    }
    else
    {
        // fprintf(stderr, "Error: ORBITAL.");
        // exit(1);
        throw Exception("Error: ORBITAL.");
    }

    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov(znucl[types[i] - 1]);

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only;  == 0: short fp only; > 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov, cutoff, lfp, sfp, print);
}

void get_fp_periodic_short(int lmax, int nat, int ntyp, int types[], double lat[3][3],
                           double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp)
{
    int i, ixyz, flag = 1;
    double rcov[nat], **lfp = NULL;
    int lseg, l;

    if (lmax == 0)
    {
        lseg = 1;
        l = 1;
    }
    else if (lmax == 1)
    {
        lseg = 4;
        l = 2;
    }
    else
    {
        // fprintf(stderr, "Error: ORBITAL.");
        // exit(1);
        throw Exception("Error: ORBITAL.");
    }

    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov(znucl[types[i] - 1]);

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only; > 0 short fp only; == 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov, cutoff, lfp, sfp);
}

void get_fp_periodic_long(int lmax, int nat, int ntyp, int types[], double lat[3][3],
                          double rxyz[][3], int znucl[], int natx, double cutoff, double **lfp)
{
    int i, ixyz, flag = -1;
    double rcov[nat], **sfp = NULL;
    int lseg, l;

    if (lmax == 0)
    {
        lseg = 1;
        l = 1;
    }
    else if (lmax == 1)
    {
        lseg = 4;
        l = 2;
    }
    else
    {
        // fprintf(stderr, "Error: ORBITAL.");
        // exit(1);
        throw Exception("Error: ORBITAL.");
    }

    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov(znucl[types[i] - 1]);

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only;  == 0: short fp only; > 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov, cutoff, lfp, sfp);
}

double get_fpdistance_periodic(int nat, int ntyp, int types[], int fp_len,
                               double **fp1, double **fp2, int f[])
{
    double fpd, cc, tt, costmp[nat][nat], *a;
    int iat, jat, ityp, i, j, k, ii, jj, inat, *ft;

    fpd = 0.0;
    inat = 0;

    for (ityp = 1; ityp <= ntyp; ityp++)
    {
        i = 0;
        for (iat = 0; iat < nat; iat++)
        {
            if (types[iat] == ityp)
            {
                i++;
                j = 0;
                for (jat = 0; jat < nat; jat++)
                {
                    if (types[jat] == ityp)
                    {
                        j++;
                        tt = 0.0;
                        for (k = 0; k < fp_len; k++)
                            tt += (fp1[iat][k] - fp2[jat][k]) * (fp1[iat][k] - fp2[jat][k]);
                        costmp[i - 1][j - 1] = sqrt(tt);
                    }
                }
            }
        }
        ft = (int *)malloc(sizeof(int) * i);
        a = (double *)malloc(sizeof(double) * i * i);
        //printf("nt %d %d\n", i,j);

        for (ii = 0; ii < i; ii++)
            for (jj = 0; jj < i; jj++)
                a[ii * i + jj] = costmp[ii][jj];
        apc(i, a, &cc, ft);

        for (k = 0; k < i; k++)
        {
            f[inat] = ft[k];
            inat++;
        }

        free(ft);
        free(a);
        ft = NULL;
        a = NULL;
        fpd += cc;
    }
    fpd /= sqrt(nat);

    return fpd;
}

#endif