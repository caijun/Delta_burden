#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void readVector(const char *path, int nrow, int vector[])
{
    FILE *fp;
    fp = fopen(path, "r");
    for (int i = 0; i < nrow; i++)
    {
        fscanf(fp, "%d", &vector[i]);
    }
    fclose(fp);
}

void readVector(const char *path, int nrow, double vector[])
{
    FILE *fp;
    fp = fopen(path, "r");
    for (int i = 0; i < nrow; i++)
    {
        fscanf(fp, "%lf", &vector[i]);
    }
    fclose(fp);
}

vector<int> readVector(const char *path, vector<int> vector)
{
    char textline[100];
    char *token;
    ifstream infile(path);
    if (!infile)
    {
        cout << "read input file fails." << endl;
    }

    while (infile)
    {
        infile.getline(textline, 100);
        if (infile)
        {
            token = strtok(textline, " ");
            int num = atoi(token);
            vector.push_back(num);
        }
    }
    infile.close();

    return (vector);
}

vector<vector<double>> readMatrix(const char *path, int nrow, int ncol, vector<vector<double>> mat)
{
    FILE *fp;
    fp = fopen(path, "r");
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            fscanf(fp, "%lf", &mat[i][j]);
        }
    }
    fclose(fp);

    return (mat);
}

vector<vector<int>> readMatrix(const char *path, int nrow, int ncol, vector<vector<int>> mat)
{
    FILE *fp;
    fp = fopen(path, "r");
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            fscanf(fp, "%d", &mat[i][j]);
        }
    }
    fclose(fp);

    return (mat);
}

double max_eigenval(vector<vector<double>> mat, int n)
{
    double data[n * n];
    int idx = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            data[idx] = mat[i][j];
            idx++;
        }
    }

    gsl_matrix_view m = gsl_matrix_view_array(data, n, n);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(n);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n, n);
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n);
    gsl_eigen_nonsymmv(&m.matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free(w);
    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

    gsl_complex eval_0 = gsl_vector_complex_get(eval, 0);
    double max_eval = GSL_REAL(eval_0);

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);

    return (max_eval);
}

int sum(int a[], int n)
{
    int sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += a[i];
    }

    return (sum);
}

int main(int argc, char *argv[])
{
    //input parameters
    char *dirName = argv[1];      //output directory
    int Tonset = atoi(argv[2]);   //date index of the epidemic onset since November 30, 2020; default value is 366 (note that index is zero-based in C/C++), representing the epidemic will start at 2021-12-01
    int strategy = atoi(argv[3]); //vaccination strategy; candidate values are 1, 2
    int capacity = atoi(argv[4]); //the fold of maximum capacity of daily vaccine doses; candidate values are 1, 2, 3, and 4, representing 1.25, 1.5, 1.75, and 2 folds respectively
    double ve2 = atof(argv[5]);   //the vaccine efficacy after 14 days of the 2nd dose; default is 80%
    double R0 = atof(argv[6]);    //R0
    int nsim = atoi(argv[7]);     //number of simulations
    int susflag = atoi(argv[8]);  //susceptibility to infection by age (1 = heterogeneous, 2 = homogeneous)
    int cmflag = atoi(argv[9]);   //in which period contact matrix will be used (1 = baseline, 2 = postlockdown)
    int ni = atoi(argv[10]);      //number of initial seed infectors
    double gt = atof(argv[11]);   //generation time (days)

    //set parameters
    int Tmax = 365;
    int zeta = 4;
    int ngr = 16;
    int delay1 = 21, delay2 = 14;
    double gamma = 1.0 / gt;
    double prob_gamma = 1.0 - exp(-gamma / zeta); //probability of transition from I to R at each time step
    double ve1 = 0.67 / 0.8 * ve2;                //the vaccine efficacy right after administration of the 2nd dose. note that the vaccine efficacy after the 1st dose is 0%.
    double beta, max_eval;
    double prop_ci[ngr];
    double prop_preg[ngr];
    double prop_im[ngr];
    vector<int> doses;
    vector<vector<double>> sus(nsim, vector<double>(ngr));
    vector<vector<double>> cm(ngr, vector<double>(ngr));
    vector<vector<int>> vax_strategy(ngr, vector<int>(4));

    //epidemiological state variables
    int N[ngr], S[ngr], C[ngr], I[ngr], R[ngr];
    int V0mat[ngr][delay1 * zeta], V1mat[ngr][delay2 * zeta];
    int V0[ngr], V1[ngr], V2[ngr];

    int newI[ngr];                                //counter of daily number of new infections by age
    int newI_V0[ngr], newI_V1[ngr], newI_V2[ngr]; //counter of number of infected people by age with vaccine administration at each day
    int nS_I[ngr], nC_I[ngr], nI_R[ngr];
    int nV0_Imat[ngr][delay1 * zeta]; //number of people progressing from V0 to I by age group and time step
    int nV1_Imat[ngr][delay2 * zeta]; //number of people progressing from V1 to I by age group and time step
    int nV0_I[ngr], nV1_I[ngr], nV2_I[ngr];
    int nS_V0[ngr], nV0_V1[ngr], nV1_V2[ngr], S_vax[ngr];
    int dose1[ngr];    //counter of number of people by age administered the 1st dose at each day
    int dose2[ngr];    //counter of number of people by age administer the 2nd dose at each day
    int dose2eff[ngr]; //counter of number of people by age with the 2nd dose becoming effective at each day

    //read data
    readVector("data/population", ngr, N);             //population by age
    readVector("data/contraindication", ngr, prop_ci); //the proportion of contraindications by age
    readVector("data/pregnant", ngr, prop_preg);       //the proportion of pregnant women by age
    readVector("data/initial_immunity", ngr, prop_im); //the initial immunity by age before epidemic onset
    char path[100];
    sprintf(path, "data/doses%d", capacity);
    doses = readVector(path, doses); //daily number of doses administered
    //relative susceptibility by age
    if (susflag == 1)
    {
        sus = readMatrix("data/susceptibility/susceptibility_heter", nsim, ngr, sus);
    }
    else if (susflag == 2)
    {
        sus = readMatrix("data/susceptibility/susceptibility_homo", nsim, ngr, sus);
    }

    vax_strategy = readMatrix("data/strategy", ngr, 4, vax_strategy); //vaccination strategy
    //output file
    char outfile[100];
    sprintf(outfile, "%s/newI_doses_R0-%3.1lf_ve2-%.3lf_Tonset-%d_strategy-%d_nsim-%d_sus-%d_cm-%d_seed-%d_gt-%3.1lf.txt", 
    dirName, R0, ve2, Tonset, strategy, nsim, susflag, cmflag, ni, gt);
    FILE *fpout;
    fpout = fopen(outfile, "w");

    //initialize gsl random generator
    unsigned int seed = 12306;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);

    for (int sim = 0; sim < nsim; sim++) //contact matrix
    {
        if (cmflag == 1)
        {
            sprintf(path, "data/cm/baseline/contact_matrix%d", sim + 1);
        }
        else if (cmflag == 2)
        {
            sprintf(path, "data/cm/postlockdown/contact_matrix%d", sim + 1);
        }
        cm = readMatrix(path, ngr, ngr, cm);

        for (int i = 0; i < ngr; i++)
        {
            for (int j = 0; j < ngr; j++)
            {
                cm[i][j] = cm[i][j] * sus[sim][i];
            }
        }

        max_eval = max_eigenval(cm, ngr);
        beta = R0 * gamma / max_eval;

        //initialize state variables
        for (int i = 0; i < ngr; i++)
        {
            S[i] = round(N[i] * (1.0 - prop_im[i] - prop_ci[i] - prop_preg[i])); //susceptibles by age who can be vaccinated
            C[i] = round(N[i] * (prop_ci[i] + prop_preg[i]));                    //contraindications and pregnant women by age who cannot be vaccinated
            I[i] = 0;
            R[i] = 0;

            for (int j = 0; j < delay1 * zeta; j++)
            {
                V0mat[i][j] = 0;
            }

            for (int j = 0; j < delay2 * zeta; j++)
            {
                V1mat[i][j] = 0;
            }

            V0[i] = 0;
            V1[i] = 0;
            V2[i] = 0;

            newI[i] = 0;
            newI_V0[i] = 0;
            newI_V1[i] = 0;
            newI_V2[i] = 0;

            dose1[i] = 0;
            dose2[i] = 0;
            dose2eff[i] = 0;
        }

        //start to simulate transmission dynamics
        int Td = 0, Td_doses, gr_vax_idx;
        for (int t = 0; t < (Tonset + Tmax) * zeta; t++)
        {
            if (strategy == 1)
            {
                if (Td < 119)
                { //first vaccinate 18-59 between 2020-11-30 and 2021-03-29

                    gr_vax_idx = 0;
                }
                else if (Td < 244)
                { //vaccinate 18+ between 2021-03-29 and 2021-08-01
                    gr_vax_idx = 1;
                }
                else
                { //vaccinate 12+ since 2021-08-01
                    gr_vax_idx = 2;
                }
            }
            else if (strategy == 2)
            {
                if (Td < 119)
                { //first vaccinate 18-59 between 2020-11-30 and 2021-03-29
                    gr_vax_idx = 0;
                }
                else if (Td < 244)
                { //vaccinate 18+ between 2021-03-29 and 2021-08-01
                    gr_vax_idx = 1;
                }
                else if (Td < 332)
                { //vaccinate 12+ between 2021-08-01 and 2021-10-28
                    gr_vax_idx = 2;
                }
                else
                { //vaccinate 3+ since 2021-10-28
                    gr_vax_idx = 3;
                }
            }

            if (t == Tonset * zeta) //Tonset == 0, i.e., 2020-11-30
            {
                //seed initial infections
                for (int i = 0; i < ni; i++)
                {
                    //randomly select age group
                    int q = gsl_rng_uniform_int(rng, ngr * 2);
                    if (q < ngr)
                    {
                        if (S[q] == 0)
                        {
                            i--;
                        }
                        else
                        {
                            S[q]--;
                            I[q]++;
                        }
                    }
                    else
                    {
                        q -= ngr;
                        if (C[q] == 0)
                        {
                            i--;
                        }
                        else
                        {
                            C[q]--;
                            I[q]++;
                        }
                    }
                }
            }

            //compute age-dependent force of infection
            for (int i = 0; i < ngr; i++)
            {
                double foi = 0;
                for (int j = 0; j < ngr; j++)
                {
                    foi += beta * cm[i][j] * (double)I[j] / (double)N[j];
                }

                double prob_foi = 1.0 - exp(-foi / zeta);
                double prob_foi_ve1 = 1.0 - exp(-(1.0 - ve1) * foi / zeta);
                double prob_foi_ve2 = 1.0 - exp(-(1.0 - ve2) * foi / zeta);

                //epidemiological transitions
                //new infections
                nS_I[i] = gsl_ran_binomial(rng, prob_foi, S[i]);
                nC_I[i] = gsl_ran_binomial(rng, prob_foi, C[i]);
                //among vaccinated V0
                nV0_I[i] = 0;
                for (int j = 0; j < delay1 * zeta; j++)
                {
                    nV0_Imat[i][j] = gsl_ran_binomial(rng, prob_foi, V0mat[i][j]);
                    nV0_I[i] += nV0_Imat[i][j];
                }
                //among vaccinated V1
                nV1_I[i] = 0;
                for (int j = 0; j < delay2 * zeta; j++)
                {
                    nV1_Imat[i][j] = gsl_ran_binomial(rng, prob_foi_ve1, V1mat[i][j]);
                    nV1_I[i] += nV1_Imat[i][j];
                }
                //among vaccinated V2
                nV2_I[i] = gsl_ran_binomial(rng, prob_foi_ve2, V2[i]);
                //progression from infection to recovered
                nI_R[i] = gsl_ran_binomial(rng, prob_gamma, I[i]);

                //update state variables due to transmission
                S[i] -= nS_I[i];
                C[i] -= nC_I[i];
                V0[i] -= nV0_I[i];
                V1[i] -= nV1_I[i];
                V2[i] -= nV2_I[i];
                I[i] += nS_I[i] + nC_I[i] + nV0_I[i] + nV1_I[i] + nV2_I[i] - nI_R[i];
                R[i] += nI_R[i];
                //update V0mat and V1mat
                for (int j = 0; j < delay1 * zeta; j++)
                {
                    V0mat[i][j] -= nV0_Imat[i][j];
                }

                for (int j = 0; j < delay2 * zeta; j++)
                {
                    V1mat[i][j] -= nV1_Imat[i][j];
                }

                newI[i] += nS_I[i] + nC_I[i] + nV0_I[i] + nV1_I[i] + nV2_I[i];
                newI_V0[i] += nV0_I[i];
                newI_V1[i] += nV1_I[i];
                newI_V2[i] += nV2_I[i];
            }

            //vaccination at each time step
            //doses available for allocation at day Td
            if (Td >= doses.size())
            {
                Td_doses = doses.back();
            }
            else
            {
                Td_doses = doses[Td];
            }

            //first allocate doses to those people who need to vaccinate the 2nd dose
            //number of people progressing from V0 to V1 at each time step, namely the number of people administered the 2nd dose at each time step
            for (int i = 0; i < ngr; i++)
            {
                nV0_V1[i] = V0mat[i][delay1 * zeta - 1];
                dose2[i] += nV0_V1[i];
            }
            int nd2 = sum(nV0_V1, ngr);
            //then allocate the 1st doses
            int nd1 = round((double)Td_doses / (double)zeta) - nd2;
            if (nd1 > 0)
            {
                //number of susceptibles can be given the 1st dose
                for (int i = 0; i < ngr; i++)
                {
                    S_vax[i] = S[i] * vax_strategy[i][gr_vax_idx];
                }
                int total_S_vax = sum(S_vax, ngr);
                if (nd1 <= total_S_vax)
                {
                    for (int i = 0; i < ngr; i++)
                    {
                        nS_V0[i] = round((double)nd1 * (double)S_vax[i] / (double)total_S_vax);
                    }
                }
                else
                {
                    for (int i = 0; i < ngr; i++)
                    {
                        nS_V0[i] = S_vax[i];
                    }
                }
            }
            else
            {
                for (int i = 0; i < ngr; i++)
                {
                    nS_V0[i] = 0;
                }
            }

            for (int i = 0; i < ngr; i++)
            {
                dose1[i] += nS_V0[i];
            }

            //number of people progressing from V1 to V2 at each time step
            for (int i = 0; i < ngr; i++)
            {
                nV1_V2[i] = V1mat[i][delay2 * zeta - 1];
                dose2eff[i] += nV1_V2[i];
            }

            //update state variables due to vaccination
            for (int i = 0; i < ngr; i++)
            {
                S[i] -= nS_V0[i];
                V0[i] = V0[i] + nS_V0[i] - nV0_V1[i];
                V1[i] = V1[i] + nV0_V1[i] - nV1_V2[i];
                V2[i] = V2[i] + nV1_V2[i];

                //update V0mat
                for (int j = (delay1 * zeta - 1); j > 0; j--)
                {
                    V0mat[i][j] = V0mat[i][j - 1];
                }
                V0mat[i][0] = nS_V0[i];
                //update V1mat
                for (int j = (delay2 * zeta - 1); j > 0; j--)
                {
                    V1mat[i][j] = V1mat[i][j - 1];
                }
                V1mat[i][0] = nV0_V1[i];
            }

            //check state variables
            for (int i = 0; i < ngr; i++)
            {
                if (S[i] < 0 || C[i] < 0 || V0[i] < 0 || V1[i] < 0 || V2[i] < 0 || I[i] < 0 || R[i] < 0)
                {
                    cout << "Negative population in S/C/V0/V1/V2/I/R" << endl;
                    cout << "t = " << t << "\t" << S[i] << "\t" << C[i] << "\t" << V0[i] << "\t" << V1[i] << "\t" << V2[i] << "\t" << I[i] << "\t" << R[i] << endl;
                }
            }

            //output daily number of new infections, 1st dose, and 2nd dose by age at each day
            if (t % zeta == 0)
            {
                if (sim == 0 && t == 0) //header
                {
                    fprintf(fpout, "sim t");
                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " newI.%d", i + 1);
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " newI_V0.%d", i + 1);
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " newI_V1.%d", i + 1);
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " newI_V2.%d", i + 1);
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " dose1.%d", i + 1);
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " dose2.%d", i + 1);
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " dose2eff.%d", i + 1);
                    }
                    fprintf(fpout, "\n");
                }

                fprintf(fpout, "%d %d", sim, Td);
                for (int i = 0; i < ngr; i++)
                {
                    fprintf(fpout, " %d", newI[i]);
                    newI[i] = 0;
                }

                for (int i = 0; i < ngr; i++)
                {
                    fprintf(fpout, " %d", newI_V0[i]);
                    newI_V0[i] = 0;
                }

                for (int i = 0; i < ngr; i++)
                {
                    fprintf(fpout, " %d", newI_V1[i]);
                    newI_V1[i] = 0;
                }

                for (int i = 0; i < ngr; i++)
                {
                    fprintf(fpout, " %d", newI_V2[i]);
                    newI_V2[i] = 0;
                }

                for (int i = 0; i < ngr; i++)
                {
                    fprintf(fpout, " %d", dose1[i]);
                    dose1[i] = 0;
                }

                for (int i = 0; i < ngr; i++)
                {
                    fprintf(fpout, " %d", dose2[i]);
                    dose2[i] = 0;
                }

                for (int i = 0; i < ngr; i++)
                {
                    fprintf(fpout, " %d", dose2eff[i]);
                    dose2eff[i] = 0;
                }
                fprintf(fpout, "\n");

                Td++;
            }
        }
    }

    fclose(fpout);

    //check data and parameters
    cout << "outdir: " << dirName << endl;
    cout << "epidemic starts:" << Tonset << endl;
    cout << "vaccination strategy: " << strategy << endl;
    cout << "fold in daily doses: " << capacity << endl;
    cout << "vaccine efficacy: " << ve2 << endl;
    cout << "R0: " << R0 << endl;
    cout << "nsim: " << nsim << endl;
    cout << "susceptibility to infection: " << susflag << endl;
    cout << "contact matrix: " << cmflag << endl;
    cout << "number of seed infectiors: " << ni << endl;
    cout << "generation time: " << gt << endl;

    return 0;
}
