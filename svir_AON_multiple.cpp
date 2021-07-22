#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <numeric>
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

vector<double> readVector(const char *path, vector<double> vector)
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
            double num = atof(token);
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
    //input parameters from command line
    char *dirName = argv[1];      //output directory
    int Tonset = atoi(argv[2]);   //date index of the epidemic onset since November 30, 2020; default value is 275 (note that index is zero-based in C/C++), representing the epidemic will start at 2021-09-01
    int strategy = atoi(argv[3]); //vaccination strategy; candidate values are 1, 2
    int capacity = atoi(argv[4]); //the fold of maximum capacity of daily vaccine doses; candidate values are 1, 2, 3, and 4, representing 1.25, 1.5, 1.75, and 2 folds respectively
    int nsim = atoi(argv[5]);     //number of simulations
    //input parameters from reading files
    vector<double> ve2; //the vaccine efficacy after 14 days of the 2nd dose for wild-type and 4 variant strains; default is 80%
    ve2 = readVector("data/multiple/param/ve2", ve2);
    vector<double> R0; //R0 for wild-type and 4 variant strains (alpha, beta, gamma, delta)
    R0 = readVector("data/multiple/param/R0", R0);
    vector<double> init_f; //initial fractions of wild-type and 4 variant strains
    init_f = readVector("data/multiple/param/init_f", init_f);

    //set parameters
    int nvar = init_f.size();
    int Tmax = 365;
    int zeta = 4;
    int ni = 40;
    int ni_vec[nvar];
    double ve1[nvar]; //the vaccine efficacy right after administration of the 2nd dose. note that the vaccine efficacy after the 1st dose is 0%.
    for (int i = 0; i < nvar; i++)
    {
        ni_vec[i] = round(ni * init_f[i]);
        ve1[i] = 0.67 / 0.8 * ve2[i];
    }
    const char *variants[nvar] = {"Wild-type", "Alpha", "Beta", "Gamma", "Delta"};

    int ngr = 16;
    int delay1 = 21, delay2 = 14;
    double gamma = 1.0 / 7;
    double prob_gamma = 1.0 - exp(-gamma / zeta); //probability of transition from I to R at each time step
    double beta[nvar], max_eval;
    double prop_ci[ngr];
    double prop_preg[ngr];
    double prop_im[ngr];
    vector<int> doses;
    vector<vector<double>> sus(nsim, vector<double>(ngr));
    vector<vector<double>> cm(ngr, vector<double>(ngr));
    vector<vector<int>> vax_strategy(ngr, vector<int>(3));

    //epidemiological state variables
    int N[ngr], S[ngr], C[ngr], I[ngr][nvar], R[ngr][nvar];
    int V0mat[ngr][delay1 * zeta], V1mat[ngr][delay2 * zeta];
    int V0[ngr], V1[ngr], V2[ngr];

    double foi[ngr][nvar];
    int newI[ngr][nvar], Rv[ngr];                                   //counter of daily number of new infections by age
    int newI_V0[ngr][nvar], newI_V1[ngr][nvar], newI_V2[ngr][nvar]; //counter of number of infected people by age with vaccine administration at each day
    int nS_I[ngr][nvar], nC_I[ngr][nvar], nI_R[ngr][nvar];
    int tot_nS_I[ngr], tot_nC_I[ngr];
    int nV0_Imat[ngr][delay1 * zeta][nvar]; //number of people progressing from V0 to I by age group and time step
    int nV1_Imat[ngr][delay2 * zeta][nvar]; //number of people progressing from V1 to I by age group and time step
    int tot_nV0_Imat[ngr][delay1 * zeta], tot_nV1_Imat[ngr][delay2 * zeta];
    int nV0_I[ngr][nvar], nV1_I[ngr][nvar], nV2_I[ngr][nvar];
    int tot_nV0_I[ngr], tot_nV1_I[ngr], tot_nV2_I[ngr];
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
    doses = readVector(path, doses);                                  //daily number of doses administered
    sus = readMatrix("data/susceptibility", nsim, ngr, sus);          //relative susceptibility by age
    vax_strategy = readMatrix("data/strategy", ngr, 3, vax_strategy); //vaccination strategy
    //output file
    char outfile[100];
    sprintf(outfile, "%s/newI_doses_R0-%3.1lf-%3.1lf-%3.1lf-%3.1lf-%3.1lf_ve2-%.3lf-%.3lf-%.3lf-%.3lf-%.3lf_init-f-%.3lf-%.3lf-%.3lf-%.3lf-%.3lf_Tonset-%d_strategy-%d_nsim-%d.txt",
            dirName, R0[0], R0[1], R0[2], R0[3], R0[4], ve2[0], ve2[1], ve2[2], ve2[3], ve2[4], init_f[0], init_f[1], init_f[2], init_f[3], init_f[4], Tonset, strategy, nsim);
    FILE *fpout;
    fpout = fopen(outfile, "w");

    //initialize gsl random generator
    // unsigned int seed = time(NULL);
    unsigned int seed = 12306;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);

    for (int sim = 0; sim < nsim; sim++) //contact matrix
    {
        sprintf(path, "data/cm/contact_matrix%d", sim + 1);
        cm = readMatrix(path, ngr, ngr, cm);

        for (int i = 0; i < ngr; i++)
        {
            for (int j = 0; j < ngr; j++)
            {
                cm[i][j] = cm[i][j] * sus[sim][i];
            }
        }

        max_eval = max_eigenval(cm, ngr);
        for (int i = 0; i < nvar; i++)
        {
            beta[i] = R0[i] * gamma / max_eval;
        }

        //initialize state variables
        for (int i = 0; i < ngr; i++)
        {
            S[i] = round(N[i] * (1.0 - prop_im[i] - prop_ci[i] - prop_preg[i])); //susceptibles by age who can be vaccinated
            C[i] = round(N[i] * (prop_ci[i] + prop_preg[i]));                    //contraindications and pregnant women by age who cannot be vaccinated
            for (int j = 0; j < nvar; j++)
            {
                I[i][j] = 0;
                R[i][j] = 0;
            }

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

            for (int j = 0; j < nvar; j++)
            {
                newI[i][j] = 0;
                newI_V0[i][j] = 0;
                newI_V1[i][j] = 0;
                newI_V2[i][j] = 0;
            }

            Rv[i] = 0;

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
                else
                { //vaccinate 18-59 and 60+ since 2021-03-29
                    gr_vax_idx = 1;
                }
            }
            else if (strategy == 2)
            {
                if (Td < 119)
                { //first vaccinate 18-59 between 2020-11-30 and 2021-03-29
                    gr_vax_idx = 0;
                }
                else if (Td < 275)
                { //vaccinate 18-59 and 60+ between 2021-03-29 and 2021-09-01
                    gr_vax_idx = 1;
                }
                else
                {
                    gr_vax_idx = 2;
                }
            }

            if (t == Tonset * zeta) //Tonset == 0, i.e., 2020-11-30
            {
                //seed initial infections
                for (int k = 0; k < nvar; k++)
                {
                    for (int i = 0; i < ni_vec[k]; i++)
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
                                I[q][k]++;
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
                                I[q][k]++;
                            }
                        }
                    }
                }
            }

            //compute age-dependent force of infection
            for (int i = 0; i < ngr; i++)
            {
                tot_nS_I[i] = 0;
                tot_nC_I[i] = 0;
                tot_nV0_I[i] = 0;
                tot_nV1_I[i] = 0;
                tot_nV2_I[i] = 0;
                for (int k = 0; k < nvar; k++)
                {
                    foi[i][k] = 0.0;
                    for (int j = 0; j < ngr; j++)
                    {
                        foi[i][k] += beta[k] * cm[i][j] * (double)I[j][k] / (double)N[j];
                    }

                    double prob_foi = 1.0 - exp(-foi[i][k] / zeta);

                    //epidemiological transitions
                    //new infections
                    nS_I[i][k] = gsl_ran_binomial(rng, prob_foi, S[i]);
                    nC_I[i][k] = gsl_ran_binomial(rng, prob_foi, C[i]);
                    //among vaccinated V0
                    nV0_I[i][k] = 0;
                    for (int j = 0; j < delay1 * zeta; j++)
                    {
                        nV0_Imat[i][j][k] = gsl_ran_binomial(rng, prob_foi, V0mat[i][j]);
                        nV0_I[i][k] += nV0_Imat[i][j][k];
                    }
                    //among vaccinated V1
                    nV1_I[i][k] = 0;
                    for (int j = 0; j < delay2 * zeta; j++)
                    {
                        nV1_Imat[i][j][k] = gsl_ran_binomial(rng, prob_foi, V1mat[i][j]);
                        nV1_I[i][k] += nV1_Imat[i][j][k];
                    }
                    //among vaccinated V2
                    nV2_I[i][k] = gsl_ran_binomial(rng, prob_foi, V2[i]);
                    //progression from infection to recovered
                    nI_R[i][k] = gsl_ran_binomial(rng, prob_gamma, I[i][k]);

                    tot_nS_I[i] += nS_I[i][k];
                    tot_nC_I[i] += nC_I[i][k];
                    tot_nV0_I[i] += nV0_I[i][k];
                    tot_nV1_I[i] += nV1_I[i][k];
                    tot_nV2_I[i] += nV2_I[i][k];
                }

                //update state variables due to transmission
                S[i] -= tot_nS_I[i];
                C[i] -= tot_nC_I[i];
                V0[i] -= tot_nV0_I[i];
                V1[i] -= tot_nV1_I[i];
                V2[i] -= tot_nV2_I[i];
                for (int k = 0; k < nvar; k++)
                {
                    I[i][k] += nS_I[i][k] + nC_I[i][k] + nV0_I[i][k] + nV1_I[i][k] + nV2_I[i][k] - nI_R[i][k];
                    R[i][k] += nI_R[i][k];
                }
                //update V0mat and V1mat
                for (int j = 0; j < delay1 * zeta; j++)
                {
                    tot_nV0_Imat[i][j] = 0;
                    for (int k = 0; k < nvar; k++)
                    {
                        tot_nV0_Imat[i][j] += nV0_Imat[i][j][k];
                    }
                    V0mat[i][j] -= tot_nV0_Imat[i][j];
                }

                for (int j = 0; j < delay2 * zeta; j++)
                {
                    tot_nV1_Imat[i][j] = 0;
                    for (int k = 0; k < nvar; k++)
                    {
                        tot_nV1_Imat[i][j] += nV1_Imat[i][j][k];
                    }
                    V1mat[i][j] -= tot_nV1_Imat[i][j];
                }

                for (int k = 0; k < nvar; k++)
                {
                    newI[i][k] += nS_I[i][k] + nC_I[i][k] + nV0_I[i][k] + nV1_I[i][k] + nV2_I[i][k];
                    newI_V0[i][k] += nV0_I[i][k];
                    newI_V1[i][k] += nV1_I[i][k];
                    newI_V2[i][k] += nV2_I[i][k];
                }
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
                double mean_ve2 = 0.0;
                double numerator = 0.0;
                double denominator = 0.0;
                for (int k = 0; k < nvar; k++)
                {
                    numerator += foi[i][k] * ve2[k];
                    denominator += foi[i][k];
                }

                if (denominator > 0) // use average ve2 weighted by foi if foi is not equal to 0
                {
                    mean_ve2 = numerator / denominator;
                }
                else // use average ve2 weighted by initial fractions
                {
                    for (int k = 0; k < nvar; k++)
                    {
                        mean_ve2 += init_f[k] * ve2[k];
                    }
                }

                double mean_ve1 = 0.67 / 0.8 * mean_ve2;

                S[i] -= nS_V0[i];
                V0[i] = V0[i] + nS_V0[i] - nV0_V1[i];
                V1[i] = V1[i] + round((1 - mean_ve1) * nV0_V1[i]) - nV1_V2[i];
                V2[i] = V2[i] + round((1 - (mean_ve2 - mean_ve1)) * nV1_V2[i]);
                Rv[i] = Rv[i] + round(mean_ve1 * nV0_V1[i]) + round((mean_ve2 - mean_ve1) * nV1_V2[i]);

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
                V1mat[i][0] = round((1 - mean_ve1) * nV0_V1[i]);
            }

            //check state variables
            for (int i = 0; i < ngr; i++)
            {
                if (S[i] < 0 || C[i] < 0 || V0[i] < 0 || V1[i] < 0 || V2[i] < 0)
                {
                    cout << "Negative population in S/C/V0/V1/V2" << endl;
                    cout << "t = " << t << "\t" << S[i] << "\t" << C[i] << "\t" << V0[i] << "\t" << V1[i] << "\t" << V2[i] << endl;
                }

                for (int k = 0; k < nvar; k++)
                {
                    if (I[i][k] < 0 || R[i][k] < 0)
                    {
                        cout << "Negative population in I/R" << endl;
                        cout << "t = " << t << "\t" << I[i][k] << "\t" << R[i][k] << endl;
                    }
                }
            }

            //output daily number of new infections, 1st dose, and 2nd dose by age at each day
            if (t % zeta == 0)
            {
                if (sim == 0 && t == 0) //header
                {
                    fprintf(fpout, "sim t");
                    for (int k = 0; k < nvar; k++)
                    {
                        for (int i = 0; i < ngr; i++)
                        {
                            fprintf(fpout, " %s_newI.%d", variants[k], i + 1);
                        }

                        for (int i = 0; i < ngr; i++)
                        {
                            fprintf(fpout, " %s_newI_V0.%d", variants[k], i + 1);
                        }

                        for (int i = 0; i < ngr; i++)
                        {
                            fprintf(fpout, " %s_newI_V1.%d", variants[k], i + 1);
                        }

                        for (int i = 0; i < ngr; i++)
                        {
                            fprintf(fpout, " %s_newI_V2.%d", variants[k], i + 1);
                        }
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
                for (int k = 0; k < nvar; k++)
                {
                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " %d", newI[i][k]);
                        newI[i][k] = 0;
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " %d", newI_V0[i][k]);
                        newI_V0[i][k] = 0;
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " %d", newI_V1[i][k]);
                        newI_V1[i][k] = 0;
                    }

                    for (int i = 0; i < ngr; i++)
                    {
                        fprintf(fpout, " %d", newI_V2[i][k]);
                        newI_V2[i][k] = 0;
                    }
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
    cout << dirName << endl;
    cout << Tonset << endl;
    cout << strategy << endl;
    cout << capacity << endl;
    cout << nsim << endl;

    for (int k = 0; k < nvar; k++)
    {
        cout << ve2[k] << endl;
        cout << ve1[k] << endl;
    }

    for (int k = 0; k < nvar; k++)
    {
        cout << R0[k] << endl;
    }

    for (int k = 0; k < nvar; k++)
    {
        cout << init_f[k] << endl;
    }

    for (int k = 0; k < nvar; k++)
    {
        cout << variants[k] << endl;
    }

    return 0;
}
