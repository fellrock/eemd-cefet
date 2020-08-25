#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>
#include <map>

using namespace std;

void printInFile(const string method,
    const int numberOfIterations,
    const double *mse,
    const double *msd) {

    ofstream myfile1;
    string mseName = method + "_MSE.txt";
    myfile1.open(mseName, ofstream::trunc);

    ofstream myfile2;
    string msdName = method + "_MSD.txt";
    myfile2.open(msdName, ofstream::trunc);

    for (unsigned int k = 0; k < numberOfIterations; ++k) {
        myfile1 << mse[k] << endl;
        myfile2 << msd[k] << endl;
    }

    myfile1.close();
    myfile2.close();
}

void MGBlockSparsity(const int N, double *w) {

    const double p1 = 0.990;
    const double p2 = 0.955;
    const double sigmas2 = 1;

    default_random_engine generator;
    generator.seed( time( NULL ) );

    uniform_int_distribution<>  i_dist(0, 1);
    normal_distribution<double> s_dist(0.0, sqrt(sigmas2));

    int currentState = 1; // Coeficiente nao nulo
    int prob = i_dist ( generator );
    if (prob < ( 1 - p2 ) / ( 2 - p1 - p2 )) {
        currentState = 0; // Coeficiente nulo
    }

    for (int i = 0; i < N; ++i) {

        if (currentState == 1) {
            w[i] = s_dist ( generator );
        } else {
            w[i] = 0.0;
        }

        prob = i_dist ( generator );
        if (currentState == 0) { // Se coeficiente atual for nulo
            if (prob > p1) {
                currentState = 1;
            }
        } else {                 // Se coeficiente atual nao for nulo
            if (prob > p2) {
                currentState = 0;
            }
        }
    }
}

void run(const int N,
    const unsigned int numberOfRepeats,
    const double sigmau2,
    const double sigmanu2,
    const double sigmaeta2,
    const double betaLMS,
    const int numberOfIterations,
    const int P,
    const int B) {

    default_random_engine generator;
    generator.seed( time( NULL ) );

    normal_distribution<double> u_dist(0.0,   sqrt(sigmau2));
    normal_distribution<double> nu_dist(0.0,  sqrt(sigmanu2));
    normal_distribution<double> eta_dist(0.0, sqrt(sigmaeta2));

    const unsigned int tam = numberOfIterations + N - 1;

    double *LMSmse         = (double*)calloc(numberOfIterations, sizeof(double));
    double *LMSmsd         = (double*)calloc(numberOfIterations, sizeof(double));

    double *u   = (double*)calloc(tam, sizeof(double));
    double *nu  = (double*)calloc(tam, sizeof(double));
    double *eta = (double*)calloc(tam, sizeof(double));
    double *x   = (double*)calloc(tam, sizeof(double));
    double *d   = (double*)calloc(tam, sizeof(double));

    double *w              = (double*)calloc(N, sizeof(double));
    double *LMSwk          = (double*)calloc(N, sizeof(double));

    double *listNorm2 = (double*)calloc(N, sizeof(double));
    vector< double/* norm */> listSortByNorm2;
    int finalBlockPos = (B * P);

    //cout << "repeat\ttime" << endl;
    for (unsigned int repeat = 0; repeat < numberOfRepeats; ++repeat) {

        //const clock_t begin_time = clock();

        // Define w* - using MG (Block-Sparsity) model
        MGBlockSparsity(N, w);

        for (unsigned int k = 0; k < tam; ++k) {
            u[k]   = u_dist  ( generator );
            nu[k]  = nu_dist ( generator );
            eta[k] = eta_dist( generator );
            x[k]   = u[k] + eta[k];
            d[k]   = 0.0;
        }
        
        for (unsigned int k = N - 1; k < tam; ++k) {
            for (int i = 0, it = k; i < N; ++i, --it) {
                d[k] += w[i] * u[it];
            }
            d[k] += nu[k];
        }

        for (int i = 0; i < N; ++i) {
            LMSwk[i]         = 0.0;
        }

        for (unsigned int k = N - 1; k < tam; ++k) {
            
            // Define 'yk'
            double LMSyk         = 0.0;
            double BCLMSyk       = 0.0;
            double NLMSyk        = 0.0;
            double BSLMSyk       = 0.0;
            double BSNLMSyk      = 0.0;
            double BCBSNLMSyk    = 0.0;
            double SPUBCBSNLMSyk = 0.0;
            for (int i = 0, it = k; i < N; ++i, --it) {
                LMSyk         += LMSwk[i]         * x[it];
            }

            // Define 'ek'
            const double LMSek         = d[k] - LMSyk;
            const double BCLMSek       = d[k] - BCLMSyk;
            const double NLMSek        = d[k] - NLMSyk;
            const double BSLMSek       = d[k] - BSLMSyk;
            const double BSNLMSek      = d[k] - BSNLMSyk;
            const double BCBSNLMSek    = d[k] - BCBSNLMSyk;
            const double SPUBCBSNLMSek = d[k] - SPUBCBSNLMSyk;


            // Update 'wk'
            double norm_x2            = 0.0;
            double block_norm_x2      = 0.0;
            double LMSnorm_w2         = 0.0;
            listSortByNorm2.clear();
            int idx = 0;
            for (int i = 0, it = k, l = 1; i < N; ++i, --it, ++l) {

                norm_x2       += (x[it] * x[it]);
                block_norm_x2 += (x[it] * x[it]);
                
                // LMS
                LMSwk[i] += (betaLMS * LMSek * x[it]);
                const double LMSwki = LMSwk[i] - w[i];
                LMSnorm_w2 += (LMSwki * LMSwki);

                // SPU-BC-BS-NLMS : norm by block (each one with size P)
                // SPU-BC-BS-NLMS : calculates the norm for n-blocks, to consider the B-greatestes only
                if (l % P == 0) {

                    listSortByNorm2.push_back(block_norm_x2);
                    for (; l > 0; --l) {
                        listNorm2[idx++] = block_norm_x2;
                    }
                    block_norm_x2 = 0.0;
                }
            }

            sort(listSortByNorm2.begin(), listSortByNorm2.end());

            /* FOR DEBUG ONLY

            cout << listSortByNorm2.size() << endl;
            for (int j = 1, i = listSortByNorm2.size() - finalBlockPos; i < listSortByNorm2.size(); ++i, ++j) {
                cout << j << " .. " << listSortByNorm2[i] << endl;
            }
            cout << "minBlockNorm2 = " << listSortByNorm2[listSortByNorm2.size() - finalBlockPos] << endl;
            cout << "maxBlockNorm2 = " << listSortByNorm2[listSortByNorm2.size() - 1] << endl;
            */

            double minBlockNorm2 = listSortByNorm2[listSortByNorm2.size() - finalBlockPos];
            LMSmse[k - N + 1]         += (LMSek * LMSek)                 / numberOfRepeats;
            LMSmsd[k - N + 1]         += LMSnorm_w2                      / numberOfRepeats;
        }

        //cout << repeat << "\t" << float( clock () - begin_time ) / CLOCKS_PER_SEC << endl;
    }

    printInFile("out-LMS",            numberOfIterations, LMSmse,         LMSmsd);

    free(LMSmse);
    free(LMSmsd);
    free(LMSwk);

    free(u);
    free(nu);
    free(eta);
    free(x);
    free(d);
    free(w);
}

int main() {

    // General variables
    const int N            = 800;
    const double sigmau2   = 1;
    const double sigmaeta2 = pow(10, -6);
    const double sigmanu2  = pow(10, -3);
    const double sigmax2   = sigmau2 + sigmaeta2;
    const int numberOfRepeats = 200;
    const int numberOfIterations = 24000;

    // LMS variables
    const double betaLMS = 0.5 / (N * sigmax2);

    const int    P         = 4; // use to create a matrix with n-blocks (each one with size P)


    // SPU-BC-BS-NLMS variables
    //const int P = 4; // use to create a matrix with n-blocks (each one with size P)
    const int B = 4;   // after calculates the norm for n-blocks, considerates the B-greatestes only.

    const clock_t begin_time = clock();
    run(N, numberOfRepeats, sigmau2, sigmanu2, sigmaeta2, betaLMS, numberOfIterations, P, B);
    cout << "Time : " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds." << endl;

    return EXIT_SUCCESS;
}
