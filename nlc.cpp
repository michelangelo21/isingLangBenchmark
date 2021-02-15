//Weronika
#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <random>

#define M_PI 3.14159265358979323846

using namespace std;

int main()
{
    ofstream energia("E.txt");
    ofstream neff("neff.txt");

    const int L = 20;
    const int N = L * L;
    const double T = 2;
    const double ksi = 20.0;
    const int dPhi = 10;
    const int mcs = 230000;
    const double n0 = 1.5;
    const double ne = 1.7;

    //wielomian Legendre'a P2
    double P2[361];
    for (int i = 0; i < 361; i++)
    {
        P2[i] = 1.5 * pow(cos((i - 180) * M_PI / 180), 2) - 0.5;
    }

    //natężenie
    const int s = 20;
    const double dE = 0.2;
    double E[s];
    E[0] = 0;
    energia << E[0] << endl;
    for (int i = 1; i < s; i++)
    {
        E[i] = E[i - 1] + dE;
        energia << E[i] << endl;
    }

    //pierwsza konfiguracja kątów
    int phi[L][L];
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            phi[i][j] = 0;
        }
    }

    //tablica najbliższych sąsiadów
    int in[L], ip[L];
    for (int i = 0; i < L; i++)
    {
        in[i] = i + 1;
        ip[i] = i - 1;
    }
    in[L - 1] = 0;
    ip[0] = L - 1;

    //pętla MC
    int phi_trial;
    double U_0, U_new, dU;
    double n_eff = 0;
    double n = 0;

    for (int l = 0; l < s; l++)
    {
        cout << l << endl;
        for (int k = 0; k < mcs; k++)
        {
            for (int i = 1; i < (L - 1); i++)
            {
                for (int j = 0; j < L; j++)
                {
                    //nowe kąty
                    double R = (double)rand() / RAND_MAX;
                    phi_trial = phi[i][j] + int((R - 0.5) * dPhi);

                    if (phi_trial > 90)
                        phi_trial = phi_trial - 180;
                    else if (phi_trial < (-90))
                        phi_trial = phi_trial + 180;

                    //zmiana energii
                    U_0 = -ksi * (P2[phi[i][j] - phi[in[i]][j] + 180] + P2[phi[i][j] - phi[ip[i]][j] + 180] + P2[phi[i][j] - phi[i][in[j]] + 180] + P2[phi[i][j] - phi[i][ip[j]] + 180]) - pow(E[l], 2) * P2[270 - phi[i][j]];
                    U_new = -ksi * (P2[phi_trial - phi[in[i]][j] + 180] + P2[phi_trial - phi[ip[i]][j] + 180] + P2[phi_trial - phi[i][in[j]] + 180] + P2[phi_trial - phi[i][ip[j]] + 180]) - pow(E[l], 2) * P2[270 - phi_trial];

                    dU = U_new - U_0;
                    if (dU < 0)
                        phi[i][j] = phi_trial;
                    else
                    {
                        double p = exp(-dU);
                        double r = (double)rand() / RAND_MAX;
                        if (r <= p)
                            phi[i][j] = phi_trial;
                    }
                }
            }
            if ((k > 30000) && (k % 100 == 0))
            {
                for (int i = 0; i < L; i++)
                {
                    for (int j = 0; j < L; j++)
                    {
                        double Cos = cos(phi[i][j] * M_PI / 180);
                        double Sin = sin(phi[i][j] * M_PI / 180);
                        n = (n0 * ne) / sqrt(pow(n0 * Cos, 2) + pow(ne * Sin, 2));
                        n_eff += n;
                    }
                }
            }
        }
        n_eff = n_eff / 800000;
        neff << n_eff << endl;
        n_eff = 0;
    }



}
