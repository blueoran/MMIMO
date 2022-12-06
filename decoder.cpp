#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstring>
#include <complex>
#include <vector>

using namespace std;

complex<double> **Inverse(complex<double> **M, int n)
{
    complex<double> **A = new complex<double> *[n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new complex<double>[n];
    }
    complex<double> **B = new complex<double> *[n];
    for (int i = 0; i < n; i++)
    {
        B[i] = new complex<double>[n];
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = M[i][j];
            B[i][j] = 0;
        }
    }
    for (int i = 0; i < n; i++)
    {
        B[i][i] = 1;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                complex<double> temp = A[j][i] / A[i][i];
                for (int k = 0; k < n; k++)
                {
                    A[j][k] -= temp * A[i][k];
                    B[j][k] -= temp * B[i][k];
                }
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            B[i][j] /= A[i][i];
        }
    }
    return B;
}

complex<double> *single_Decoder(int mod_order, int num_sender, int num_receiver, complex<double> **H, complex<double> *Y)
{
    complex<double> *X = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++)
    {
        X[i] = 0;
    }
    complex<double> **Ht = new complex<double> *[num_sender];
    for (int i = 0; i < num_sender; i++)
    {
        Ht[i] = new complex<double>[num_receiver];
    }
    for (int i = 0; i < num_sender; i++)
    {
        for (int j = 0; j < num_receiver; j++)
        {
            Ht[i][j] = H[j][i];
        }
    }
    complex<double> **HtH = new complex<double> *[num_sender];
    for (int i = 0; i < num_sender; i++)
    {
        HtH[i] = new complex<double>[num_sender];
    }
    for (int i = 0; i < num_sender; i++)
    {
        for (int j = 0; j < num_sender; j++)
        {
            HtH[i][j] = 0;
            for (int k = 0; k < num_receiver; k++)
            {
                HtH[i][j] += Ht[i][k] * H[k][j];
            }
        }
    }
    complex<double> *HtY = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++)
    {
        HtY[i] = 0;
        for (int j = 0; j < num_receiver; j++)
        {
            HtY[i] += Ht[i][j] * Y[j];
        }
    }
    complex<double> **HtH_inv = new complex<double> *[num_sender];
    for (int i = 0; i < num_sender; i++)
    {
        HtH_inv[i] = new complex<double>[num_sender];
    }
    HtH_inv = Inverse(HtH, num_sender);
    for (int i = 0; i < num_sender; i++)
    {
        for (int j = 0; j < num_sender; j++)
        {
            X[i] += HtH_inv[i][j] * HtY[j];
        }
    }
    return X;
}

complex<double> **Decoder(int mod_order, int num_sender, int num_receiver, int num_ofdm_sym, complex<double> **H, complex<double> **Y)
{
    complex<double> **X = new complex<double> *[num_ofdm_sym];

    for (int i = 0; i < num_ofdm_sym; i++)
    {
        X[i] = single_Decoder(mod_order, num_sender, num_receiver, H, Y[i]);
    }
    return X;
}
