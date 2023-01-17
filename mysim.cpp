#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>
#include <utility>
#include "decoder.cpp"

using namespace std;
using namespace Eigen;

int main()
{
    srand(time(0));
    int mod_order = 16;
    int times = 1;
    int sender = 8;
    int receiver = 6;
    complex<double> *symbols = gen_symbols(mod_order);
    Matrix<complex<double>, Dynamic, 1> S(sender);
    Matrix<complex<double>, Dynamic, 1> Y(receiver);
    Matrix<complex<double>, Dynamic, 1> w(receiver);
    Matrix<complex<double>, Dynamic, Dynamic> H(receiver, sender);
    complex<double> **X;

    H.setRandom();

    complex<double> **HH = new complex<double> *[receiver];
    for (int i = 0; i < receiver; ++i)
    {
        HH[i] = new complex<double>[sender];
    }

    for (int i = 0; i < sender; ++i)
    {
        for (int j = 0; j < receiver; ++j)
        {
            HH[j][i] = H(j, i);
        }
    }
    complex<double> *YY = new complex<double>[receiver];
    complex<double> *ww = new complex<double>[receiver];

    for (int t = 1; t <= times; ++t)
    {
        cout << "time: " << t << endl;
        w.setRandom();
        for (int i = 0; i < sender; ++i)
        {
            S(i, 0) = symbols[rand() % mod_order];
        }
        Y = w + H * S;
        cout << "S" << endl
             << S << endl
             << "w" << endl
             << w << endl
             << "Y" << endl
             << Y << endl;
        for (int i = 0; i < receiver; ++i)
        {
            YY[i] = Y(i, 0);
            ww[i] = w(i, 0);
        }
        X = Decoder(mod_order, sender, receiver, 1, HH, &YY, ww);

        for (int i = 0; i < sender; ++i)
        {
            cout << S(i, 0) << " " << X[0][i] << endl;
        }
        cout << endl;
    }
    delete[] ww;
    delete[] YY;
    for (int i = 0; i < receiver; ++i)
    {
        delete[] HH[i];
    }
    delete[] HH;

    return 0;
}