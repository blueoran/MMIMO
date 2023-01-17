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

int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        cout << "Usage:" << endl
             << argv[0]
             << " <module num> <sim times> <sender num> <receiver num> <add noise>" << endl
             << "    Where <module num> should be in {4, 16, 32, 64}," << endl
             << "    <add noise> should be 0/1, 0 representing no noise." << endl;
        exit(-1);
    }
    srand(time(0));
    int mod_order = atoi(argv[1]);
    int times = atoi(argv[2]);
    int sender = atoi(argv[3]);
    int receiver = atoi(argv[4]);
    int add_noise = atoi(argv[5]);

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

    HouseholderQR<Matrix<complex<double>, Dynamic, Dynamic>> QR(H);
    Matrix<complex<double>, Dynamic, Dynamic> Q, R;
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();

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
        Y = H * S;
        if (add_noise != 0)
            Y += w;
        for (int i = 0; i < receiver; ++i)
        {
            YY[i] = Y(i, 0);
            ww[i] = w(i, 0);
        }
        X = Decoder(mod_order, sender, receiver, 1, HH, &YY, ww);

        int error_num = 0;
        cout << "  send  |  decode" << endl;
        for (int i = 0; i < sender; ++i)
        {
            cout << setw(7) << S(i, 0) << " | " << X[0][i] << endl;
            if (S(i, 0) != X[0][i])
                error_num++;
        }
        cout.width(0);
        cout << "Error symbol num:" << error_num << endl;
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