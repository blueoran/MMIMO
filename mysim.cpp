#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <random>

#include "decoder.cpp"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {
    if (argc != 7) {
        cout << "Usage:" << endl
             << argv[0]
             << " <module num> <sim times> <sender num> <receiver num> <add "
                "noise> <log10-SNR>"
             << endl << endl
             << "    - <module num> should be in {4, 16, 32, 64}," << endl << endl
             << "    - <add noise> should be 0/1, 0 representing no noise." << endl << endl
             << "    - <log10-SNR> should be a real number," << endl
             << "       for instance, if a SNR = 20dB is wanted (S/N = 100), " << endl
             << "       then <log10-SNR> = 2." << endl << endl
             << "    - S/N = E[H(i,i)]/(2 * s^2), " << endl
             << "       where H is the channel matrix, s^2 is the variance of normal noise." << endl;
        exit(-1);
    }

    // precise random seed
    struct timeval time;
    gettimeofday(&time, NULL);
    srand(time.tv_sec * 1000 + time.tv_usec / 1000);

    int mod_order = atoi(argv[1]);
    int times = atoi(argv[2]);
    int sender = atoi(argv[3]);
    int receiver = atoi(argv[4]);
    int add_noise = atoi(argv[5]);
    double SNR = atof(argv[6]);

    complex<double> *symbols = gen_symbols(mod_order);
    Matrix<complex<double>, Dynamic, 1> S(sender);
    Matrix<complex<double>, Dynamic, 1> Y(receiver);
    Matrix<complex<double>, Dynamic, 1> w(receiver);
    Matrix<complex<double>, Dynamic, Dynamic> H(receiver, sender);
    complex<double> **X;

    double pSNR = pow(10.0, SNR);
    normal_distribution<double> normal(0, 1.0);
    uniform_real_distribution<double> uniform(-2 * pSNR, 2 * pSNR);
    default_random_engine e;

    e.seed(rand());

    // initialize
    for (int i = 0; i < receiver; ++i) {
        for(int j = 0;j < sender; ++j) {
            H(i, j) = complex<double>(uniform(e), uniform(e));
        }
    }

    complex<double> **HH = new complex<double> *[receiver];
    for (int i = 0; i < receiver; ++i) {
        HH[i] = new complex<double>[sender];
    }

    for (int i = 0; i < sender; ++i) {
        for (int j = 0; j < receiver; ++j) {
            HH[j][i] = H(j, i);
        }
    }

    // QR, for evaluation
    HouseholderQR<Matrix<complex<double>, Dynamic, Dynamic>> QR(H);
    Matrix<complex<double>, Dynamic, Dynamic> Q, R;
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();

    complex<double> *YY = new complex<double>[receiver];
    complex<double> *ww = new complex<double>[receiver];

    double total_symbol = 0;
    double total_error_symbol = 0;
    double error_time = 0;

    clock_t start, end;
    double send_gen_time = 0;
    double decode_time = 0;

    for (int t = 1; t <= times; ++t) {
        start = clock();

#ifdef VERBOSE
        cout << "time: " << t << endl;
#endif

        e.seed(rand());
        for(int i = 0; i < receiver; ++i) {
            w(i, 0) = normal(e);
        }

        for (int i = 0; i < sender; ++i) {
            // random select symbols to send
            S(i, 0) = symbols[rand() % mod_order];
        }
        Y = H * S;
        if (add_noise != 0) Y += w;
        for (int i = 0; i < receiver; ++i) {
            // deep copy
            YY[i] = Y(i, 0);
            ww[i] = w(i, 0);
        }
        end = clock();
        send_gen_time += double(end - start) / CLOCKS_PER_SEC * 1000;

        start = clock();
        X = Decoder(mod_order, sender, receiver, 1, HH, &YY,
                    add_noise != 0 ? ww : nullptr);
        end = clock();
        decode_time += double(end - start) / CLOCKS_PER_SEC * 1000;

        int error_num = 0;
#ifdef VERBOSE
        cout << "  send  |  decode" << endl;
#endif
        for (int i = 0; i < sender; ++i) {
#ifdef VERBOSE
            cout << setw(7) << S(i, 0) << " | " << X[0][i] << endl;
#endif
            if (abs(S(i, 0).real() - X[0][i].real()) > 1e-6 ||
                abs(S(i, 0).imag() - X[0][i].imag()) > 1e-6) {
                error_num++;
            }
        }
        cout.width(0);
#ifdef VERBOSE
        cout << "Error symbol num:" << error_num << endl;
#endif

        total_symbol += sender;
        total_error_symbol += error_num;
        if (error_num != 0) error_time++;
    }
#ifdef VERBOSE
    cout << endl;
#endif

#ifdef VERBOSE
    // result output
    cout << "---------- SimResult-Verbose ----------" << endl
         << "[Verbose] Modulation Order: " << mod_order << endl
         << "[Verbose] Sender Num: " << sender << endl
         << "[Verbose] Receiver Num:" << receiver << endl
         << "[Verbose] Simulation times: " << times << endl
         << "[Verbose] SNR: " << SNR * 10 << "dB" << endl
         << "[Verbose] Symbol Error Rate: " << total_error_symbol << "/"
         << total_symbol << " = " << total_error_symbol / total_symbol << endl
         << "[Verbose] Error Decoding Times: " << error_time << "/" << times
         << " = " << error_time / times << endl
         << "[Verbose] Decode Time: " << decode_time << "ms" << endl;
#endif

    cout << "SimResult," << mod_order << "," << sender << "," << receiver << ",";
    if (add_noise != 0)
        cout << SNR * 10 << "," << pSNR << ",";
    cout << times << "," << total_error_symbol / total_symbol << ","
         << error_time / times << "," << decode_time << endl;

    delete[] ww;
    delete[] YY;
    for (int i = 0; i < receiver; ++i) {
        delete[] HH[i];
    }
    delete[] HH;

    return 0;
}