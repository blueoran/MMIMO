#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;
using namespace Eigen;

// ------------------------------ Symbol Generator ------------------------------

complex<double> *gen_symbols(int mod_order);

// ------------------------------ Decoders declare ------------------------------

/// @brief Zero Forcing Single Decoder
complex<double> *single_Decoder(int mod_order, int num_sender, int num_receiver,
                                complex<double> **H, complex<double> *Y);

/// @brief Sphere Single Decoder
complex<double> *sphere_single_Decoder(int mod_order, int num_sender, int num_receiver,
                                       complex<double> **H, complex<double> *Y);

// ------------------------------ Interface ------------------------------

/// @brief Decoder
/// @param mod_order QFDM scale (per dim)
/// @param num_sender sender antenna num
/// @param num_receiver receiver antenna num
/// @param H signal distortion
/// @param Y received signal
/// @param w noise
/// @return decoded vectors X
complex<double> **Decoder(int mod_order, int num_sender, int num_receiver,
                          int num_ofdm_sym, complex<double> **H,
                          complex<double> **Y, complex<double> *w = nullptr)
{
    complex<double> **X = new complex<double> *[num_ofdm_sym];

    for (int i = 0; i < num_ofdm_sym; i++)
    {
        X[i] = sphere_single_Decoder(mod_order, num_sender, num_receiver, H, Y[i]);
    }
    return X;
}

// ------------------------------ Sphere ------------------------------

void dfs_search(int mod_order, int num_sender, int num_receiver,
                complex<double> *symbols,
                complex<double> *y,
                complex<double> *R,
                int cur_layer,
                complex<double> *cur_s,
                double cur_partial_dis,
                complex<double> cur_temp_accumu_for_cost_calc,
                double &cur_radius_square,
                complex<double> *X)
{
    // cout << cur_layer << " " << cur_radius_square << " " << cur_partial_dis << " " << cur_temp_accumu_for_cost_calc << endl;
    if (cur_layer <= 0)
    {
        if (cur_partial_dis < cur_radius_square)
        {
            cur_radius_square = cur_partial_dis;
            for (int i = 0; i < num_sender; ++i)
                X[i] = cur_s[i];
        }
        return;
    }
    pair<int, double> symbol_dis[mod_order];
    for (int i = 0; i < mod_order; ++i)
    {
        symbol_dis[i].first = i;
        complex<double> temp = cur_temp_accumu_for_cost_calc;
        if (cur_layer <= num_receiver)
            temp += (*(R + (cur_layer - 1) * num_sender + (cur_layer - 1))) * symbols[i];
        temp = -temp;
        if (cur_layer <= num_receiver)
            temp += y[cur_layer - 1];
        symbol_dis[i].second = temp.real() * temp.real() + temp.imag() * temp.imag() + cur_partial_dis;
    }
    sort(symbol_dis, symbol_dis + mod_order, [](pair<int, double> a, pair<int, double> b)
         { return a.second < b.second; });
    for (int i = 0; i < mod_order; ++i)
    {
        if (symbol_dis[i].second >= cur_radius_square)
        {
            return;
        }
        cur_s[cur_layer - 1] = symbols[symbol_dis[i].first];
        complex<double> temp = cur_temp_accumu_for_cost_calc + (*(R + (cur_layer - 1) * num_sender + (cur_layer - 1))) * symbols[cur_layer - 1];
        dfs_search(mod_order, num_sender, num_receiver, symbols, y, R, cur_layer - 1, cur_s, symbol_dis[i].second, temp, cur_radius_square, X);
    }
}

complex<double> *sphere_single_Decoder(int mod_order, int num_sender, int num_receiver,
                                       complex<double> **H, complex<double> *Y)
{
    complex<double> *symbols = gen_symbols(mod_order);

    complex<double> *X = new complex<double>[num_sender];
    complex<double> *s = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++)
        X[i] = s[i] = 0;

    Matrix<complex<double>, Dynamic, Dynamic> HH(num_receiver, num_sender);
    for (int i = 0; i < num_sender; ++i)
        for (int j = 0; j < num_receiver; ++j)
            HH(j, i) = H[j][i];

    Matrix<complex<double>, Dynamic, 1> y = Map<Matrix<complex<double>, Dynamic, 1>>(Y, num_receiver, 1);

    // Shapes:
    // H: r * s
    // Q: r * r
    // R: r * s
    // y: r * 1
    // s: s * 1
    // Formula:
    // H = Q * R
    // y_hat = Q^* * y
    HouseholderQR<Matrix<complex<double>, Dynamic, Dynamic>> QR(HH);
    Matrix<complex<double>, Dynamic, Dynamic> Q, R;
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();

    Matrix<complex<double>, Dynamic, 1> y_hat(num_receiver);
    y_hat = Q.adjoint() * y;

    double radius_square = 10000000000.0;
    dfs_search(mod_order, num_sender, num_receiver, symbols, y_hat.data(), R.data(), num_sender, s, 0, 0, radius_square, X);

    delete[] symbols;
    delete[] s;

    return X;
}

// ------------------------------ Zero Forcing ------------------------------

complex<double> **zeroforcing_Inverse(complex<double> **M, int n)
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

complex<double> *single_Decoder(int mod_order, int num_sender, int num_receiver,
                                complex<double> **H, complex<double> *Y)
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
    HtH_inv = zeroforcing_Inverse(HtH, num_sender);
    for (int i = 0; i < num_sender; i++)
    {
        for (int j = 0; j < num_sender; j++)
        {
            X[i] += HtH_inv[i][j] * HtY[j];
        }
        X[i].imag(X[i].imag());
        X[i].real(X[i].real());
    }
    return X;
}

// ------------------------------ Symbol Generator ------------------------------

complex<double> *gen_symbols(int mod_order)
{
    complex<double> *symbols = new complex<double>[mod_order];
    int num = 0;
    if (mod_order >= 4)
    {
        for (int i = -1; i <= 1; i += 2)
            for (int j = -1; j <= 1; j += 2)
                symbols[num++] = complex<double>(i, j);
    }
    if (mod_order >= 16)
    {
        for (int j = -1; j <= 1; j += 2)
            symbols[num++] = complex<double>(-3, j);
        for (int j = -1; j <= 1; j += 2)
            symbols[num++] = complex<double>(3, j);
        for (int i = -1; i <= 1; i += 2)
            symbols[num++] = complex<double>(i, -3);
        for (int i = -1; i <= 1; i += 2)
            symbols[num++] = complex<double>(i, 3);
        for (int i = -3; i <= 3; i += 6)
            for (int j = -3; j <= 3; j += 6)
                symbols[num++] = complex<double>(i, j);
    }
    if (mod_order >= 32)
    {
        for (int j = -3; j <= 3; j += 2)
            symbols[num++] = complex<double>(-5, j);
        for (int j = -3; j <= 3; j += 2)
            symbols[num++] = complex<double>(5, j);
        for (int i = -3; i <= 3; i += 2)
            symbols[num++] = complex<double>(i, -5);
        for (int i = -3; i <= 3; i += 2)
            symbols[num++] = complex<double>(i, 5);
    }
    if (mod_order >= 64)
    {
        for (int i = -5; i <= 5; i += 10)
            for (int j = -5; j <= 5; j += 10)
                symbols[num++] = complex<double>(i, j);
        for (int j = -5; j <= 5; j += 2)
            symbols[num++] = complex<double>(-7, j);
        for (int j = -5; j <= 5; j += 2)
            symbols[num++] = complex<double>(7, j);
        for (int i = -5; i <= 5; i += 2)
            symbols[num++] = complex<double>(i, -7);
        for (int i = -5; i <= 5; i += 2)
            symbols[num++] = complex<double>(i, 7);
        for (int i = -7; i <= 7; i += 14)
            for (int j = -7; j <= 7; j += 14)
                symbols[num++] = complex<double>(i, j);
    }
    return symbols;
}