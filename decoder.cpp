#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <set>

using namespace std;
using namespace Eigen;

// ------------------------------ Symbol Generator
// ------------------------------

complex<double> *gen_symbols(int mod_order);

// ------------------------------ Decoders declare
// ------------------------------

/// @brief Zero Forcing Single Decoder
complex<double> *single_Decoder(int mod_order, int num_sender, int num_receiver,
                                complex<double> **H, complex<double> *Y,
                                complex<double> *w);

/// @brief Sphere Single Decoder
complex<double> *sphere_single_Decoder(int mod_order, int num_sender,
                                       int num_receiver, complex<double> **H,
                                       complex<double> *Y, complex<double> *w);

/// @brief Fixed-Sphere Single Decoder
complex<double> *FSD_single_Decoder(int mod_order, int num_sender,
                                    int num_receiver, complex<double> **H,
                                    complex<double> *Y, complex<double> *w);

/// @brief K-Best Sphere Single Decoder
complex<double> *Kbest_sphere_single_Decoder(int mod_order, int num_sender,
                                    int num_receiver, complex<double> **H,
                                    complex<double> *Y, complex<double> *w);


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
                          complex<double> **Y, complex<double> *w = nullptr) {
    complex<double> **X = new complex<double> *[num_ofdm_sym];

    for (int i = 0; i < num_ofdm_sym; i++) {
#ifdef ZF
        X[i] = single_Decoder(mod_order, num_sender, num_receiver, H, Y[i], w);
#elif defined SP
        X[i] = sphere_single_Decoder(mod_order, num_sender, num_receiver, H,
                                     Y[i], w);
#elif defined FSD
        X[i] = FSD_single_Decoder(mod_order, num_sender, num_receiver, H,
                                     Y[i], w);
#elif defined KSD
        X[i] = Kbest_sphere_single_Decoder(mod_order, num_sender, num_receiver, H, 
                                           Y[i], w);
#else
        X[i] = single_Decoder(mod_order, num_sender, num_receiver, H, Y[i], w);
#endif
    }
    return X;
}

// ------------------------------ K-best ------------------------------

#define bestCount 16

typedef pair<double, vector<int> > setNode;

void Kbest_bfs_search(int mod_order, int num_sender, int num_receiver,
                complex<double> *symbols, complex<double> *y,
                complex<double> *R, double &cur_radius_square,
                complex<double> *X) {

    // use heap to maintain the K-Best results rather than sort after calculating all the values
    // use set instead of priority_queue (set is faster when we have O3)
    // use two sets so we can swap these sets when we are getting to the next layer
    set<setNode> cur, nxt;

    // init
    cur.emplace(0, vector<int>());

    while (cur.size()) { // traverse the states of cur_layer

        // use the naming from the original code to ensure code consistency
        double cur_partial_dis = cur.begin()->first;
        vector<int> cur_s(cur.begin()->second);

        cur.erase(cur.begin());

        int cur_layer = num_sender - cur_s.size();
        
        if (cur_layer) { // not the bottom layer

            for (int i = 0; i < mod_order; ++i) { // use the original code to calculate the distance
                // try for each symbol, calculate the cost of this symbol on this layer
                
                complex<double> temp = 0;
                // if cur_layer > num_receiver, calculation cannot be finished,
                // just let cost = 0
                if (cur_layer <= num_receiver) {
                    // calculating dot_product(R(row[cur_layer]),
                    // current_seleceted_symbol) first, calculate the symbols selected
                    // by upper layer
                    for (int i = cur_layer + 1;
                        (i <= num_receiver) && (i <= num_sender); ++i) {
                        // note: Eigen save matrix as col form, so R(i, j) actually
                        // locates at *(R + j * col_len + i)
                        temp -= (*(R + (i - 1) * num_receiver + (cur_layer - 1))) *
                                symbols[cur_s[num_sender - i]];
                    }
                    // second, calculate this possible symbol[i]
                    temp -= (*(R + (cur_layer - 1) * num_receiver + (cur_layer - 1))) *
                            symbols[i];
                    // temp = y[l] - dot_product
                    temp += y[cur_layer - 1];
                }
                // cost = |temp|^2
                // dis(l-1, i) = dis(l) + cost(i)

                // update the state
                cur_partial_dis += temp.real() * temp.real() +
                                   temp.imag() * temp.imag();
                cur_s.push_back(i); // note that we use 'push_back', so the order of cur_s is different from the original code

                // try to insert the state into the next set
                if (nxt.size() < bestCount || cur_partial_dis < nxt.rbegin()->first) {
                    if (nxt.size() == bestCount)
                        nxt.erase(--nxt.end());
                    nxt.emplace(cur_partial_dis, cur_s);
                }
                
                // revoke the update
                cur_partial_dis -= temp.real() * temp.real() +
                                   temp.imag() * temp.imag();
                cur_s.pop_back();
            }

            // not the bottom layer
            // swap the two sets, means we are getting to the next layer
            if (cur.size() == 0)
                std::swap(cur, nxt);
        }else if (cur_partial_dis < cur_radius_square) { // the bottom layer, update the answer
            cur_radius_square = cur_partial_dis;
            for (int i = 0; i < num_sender; i++)
                X[i] = symbols[cur_s[num_sender - i - 1]]; // note that we use 'push_back', so the order of cur_s is different from the original code
        }
    }

}

complex<double> *Kbest_sphere_single_Decoder(int mod_order, int num_sender,
                                       int num_receiver, complex<double> **H,
                                       complex<double> *Y, complex<double> *w) {
    // The only difference between KSD and SD is the search procedure,
    // for coding simplicity, we copy the old codes.
    
    // get the reference symbols
    complex<double> *symbols = gen_symbols(mod_order);

    complex<double> *X = new complex<double>[num_sender];
    complex<double> *s = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++) X[i] = s[i] = 0;

    Matrix<complex<double>, Dynamic, Dynamic> HH(num_receiver, num_sender);
    for (int i = 0; i < num_sender; ++i)
        for (int j = 0; j < num_receiver; ++j) HH(j, i) = H[j][i];

    // map C array to Eigen Matrix
    Matrix<complex<double>, Dynamic, 1> y =
        Map<Matrix<complex<double>, Dynamic, 1>>(Y, num_receiver, 1);

    // Shapes:
    // H: r * s
    // Q: r * r
    // R: r * s
    // y: r * 1
    // s: s * 1
    // Formula:
    // H = Q * R
    // y_hat = Q^* * y

    // Do QR decomposition
    HouseholderQR<Matrix<complex<double>, Dynamic, Dynamic>> QR(HH);
    Matrix<complex<double>, Dynamic, Dynamic> Q, R;
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();

    Matrix<complex<double>, Dynamic, 1> y_hat(num_receiver);
    y_hat = Q.adjoint() * y;

    // start search
    double radius_square = 1e10;

    Kbest_bfs_search(mod_order, num_sender, num_receiver, symbols, y_hat.data(),
               R.data(), radius_square, X);

    delete[] symbols;
    delete[] s;

    return X;
}

// ------------------------------ FSD ------------------------------

#define FSD_FE_DEPTH 2

/// @brief Fixed-Sphere Decoder DFS.
///
/// Defference between original: search all the first 2 or 3 layers,
/// (Full-Extension, FE)
/// and only search one path to the leaf nodes.
void fsd_dfs_search(int mod_order, int num_sender, int num_receiver,
                    complex<double> *symbols, complex<double> *y,
                    complex<double> *R, int cur_layer, complex<double> *cur_s,
                    double cur_partial_dis, double &cur_radius_square,
                    complex<double> *X) {
    if (cur_layer <= 0) {
        // reach the tree leaf, test if the dis is closer,
        // if so, update answer(X) and dis (radius)
        if (cur_partial_dis < cur_radius_square) {
            cur_radius_square = cur_partial_dis;
            for (int i = 0; i < num_sender; ++i) X[i] = cur_s[i];
        }
        return;
    }

    if (num_sender - cur_layer <= FSD_FE_DEPTH) {
        // At Full-Extension stage, search every symbol.
        for (int i = 0; i < mod_order; ++i) {
            complex<double> temp = 0;
            // if cur_layer > num_receiver, calculation cannot be finished,
            // just let cost = 0
            if (cur_layer <= num_receiver) {
                // calculating dot_product(R(row[cur_layer]),
                // current_seleceted_symbol) first, calculate the symbols
                // selected by upper layer
                for (int i = cur_layer + 1;
                     (i <= num_receiver) && (i <= num_sender); ++i) {
                    // note: Eigen save matrix as col form, so R(i, j) actually
                    // locates at *(R + j * col_len + i)
                    temp -= (*(R + (i - 1) * num_receiver + (cur_layer - 1))) *
                            cur_s[i - 1];
                }
                // second, calculate this possible symbol[i]
                temp -=
                    (*(R + (cur_layer - 1) * num_receiver + (cur_layer - 1))) *
                    symbols[i];
                // temp = y[l] - dot_product
                temp += y[cur_layer - 1];
            }
            double dis = cur_partial_dis + temp.real() * temp.real() +
                         temp.imag() * temp.imag();

            cur_s[cur_layer - 1] = symbols[i];
            fsd_dfs_search(mod_order, num_sender, num_receiver, symbols, y, R,
                           cur_layer - 1, cur_s, dis, cur_radius_square, X);
        }
    } else {
        // Below FE stage, only search the symbol with least partial distance.
        // [WARN] ignore cur_layer > num_receiver: we won't test sender >
        // receiver.
        double min_cost = -1;
        int selected = -1;
        for (int i = 0; i < mod_order; ++i) {
            complex<double> temp = 0;
            // if cur_layer > num_receiver, calculation cannot be finished,
            // just let cost = 0
            if (cur_layer <= num_receiver) {
                // calculating dot_product(R(row[cur_layer]),
                // current_seleceted_symbol) first, calculate the symbols
                // selected by upper layer
                for (int i = cur_layer + 1;
                     (i <= num_receiver) && (i <= num_sender); ++i) {
                    // note: Eigen save matrix as col form, so R(i, j) actually
                    // locates at *(R + j * col_len + i)
                    temp -= (*(R + (i - 1) * num_receiver + (cur_layer - 1))) *
                            cur_s[i - 1];
                }
                // second, calculate this possible symbol[i]
                temp -=
                    (*(R + (cur_layer - 1) * num_receiver + (cur_layer - 1))) *
                    symbols[i];
                // temp = y[l] - dot_product
                temp += y[cur_layer - 1];
            }
            double cost = temp.real() * temp.real() + temp.imag() * temp.imag();
            if (cost < min_cost || min_cost < 0) {
                min_cost = cost;
                selected = i;
            }
        }

        if (cur_partial_dis + min_cost > cur_radius_square) {
            return;
        }

        cur_s[cur_layer - 1] = symbols[selected];
        fsd_dfs_search(mod_order, num_sender, num_receiver, symbols, y, R,
                       cur_layer - 1, cur_s, cur_partial_dis + min_cost,
                       cur_radius_square, X);
    }
}

complex<double> *FSD_single_Decoder(int mod_order, int num_sender,
                                    int num_receiver, complex<double> **H,
                                    complex<double> *Y, complex<double> *w) {
    // The only difference between FSD and SD is the dfs procedure,
    // for coding simplicity, we copy the old codes.

    // get the reference symbols
    complex<double> *symbols = gen_symbols(mod_order);

    complex<double> *X = new complex<double>[num_sender];
    complex<double> *s = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++) X[i] = s[i] = 0;

    Matrix<complex<double>, Dynamic, Dynamic> HH(num_receiver, num_sender);
    for (int i = 0; i < num_sender; ++i)
        for (int j = 0; j < num_receiver; ++j) HH(j, i) = H[j][i];

    // map C array to Eigen Matrix
    Matrix<complex<double>, Dynamic, 1> y =
        Map<Matrix<complex<double>, Dynamic, 1>>(Y, num_receiver, 1);

    // Shapes:
    // H: r * s
    // Q: r * r
    // R: r * s
    // y: r * 1
    // s: s * 1
    // Formula:
    // H = Q * R
    // y_hat = Q^* * y

    // Do QR decomposition
    HouseholderQR<Matrix<complex<double>, Dynamic, Dynamic>> QR(HH);
    Matrix<complex<double>, Dynamic, Dynamic> Q, R;
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();

    Matrix<complex<double>, Dynamic, 1> y_hat(num_receiver);
    y_hat = Q.adjoint() * y;

    // start search
    double radius_square = 1e10;

#ifndef SP_RADIUS_OPT
    radius_square = 1e10;
#else
    if (w == nullptr) {
        radius_square = 1e10;
    } else {
        complex<double> mean = 0;
        double variance = 0;
        for (int i = 0; i < num_receiver; ++i) {
            mean += w[i];
        }
        mean /= num_receiver;
        for (int i = 0; i < num_receiver; ++i) {
            complex<double> sub = w[i] - mean;
            variance += sub.real() * sub.real() + sub.imag() * sub.imag();
        }
        if (num_receiver > 1)
            variance /= (num_receiver - 1);
        else
            variance = 1;
        radius_square = mod_order * mod_order * num_receiver * variance;
        cout << "searching radius^2: " << radius_square << endl;
    }
#endif

    fsd_dfs_search(mod_order, num_sender, num_receiver, symbols, y_hat.data(),
                   R.data(), num_sender, s, 0, radius_square, X);

    delete[] symbols;
    delete[] s;

    return X;
}

// ------------------------------ Sphere ------------------------------

void dfs_search(int mod_order, int num_sender, int num_receiver,
                complex<double> *symbols, complex<double> *y,
                complex<double> *R, int cur_layer, complex<double> *cur_s,
                double cur_partial_dis, double &cur_radius_square,
                complex<double> *X) {
    if (cur_layer <= 0) {
        // reach the tree leaf, test if the dis is closer,
        // if so, update answer(X) and dis (radius)
        if (cur_partial_dis < cur_radius_square) {
            cur_radius_square = cur_partial_dis;
            for (int i = 0; i < num_sender; ++i) X[i] = cur_s[i];
        }
        return;
    }
    pair<int, double> symbol_dis[mod_order];

    for (int i = 0; i < mod_order; ++i) {
        // try for each symbol, calculate the cost of this symbol on this layer
        symbol_dis[i].first = i;
        complex<double> temp = 0;
        // if cur_layer > num_receiver, calculation cannot be finished,
        // just let cost = 0
        if (cur_layer <= num_receiver) {
            // calculating dot_product(R(row[cur_layer]),
            // current_seleceted_symbol) first, calculate the symbols selected
            // by upper layer
            for (int i = cur_layer + 1;
                 (i <= num_receiver) && (i <= num_sender); ++i) {
                // note: Eigen save matrix as col form, so R(i, j) actually
                // locates at *(R + j * col_len + i)
                temp -= (*(R + (i - 1) * num_receiver + (cur_layer - 1))) *
                        cur_s[i - 1];
            }
            // second, calculate this possible symbol[i]
            temp -= (*(R + (cur_layer - 1) * num_receiver + (cur_layer - 1))) *
                    symbols[i];
            // temp = y[l] - dot_product
            temp += y[cur_layer - 1];
        }
        // cost = |temp|^2
        // dis(l-1, i) = dis(l) + cost(i)
        symbol_dis[i].second = temp.real() * temp.real() +
                               temp.imag() * temp.imag() + cur_partial_dis;
    }
    // sort the possible symbols by dis
    sort(symbol_dis, symbol_dis + mod_order,
         [](pair<int, double> a, pair<int, double> b) {
             return a.second < b.second;
         });

    // try them at next layer
    for (int i = 0; i < mod_order; ++i) {
        if (symbol_dis[i].second >= cur_radius_square) {
            // if the dis is not optimal, return (symbols ar sorted by dis)
            return;
        }
        // save currently selected symbol
        cur_s[cur_layer - 1] = symbols[symbol_dis[i].first];
        dfs_search(mod_order, num_sender, num_receiver, symbols, y, R,
                   cur_layer - 1, cur_s, symbol_dis[i].second,
                   cur_radius_square, X);
    }
}

complex<double> *sphere_single_Decoder(int mod_order, int num_sender,
                                       int num_receiver, complex<double> **H,
                                       complex<double> *Y, complex<double> *w) {
    // get the reference symbols
    complex<double> *symbols = gen_symbols(mod_order);

    complex<double> *X = new complex<double>[num_sender];
    complex<double> *s = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++) X[i] = s[i] = 0;

    Matrix<complex<double>, Dynamic, Dynamic> HH(num_receiver, num_sender);
    for (int i = 0; i < num_sender; ++i)
        for (int j = 0; j < num_receiver; ++j) HH(j, i) = H[j][i];

    // map C array to Eigen Matrix
    Matrix<complex<double>, Dynamic, 1> y =
        Map<Matrix<complex<double>, Dynamic, 1>>(Y, num_receiver, 1);

    // Shapes:
    // H: r * s
    // Q: r * r
    // R: r * s
    // y: r * 1
    // s: s * 1
    // Formula:
    // H = Q * R
    // y_hat = Q^* * y

    // Do QR decomposition
    HouseholderQR<Matrix<complex<double>, Dynamic, Dynamic>> QR(HH);
    Matrix<complex<double>, Dynamic, Dynamic> Q, R;
    Q = QR.householderQ();
    R = QR.matrixQR().triangularView<Upper>();

    Matrix<complex<double>, Dynamic, 1> y_hat(num_receiver);
    y_hat = Q.adjoint() * y;

    // start search
    double radius_square = 1e10;

#ifndef SP_RADIUS_OPT
    radius_square = 1e10;
#else
    if (w == nullptr) {
        radius_square = 1e10;
    } else {
        complex<double> mean = 0;
        double variance = 0;
        for (int i = 0; i < num_receiver; ++i) {
            mean += w[i];
        }
        mean /= num_receiver;
        for (int i = 0; i < num_receiver; ++i) {
            complex<double> sub = w[i] - mean;
            variance += sub.real() * sub.real() + sub.imag() * sub.imag();
        }
        if (num_receiver > 1)
            variance /= (num_receiver - 1);
        else
            variance = 1;
        radius_square = mod_order * mod_order * num_receiver * variance;
        cout << "searching radius^2: " << radius_square << endl;
    }
#endif

    dfs_search(mod_order, num_sender, num_receiver, symbols, y_hat.data(),
               R.data(), num_sender, s, 0, radius_square, X);

    delete[] symbols;
    delete[] s;

    return X;
}

// ------------------------------ Zero Forcing ------------------------------

complex<double> **zeroforcing_Inverse(complex<double> **M, int n) {
    complex<double> **A = new complex<double> *[n];
    for (int i = 0; i < n; i++) {
        A[i] = new complex<double>[n];
    }
    complex<double> **B = new complex<double> *[n];
    for (int i = 0; i < n; i++) {
        B[i] = new complex<double>[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = M[i][j];
            B[i][j] = 0;
        }
    }
    for (int i = 0; i < n; i++) {
        B[i][i] = 1;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                complex<double> temp = A[j][i] / A[i][i];
                for (int k = 0; k < n; k++) {
                    A[j][k] -= temp * A[i][k];
                    B[j][k] -= temp * B[i][k];
                }
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[i][j] /= A[i][i];
        }
    }
    return B;
}

complex<double> *single_Decoder(int mod_order, int num_sender, int num_receiver,
                                complex<double> **H, complex<double> *Y,
                                complex<double> *w) {
    complex<double> *X = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++) {
        X[i] = 0;
    }
    complex<double> **Ht = new complex<double> *[num_sender];
    for (int i = 0; i < num_sender; i++) {
        Ht[i] = new complex<double>[num_receiver];
    }
    for (int i = 0; i < num_sender; i++) {
        for (int j = 0; j < num_receiver; j++) {
            Ht[i][j] = H[j][i];
        }
    }
    complex<double> **HtH = new complex<double> *[num_sender];
    for (int i = 0; i < num_sender; i++) {
        HtH[i] = new complex<double>[num_sender];
    }
    for (int i = 0; i < num_sender; i++) {
        for (int j = 0; j < num_sender; j++) {
            HtH[i][j] = 0;
            for (int k = 0; k < num_receiver; k++) {
                HtH[i][j] += Ht[i][k] * H[k][j];
            }
        }
    }
    complex<double> *HtY = new complex<double>[num_sender];
    for (int i = 0; i < num_sender; i++) {
        HtY[i] = 0;
        for (int j = 0; j < num_receiver; j++) {
            HtY[i] += Ht[i][j] * Y[j];
        }
    }
    complex<double> **HtH_inv = new complex<double> *[num_sender];
    for (int i = 0; i < num_sender; i++) {
        HtH_inv[i] = new complex<double>[num_sender];
    }
    HtH_inv = zeroforcing_Inverse(HtH, num_sender);
    for (int i = 0; i < num_sender; i++) {
        for (int j = 0; j < num_sender; j++) {
            X[i] += HtH_inv[i][j] * HtY[j];
        }
        X[i].imag(X[i].imag());
        X[i].real(X[i].real());

        int r_i = (int)floor(X[i].real());
        if (r_i % 2 == 0) r_i++;
        int i_i = (int)floor(X[i].imag());
        if (i_i % 2 == 0) i_i++;
        X[i].real((double)r_i);
        X[i].imag((double)i_i);
    }
    return X;
}

// ------------------------------ Symbol Generator
// ------------------------------

complex<double> *gen_symbols(int mod_order) {
    complex<double> *symbols = new complex<double>[mod_order];
    int num = 0;
    if (mod_order >= 4) {
        for (int i = -1; i <= 1; i += 2)
            for (int j = -1; j <= 1; j += 2)
                symbols[num++] = complex<double>(i, j);
    }
    if (mod_order >= 16) {
        for (int j = -1; j <= 1; j += 2)
            symbols[num++] = complex<double>(-3, j);
        for (int j = -1; j <= 1; j += 2) symbols[num++] = complex<double>(3, j);
        for (int i = -1; i <= 1; i += 2)
            symbols[num++] = complex<double>(i, -3);
        for (int i = -1; i <= 1; i += 2) symbols[num++] = complex<double>(i, 3);
        for (int i = -3; i <= 3; i += 6)
            for (int j = -3; j <= 3; j += 6)
                symbols[num++] = complex<double>(i, j);
    }
    if (mod_order >= 32) {
        for (int j = -3; j <= 3; j += 2)
            symbols[num++] = complex<double>(-5, j);
        for (int j = -3; j <= 3; j += 2) symbols[num++] = complex<double>(5, j);
        for (int i = -3; i <= 3; i += 2)
            symbols[num++] = complex<double>(i, -5);
        for (int i = -3; i <= 3; i += 2) symbols[num++] = complex<double>(i, 5);
    }
    if (mod_order >= 64) {
        for (int i = -5; i <= 5; i += 10)
            for (int j = -5; j <= 5; j += 10)
                symbols[num++] = complex<double>(i, j);
        for (int j = -5; j <= 5; j += 2)
            symbols[num++] = complex<double>(-7, j);
        for (int j = -5; j <= 5; j += 2) symbols[num++] = complex<double>(7, j);
        for (int i = -5; i <= 5; i += 2)
            symbols[num++] = complex<double>(i, -7);
        for (int i = -5; i <= 5; i += 2) symbols[num++] = complex<double>(i, 7);
        for (int i = -7; i <= 7; i += 14)
            for (int j = -7; j <= 7; j += 14)
                symbols[num++] = complex<double>(i, j);
    }
    return symbols;
}