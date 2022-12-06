#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstring>
#include <complex>
#include <vector>
#include "mex.hpp"
// #include "mex.h"
#include "mexAdapter.hpp"
#include "/home/steven/Course/Network/MMIMO/decoder.cpp"

using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function
{
public:
    /**
     * outputs: x(sent signal)
     *
     * inputs: mod_order(QAM order), num_sender(number of sender's antennas), num_receiver(number of receiver's antennas), num_ofdm_sym, H(channel matrix), y(received signal)
     *
     */
    matlab::data::ArrayFactory factory;
    void operator()(ArgumentList outputs, ArgumentList inputs)
    {
        checkArguments(outputs, inputs);
        const int mod_order = inputs[0][0];
        const int num_sender = inputs[1][0];
        const int num_receiver = inputs[2][0];
        const int num_ofdm_sym = inputs[3][0];
        // std::cerr << "mod_order: " << mod_order << std::endl;
        // std::cerr << "num_sender: " << num_sender << std::endl;
        // std::cerr << "num_receiver: " << num_receiver << std::endl;
        // std::cerr << "num_ofdm_sym: " << num_ofdm_sym << std::endl;

        TypedArray<complex<double>> arrH = std::move(inputs[4]);
        complex<double> *H = arrH.release().get();

        TypedArray<complex<double>> arry = std::move(inputs[5]);
        complex<double> *y = arry.release().get();

        TypedArray<complex<double>> x = factory.createArray<complex<double>>({(unsigned long)num_ofdm_sym, (unsigned long)num_sender});
        // complex<double> *x_ptr = new complex<double>[num_ofdm_sym * num_sender];

        decoder(mod_order, num_sender, num_receiver, num_ofdm_sym, H, y, x);
        // for (int i = 0; i < num_ofdm_sym; i++)
        // {
        //     for (int j = 0; j < num_sender; j++)
        //     {
        //         std::cerr << complex<double>(x[i][j]);
        //     }
        //     std::cerr << '\n';
        // }
        // std::cerr << "\n";
        outputs[0] = std::move(x);
    }
    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs.size() != 6)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("6 inputs required")}));
        }
    }

    void decoder(int mod_order,
                 int num_sender,
                 int num_receiver,
                 int num_ofdm_sym,
                 complex<double> *H,
                 complex<double> *y,
                 TypedArray<complex<double>> &x)
    {
        complex<double> *H_ptr[num_receiver];
        for (int i = 0; i < num_receiver; i++)
        {
            H_ptr[i] = H + i * num_sender;
        }
        complex<double> *y_ptr[num_ofdm_sym];
        for (int i = 0; i < num_ofdm_sym; i++)
        {
            y_ptr[i] = y + i * num_receiver;
        }

        complex<double> **x_ = Decoder(mod_order, num_sender, num_receiver, num_ofdm_sym, H_ptr, y_ptr);
        for (int i = 0; i < num_ofdm_sym; i++)
        {
            for (int j = 0; j < num_sender; j++)
            {
                x[i][j] = x_[i][j];
            }
        }
    }
};
