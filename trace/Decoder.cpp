#include "/home/steven/Course/Network/MMIMO/decoder.cpp"

extern "C"
{

    complex<double> **Decoder(int mod_order, int num_sender, int num_receiver, int num_ofdm_sym, complex<double> **H, complex<double> **Y)
    {
        }
    std::complex<double> sum_it_cplx(std::complex<double> *array, int size)
    {
        std::complex<double> ret(0., 0.);
        for (int i = 0; i < size; ++i)
        {
            ret += array[i];
        }
        return ret;
    }
}