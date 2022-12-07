#include "../decoder.cpp"

complex<double> *C_Decoder(int mod_order, int num_sender, int num_receiver, int num_ofdm_sym, complex<double> **H, complex<double> *Y, complex<double> *w = nullptr)
{
    // std::cerr << "mod_order: " << mod_order << std::endl;
    // std::cerr << "num_sender: " << num_sender << std::endl;
    // std::cerr << "num_receiver: " << num_receiver << std::endl;
    // std::cerr << "num_ofdm_sym: " << num_ofdm_sym << std::endl;
    // std::cerr << "H: " << std::endl;
    // for (int i = 0; i < num_receiver; i++)
    // {
    //     for (int j = 0; j < num_sender; j++)
    //     {
    //         std::cerr << *H[i] << " ";
    //     }
    //     std::cerr << std::endl;
    // }
    // std::cerr << std::endl;
    // std::cerr << "Y: " << std::endl;
    // for (int i = 0; i < num_receiver; i++)
    // {
    //     std::cerr << Y[i] << " ";
    // }
    // std::cerr << std::endl;
    // std::cerr << "w: " << std::endl;
    // for (int i = 0; i < num_receiver; i++)
    // {
    //     std::cerr << w[i] << " ";
    // }
    // std::cerr << std::endl;

    complex<double> **Y2 = new complex<double> *;
    Y2[0] = Y;
    return Decoder(mod_order, num_sender, num_receiver, num_ofdm_sym, H, Y2, w)[0];
}

extern "C"
{
    double *CDecoder(int mod_order, int num_sender, int num_receiver, int num_ofdm_sym, double *H, double *Y, double *w = nullptr)
    {
        double **H_2 = new double *[num_receiver];
        for (int i = 0; i < num_receiver; i++)
        {
            H_2[i] = H + 2 * i * num_receiver;
        }
        complex<double> **H_c = (complex<double> **)H_2;
        complex<double> *Y_c = (complex<double> *)Y;
        complex<double> *w_c = (complex<double> *)w;
        complex<double> *X_c = C_Decoder(mod_order, num_sender, num_receiver, num_ofdm_sym, H_c, Y_c, w_c);

        // for (int i = 0; i < num_receiver; i++)
        // {
        //     std::cerr << Y_c[i] << " ";
        // }
        // std::cerr << std::endl;
        // for (int i = 0; i < num_sender; i++)
        // {
        //     std::cerr << X_c[i] << " ";
        // }
        // std::cerr << std::endl;
        return (double *)X_c;
    }
}
