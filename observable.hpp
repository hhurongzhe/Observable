#pragma once
#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

#include "lib_define.hpp"
#include "mmatrix.hpp"

namespace observable
{

    constexpr std::complex<double> im_unit{0, 1};  // imaginary unit i.
    constexpr double ppi = 3.14159265358979323846; // Pi.

    std::vector<double> make_table(double xstart, double xend, int inter)
    {
        std::vector<double> temp;
        for (int i = 1; i <= inter; i = i + 1)
        {
            double tempp = xstart + i * (xend - xstart) / inter;
            temp.push_back(tempp);
        }
        return temp;
    }

    void differential_cross_section(double qq, int J_max, std::vector<double> &phase_uncoupled, std::vector<double> &phase_coupled)
    {
        double factor = 0.3894;           // transform factor from [GeV]^(-2) to mb.
        double angle_trans = 180.0 / ppi; // transform for angle to degree unit.

        std::vector<double> temp;                                                 // used to store observables calculated.
        double theta_min = 0;                                                     // theta min value: 0.
        double theta_max = ppi;                                                   // theta max value: Pi.
        int interv = 180;                                                         // theta points in calculation, needs to be an integer.
        std::vector<double> thetalist = make_table(theta_min, theta_max, interv); // a vector storing angles: [0, 1, 2, ..., 180].
        for (int i = 0; i < thetalist.size(); i = i + 1)
        {
            double t = cos(thetalist[i]); // cos(theta)
            std::complex<double> M11 = mmatrix::m11(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> M10 = mmatrix::m10(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> M01 = mmatrix::m01(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> M00 = mmatrix::m00(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> Mpm = mmatrix::mpm(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> Mss = mmatrix::mss(t, qq, J_max, phase_uncoupled, phase_coupled);
            double X = 0.5 * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss)); // differential cross section for this angle.
            temp.push_back(X);
        }

        // print angles in degree:
        for (int i = 0; i < thetalist.size(); i = i + 1)
        {
            std::cout << angle_trans * thetalist[i] << ", ";
        }
        std::cout << "\n\n";
        // print differential cross sections in [mb]:
        for (int i = 0; i < temp.size(); i = i + 1)
        {
            std::cout << factor * temp[i] << ", ";
        }
        std::cout << "\n\n";
    }

    void spin_observables(double qq, int J_max, std::vector<double> &phase_uncoupled, std::vector<double> &phase_coupled)
    {
        double factor = 0.3894;
        double angle_trans = 180.0 / ppi;

        // spin observables,
        // named in accordance with Nijmegen group,
        // can be found on website: https://nn-online.org/NN/help/index.php?page=observables
        std::vector<double> DSG, D, P, A, R, Rp, Axx, Azz, Axz;
        double theta_min = 0;
        double theta_max = ppi;
        int interv = 180;
        std::vector<double> thetalist = make_table(theta_min, theta_max, interv);
        for (int i = 0; i < thetalist.size(); i = i + 1)
        {
            double the = thetalist[i];
            double t = cos(the);
            double tt = sin(the);
            std::complex<double> M11 = mmatrix::m11(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double>
                M10 = mmatrix::m10(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> M01 = mmatrix::m01(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> M00 = mmatrix::m00(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> Mpm = mmatrix::mpm(t, qq, J_max, phase_uncoupled, phase_coupled);
            std::complex<double> Mss = mmatrix::mss(t, qq, J_max, phase_uncoupled, phase_coupled);
            double I0 = 0.5 * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss)); // I0
            DSG.push_back(I0);
            double I01D = 0.25 * norm(M11 + Mpm - Mss) + 0.25 * norm(M11 - Mpm - M00) + 0.5 * norm(M10 + M01); // I0(1-D)
            D.push_back(1 - I01D / I0);
            double I0P = sqrt(2) / 4.0 * real(im_unit * (M10 - M01) * conj(M11 - Mpm + M00)); // I0P
            P.push_back(I0P / I0);
            double I0A = -0.5 *
                         real((M00 + sqrt(2) * (t + 1) / tt * M10) * conj(M11 + Mpm + Mss) -
                              sqrt(2) / tt * (M10 + M01) * conj(M11 + Mpm)) *
                         sin(0.5 * the); // I0A
            A.push_back(I0A / I0);
            double I0R = 0.5 *
                         real((M00 + sqrt(2) * (t - 1) / tt * M10) * conj(M11 + Mpm + Mss) +
                              sqrt(2) / tt * (M10 + M01) * conj(Mss)) *
                         cos(0.5 * the); // I0R
            R.push_back(I0R / I0);
            double I0Rp = 0.5 *
                          real((M00 + sqrt(2) * (t + 1) / tt * M10) * conj(M11 + Mpm + Mss) -
                               sqrt(2) / tt * (M10 + M01) * conj(Mss)) *
                          sin(0.5 * the); // I0Rp
            Rp.push_back(I0Rp / I0);
            double I0Axx = 0.25 * norm(M00) - 0.25 * norm(Mss) - 0.5 * norm(M01) + 0.5 * norm(M10) + real(M11 * conj(Mpm));
            Axx.push_back(I0Axx / I0);
            double I0Azz =
                0.5 * norm(M11) - 0.25 * norm(M00) - 0.25 * norm(Mss) + 0.5 * norm(M01) - 0.5 * norm(M10) + 0.5 * norm(Mpm);
            Azz.push_back(I0Azz / I0);
            double I0Axz = 0.25 * tan(the) * (norm(M11 - Mpm) - norm(M00)) - 0.5 / tan(the) * (norm(M01) - norm(M10));
            Axz.push_back(I0Axz / I0);
        }
        // printing.
        std::cout << "theta:\n";
        for (int i = 0; i < thetalist.size(); i = i + 1)
        {
            std::cout << angle_trans * thetalist[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "DSG:\n";
        for (int i = 0; i < DSG.size(); i = i + 1)
        {
            std::cout << factor * DSG[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "D:\n";
        for (int i = 0; i < D.size(); i = i + 1)
        {
            std::cout << D[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "P:\n";
        for (int i = 0; i < P.size(); i = i + 1)
        {
            std::cout << P[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "A:\n";
        for (int i = 0; i < A.size(); i = i + 1)
        {
            std::cout << A[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "R:\n";
        for (int i = 0; i < R.size(); i = i + 1)
        {
            std::cout << R[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "R':\n";
        for (int i = 0; i < Rp.size(); i = i + 1)
        {
            std::cout << Rp[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "Axx:\n";
        for (int i = 0; i < Axx.size(); i = i + 1)
        {
            std::cout << Axx[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "Azz:\n";
        for (int i = 0; i < Azz.size(); i = i + 1)
        {
            std::cout << Azz[i] << ", ";
        }
        std::cout << "\n\n";

        std::cout << "Axz:\n";
        for (int i = 0; i < Axz.size(); i = i + 1)
        {
            std::cout << Axz[i] << ", ";
        }
        std::cout << "\n\n";
    }

} // namespace observable

#endif // OBSERVABLE_HPP
