#include "observable.hpp"

int main()
{
    // 4 input parameters needed:
    std::vector<double> phase_uncoupled; // phase shifts for uncoupled channels (in rad).
    std::vector<double> phase_coupled;   // phase shifts for coupled channels (in rad).
    double qq;                           // C.M. relative momentum (in GeV/c).
    int J_max;                           // max J in calculating observables.

    // [exmaple]
    // uncoupled phifts are in order:
    // 1S0, 3P0, 1P1, 3P1, 1D2, 3D2, 1F3, 3F3, 1G4, 3G4, 1H5, 3H5, 1I6, 3I6, 1J7, 3J7, 1K8, 3K8.
    phase_uncoupled = {0.641476, 0.185156, -0.186337, -0.162973, 0.0362048, 0.186437, -0.0238026, -0.0133618, 0.00323008, 0.0173149, -0.00406197, -0.00179301, 0.000479861, 0.00247554, -0.000711222, -0.000261119, 7.7545e-05, 0.000403293};
    // coupled phifts are in order:
    // 3S1, 3D1, E1, 3P2, 3F2, E2, 3D3, 3G3, E3, 3F4, 3H4, E4, 3G5, 3I5, E5, 3H6, 3J6, E6, 3I7, 3K7, E7, 3J8, 3L8, E8.
    phase_coupled = {1.01636, -0.135614, 0.0337778, 0.124484, 0.00750427, -0.0328193, 0.00732112, -0.00653954, 0.0356426, 0.00290842, 0.000559291, -0.00400598, -0.00130156, -0.000641792, 0.00516916, 0.000169245, 5.82765e-05, -0.000576849, -0.000187605, -8.09484e-05, 0.000847567, 1.81514e-05, 7.08604e-06, -8.96401e-05};
    // Tlab = 60 MeV, gives p_{cm} = 0.167774 GeV.
    qq = 0.167774;
    // max J in calculating observables, this value can be smaller than J_max in phase shifts.
    J_max = 8;

    observable::spin_observables(qq, J_max, phase_uncoupled, phase_coupled);
}