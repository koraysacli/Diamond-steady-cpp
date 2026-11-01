#pragma once
#include <cmath>
#include <stdexcept>
namespace diamond {
// User settings: SI units, angles in degrees. U_y = -U*sin(alpha).
constexpr double mach=3.0, altitudeFt=80000.0, alphaDeg=0.0;
constexpr double chord=0.2, halfThickness=0.01, span=0.001;
constexpr double xmin=-1.0, xmax=4.0, ymax=2.5;
constexpr int processes=20, maximumIterations=20000, reportEvery=20;
constexpr double gasR=287.05287, gamma=1.4;
struct Conditions { double T,p,rho,speed,ux,uy; };
inline Conditions conditions() {
    const double z=altitudeFt*0.3048, H=6356766.0*z/(6356766.0+z);
    double T,p;
    // First three layers of the standard atmosphere; error outside this range.
    if (H<0 || H>32000) throw std::runtime_error("Atmosphere model range is 0--32 km");
    if (H<=11000) {
        T=288.15-0.0065*H;
        p=101325*std::pow(T/288.15,9.80665/(gasR*0.0065));
    } else if (H<=20000) {
        T=216.65;
        p=22632.04*std::exp(-9.80665*(H-11000)/(gasR*T));
    } else {
        T=216.65+0.001*(H-20000);
        p=5474.889*std::pow(216.65/T,9.80665/(gasR*0.001));
    }
    const double speed=mach*std::sqrt(gamma*gasR*T);
    const double alpha=alphaDeg*std::acos(-1.0)/180;
    return {T,p,p/(gasR*T),speed,speed*std::cos(alpha),-speed*std::sin(alpha)};
}
}
