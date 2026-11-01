#include "CaseConfig.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
// Egik sok + Prandtl-Meyer genislemesi: 2D, inviscid, sifir hucum acisi.
double thetaBeta(double beta,double M) {
    const double g=diamond::gamma;
    return std::atan(2/std::tan(beta)*(M*M*std::pow(std::sin(beta),2)-1)/(M*M*(g+std::cos(2*beta))+2));
}
double nu(double M) {
    const double g=diamond::gamma;
    return std::sqrt((g+1)/(g-1))*std::atan(std::sqrt((g-1)/(g+1)*(M*M-1)))-std::atan(std::sqrt(M*M-1));
}
int main() {
    const auto c=diamond::conditions();
    const double g=diamond::gamma, M=diamond::mach, pi=std::acos(-1.0);
    const double theta=std::atan(2*diamond::halfThickness/diamond::chord);
    double lo=std::asin(1/M),hi=lo;
    // Zayif sok dalini ilk kokte sec; ayrilmis sok icin hata ver.
    while(hi<pi/2 && thetaBeta(hi,M)<theta) hi+=0.0001;
    if(hi>=pi/2) {std::cerr<<"Attached weak shock bulunamadi\n";return 1;}
    for(int i=0;i<100;++i){double mid=(lo+hi)/2;if(thetaBeta(mid,M)<theta)lo=mid;else hi=mid;}
    const double beta=(lo+hi)/2, Mn=M*std::sin(beta);
    const double pr=1+2*g/(g+1)*(Mn*Mn-1);
    const double rn=(g+1)*Mn*Mn/(2+(g-1)*Mn*Mn);
    const double Mn2=std::sqrt((1+(g-1)/2*Mn*Mn)/(g*Mn*Mn-(g-1)/2));
    const double M2=Mn2/std::sin(beta-theta), target=nu(M2)+2*theta;
    lo=M2;hi=100;
    for(int i=0;i<100;++i){double mid=(lo+hi)/2;if(nu(mid)<target)lo=mid;else hi=mid;}
    const double M3=(lo+hi)/2;
    const double rear=pr*std::pow((1+(g-1)/2*M2*M2)/(1+(g-1)/2*M3*M3),g/(g-1));
    const double q=0.5*c.rho*c.speed*c.speed;
    const double Cd=2*diamond::halfThickness/diamond::chord*(pr-rear)*c.p/q;
    assert(beta>std::asin(1/M) && beta<pi/2 && pr>1 && M2<M && M3>M2);
    assert(std::abs(c.speed-893.931420629)<1e-6);
    std::cout<<std::setprecision(12)<<"THEORETICAL reference, NOT CFD result\n"
      <<"theta_deg="<<theta*180/pi<<"\nbeta_deg="<<beta*180/pi<<"\np_front_over_pinf="<<pr
      <<"\np_rear_over_pinf="<<rear<<"\nT_front_over_Tinf="<<pr/rn<<"\nMach_front="<<M2
      <<"\nMach_rear="<<M3<<"\nCd_pressure="<<Cd<<"\nCl_pressure=0\n";
}
