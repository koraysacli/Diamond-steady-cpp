#include "CaseConfig.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
namespace fs=std::filesystem;
using namespace diamond;
using Face=std::array<int,4>;
std::string read(const fs::path& p) { std::ifstream in(p); if(!in) throw std::runtime_error("Dosya bulunamadi: "+p.string()); return {std::istreambuf_iterator<char>(in),{}}; }
int main(int argc,char** argv) try {
    if(argc<2) throw std::runtime_error("generateCase CASE [coarse]");
    const fs::path root=fs::absolute(argv[1]), src=fs::path(__FILE__).parent_path();
    const bool coarse=argc>2 && std::string(argv[2])=="coarse";
    if(fs::exists(root/"constant/polyMesh/points") || fs::exists(root/"processor0"))
        throw std::runtime_error("Mevcut mesh/sonuclarin ustune yazilmaz. Yeni vaka dizini kullanin.");
    for(auto d:{"0","system","constant","src"}) fs::create_directories(root/d);
    if(fs::canonical(src)!=fs::canonical(root/"src")) {
        for(const auto& entry: fs::directory_iterator(src))
            if(entry.is_regular_file()) fs::copy_file(entry.path(),root/"src"/entry.path().filename(),fs::copy_options::overwrite_existing);
        for(const auto& name:{"Allrun","monitor.py","README.md",".gitignore"})
            if(fs::exists(src.parent_path()/name)) fs::copy_file(src.parent_path()/name,root/name,fs::copy_options::overwrite_existing);
    }
    auto write=[&](const std::string& path,const std::string& body,const std::string& cls="dictionary") {
        std::ofstream f(root/path); if(!f) throw std::runtime_error("Yazilamadi: "+path);
        f<<"FoamFile { version 2.0; format ascii; class "<<cls<<"; object "<<fs::path(path).filename().string()<<"; }\n"<<body;
    };
    auto str=[](double x){std::ostringstream s;s<<std::setprecision(15)<<x;return s.str();};
    const auto c=conditions();
    const std::vector<double> xs={xmin,-0.2,0,chord/2,chord,0.6,1.5,xmax};
    std::vector<int> nx={100,150,400,400,250,300,400};
    const int ny=coarse?30:300;
    if(coarse) for(auto& n:nx) n/=10;
    std::vector<std::array<double,3>> vertices;
    std::map<std::array<double,3>,int> lookup;
    std::map<std::string,std::vector<Face>> patches;
    for(auto p:{"inlet","outlet","farfield","diamond","frontAndBack"}) patches[p]={};
    auto vid=[&](double x,double y,double z){std::array<double,3> key{x,y,z}; auto it=lookup.find(key);if(it!=lookup.end()) return it->second;int n=vertices.size();vertices.push_back(key);lookup[key]=n;return n;};
    auto h=[](double x){if(x<=0 || x>=chord) return 0.0; else if(x<=chord/2) return 2*halfThickness*x/chord;else return 2*halfThickness*(chord-x)/chord;};
    std::ostringstream blocks; long cells=0;
    for(int sign:{1,-1}) for(size_t i=0;i+1<xs.size();++i) {
        const double x0=xs[i],x1=xs[i+1];
        std::array<std::array<double,2>,4> pts={{{x0,sign*h(x0)},{x1,sign*h(x1)},{x1,sign*ymax},{x0,sign*ymax}}};
        if(sign<0) {std::swap(pts[0],pts[3]);std::swap(pts[1],pts[2]);}
        std::array<int,8> ids;
        for(int k=0;k<2;++k) for(int j=0;j<4;++j) ids[4*k+j]=vid(pts[j][0],pts[j][1],k*span);
        const double grading=coarse?20.0:500.0;
        blocks<<"hex (";for(int v:ids) blocks<<v<<' '; blocks<<") ("<<nx[i]<<' '<<ny<<" 1) simpleGrading (1 "<<std::setprecision(15)<<(sign>0?grading:1/grading)<<" 1)\n";
        cells+=nx[i]*ny;
        auto face=[&](std::array<int,4> js){Face f;for(int j=0;j<4;++j)f[j]=ids[js[j]];return f;};
        const Face low=face({0,4,5,1}),high=face({3,2,6,7});
        patches["farfield"].push_back(sign>0?high:low);
        if(x0>=0 && x1<=chord) patches["diamond"].push_back(sign>0?low:high);
        if(i==0) patches["inlet"].push_back(face({0,3,7,4}));
        if(i+2==xs.size()) patches["outlet"].push_back(face({1,5,6,2}));
        patches["frontAndBack"].push_back(face({0,3,2,1}));patches["frontAndBack"].push_back(face({4,5,6,7}));
    }
    std::ostringstream mesh;mesh<<std::setprecision(15)<<"convertToMeters 1;\nvertices (\n";
    for(auto v:vertices)mesh<<'('<<v[0]<<' '<<v[1]<<' '<<v[2]<<")\n";
    mesh<<");\nblocks (\n"<<blocks.str()<<");\nedges ();\nboundary (\n";
    for(const auto& [name,faces]:patches){mesh<<name<<" { type "<<(name=="frontAndBack"?"empty":name=="diamond"?"wall":"patch")<<"; faces (\n";for(auto f:faces){mesh<<'(';for(int v:f)mesh<<v<<' ';mesh<<")\n";}mesh<<"); }\n";}
    mesh<<");\nmergePatchPairs ();\n";write("system/blockMeshDict",mesh.str());
    const std::string vel="("+str(c.ux)+" "+str(c.uy)+" 0)";
    for(const std::string field:{"U","p","T"}) {
        const bool u=field=="U";
        const std::string value=u?vel:field=="p"?str(c.p):str(c.T);
        std::string body="dimensions "+std::string(u?"[0 1 -1 0 0 0 0]":field=="p"?"[1 -1 -2 0 0 0 0]":"[0 0 0 1 0 0 0]")+";\ninternalField uniform "+value+";\nboundaryField {\n";
        for(const std::string patch:{"inlet","farfield","outlet","diamond","frontAndBack"}) {
            body+=patch+" { ";
            if(patch=="frontAndBack") body+="type empty;";
            else if(patch=="outlet") body+="type zeroGradient;";
            else if(patch=="diamond") body+=u?"type slip;":"type zeroGradient;";
            else if(u) body+="type codedFixedValue; value uniform "+value+"; name diamondFreestream_"+patch+";\ncodeInclude #{\n#include \"CaseConfig.hpp\"\n#};\ncodeOptions #{ -I$(FOAM_CASE)/src #};\ncode #{\n"+read(src/"inletVelocity.C")+"\n#};";
            else body+="type fixedValue; value uniform "+value+";";
            body+="\n}\n";
        }
        write("0/"+field,body+"}\n",u?"volVectorField":"volScalarField");
    }
    write("constant/physicalProperties","thermoType { type hePsiThermo; mixture pureMixture; transport const; thermo hConst; equationOfState perfectGas; specie specie; energy sensibleInternalEnergy; }\nmixture { specie { molWeight "+str(8314.46261815324/gasR)+"; } thermodynamics { Cp "+str(diamond::gamma*gasR/(diamond::gamma-1))+"; hf 0; } transport { mu 0; Pr 0.72; } }\n");
    write("constant/momentumTransport","simulationType laminar;\nlaminar { model Stokes; }\n");
    write("system/fvSchemes",R"(fluxScheme Kurganov;
ddtSchemes { default localEuler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; reconstruct(rho) Minmod; reconstruct(U) MinmodV; reconstruct(T) Minmod; }
snGradSchemes { default corrected; }
)");
    write("system/fvSolution",R"(solvers {
"(rho|rhoU|rhoE).*" { solver diagonal; }
"U.*" { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-9; relTol 0.01; }
"e.*" { $U; tolerance 1e-10; relTol 0; }
}
PIMPLE { maxCo 0.2; rDeltaTSmoothingCoeff 0.02; }
)");
    write("system/decomposeParDict","numberOfSubdomains "+std::to_string(processes)+";\nmethod scotch;\n");
    std::string diag=R"(steadyDiagnostics {
type coded;
libs ("libutilityFunctionObjects.so");
executeControl timeStep;
executeInterval )"+std::to_string(reportEvery)+R"(;
writeControl none;
codeInclude #{
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "CaseConfig.hpp"
#include <fstream>
#include <iomanip>
#};
codeOptions #{ -I$(FOAM_CASE)/src #};
codeData #{ scalarField previousP_; #};
codeExecute #{
)"+read(src/"diagnostics.C")+"\n#};\n}\n";
    write("system/controlDict","application foamRun;\nsolver shockFluid;\nstartFrom latestTime;\nstartTime 0;\nstopAt endTime;\nendTime "+std::to_string(coarse?1200:maximumIterations)+";\ndeltaT 1;\nadjustTimeStep no;\nwriteControl timeStep;\nwriteInterval "+std::to_string(coarse?400:1000)+";\npurgeWrite 3;\nwriteFormat binary;\nwritePrecision 10;\nwriteCompression off;\ntimeFormat general;\ntimePrecision 12;\nrunTimeModifiable true;\nfunctions {\n"+diag+"}\n");
    std::ofstream(root/"diamond.foam");
    std::ofstream manifest(root/"case-summary.txt");
    manifest<<std::setprecision(12)<<"cells="<<cells<<"\nMach="<<mach<<"\nalpha_deg="<<alphaDeg<<"\nT_K="<<c.T<<"\np_Pa="<<c.p<<"\nrho="<<c.rho<<"\nU="<<c.speed<<"\nreference_area="<<chord*span<<"\nLTS: Time is iteration count, NOT physical seconds.\n";
    std::cout<<"Generated "<<cells<<" cells in "<<root<<'\n';
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
