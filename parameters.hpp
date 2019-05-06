#include <cmath>
#include <string>
#include <cstdlib>

const double PI = 3.141592653589793238462;
const int anzahl = 250; // number of particles involved
const int n=1;
//const int l=n;

long double beta = 1.;//(l*l/n*n); //beta = (L/N)^2
long double time_max = 15. *2.*PI * n*n*n; //how many steps the program iterates
long double delta_t, delta_t_max = 1.e-003; //how much time (arb. units) one step takes
int framerate = 1000; //for output: output is written every (framerate) steps
long double epsilon = beta * 0.001; //error parameter for delta_t
//long double q = 0.99;
long double sigma = 0.01;
//long double distanz = pow(1/((1/q)-1), 1./3) * n*n;
//long double omega = sqrt(1/((1-q)*distanz*distanz*distanz));

bool calcprotons = false;
bool changingtime = true;
bool centerofmass = false;
bool shouldoutput = true;
bool thurukoutput = false;
bool fileoutput = true;
bool cloudoutput = false; //only the cloud as single particles; average is governed by fileoutput
bool rotatingframe = false;

struct particle {
std::string name;
long double r[3];
long double p[3];
long double dummyr[3];
long double dummyp[3];
long double charge;
long double mass;
int index;
bool dispersed;
};



long double com[3];
//long double eps_0 = 8.85418781762e+008;
//long double pi = 3.14159265358979;
long double k = 1;//(4*pi*eps_0);
//long double hbar = 1;//.05457148e-034;
//long double r_0 = 0.52917721092e-010; //Bohr Radius
long double rotmatrix[3][3];
