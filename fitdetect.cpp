#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <dirent.h>

using namespace std;


int main(int argc, const char *argv[]) {

stringstream outname, inname, inname2;
std::ifstream in, in2;
std::ofstream out;
long double distance, c, l, r_anf;

long double dummy[7];

inname.str("");
inname<<"./"<<argv[1]<<"/fitting_parameters.txt";

outname.str("");
outname<<"./"<<argv[1]<<"/disptime_fit.txt";

in.open(inname.str().c_str());
if (!in.is_open()) {
	cerr << "Couldn't open in!"<<endl;
	return 0;
}

out.open(outname.str().c_str());
if (!out.is_open()) {
	cerr << "Couldn't open out!"<<endl;
	return 0;
}

while (in >> distance >> c >> l) {

	inname2.str("");
	inname2<<"./"<<argv[1]<<"/avg_cloud_average_"<<distance<<".txt";

	in2.open(inname.str().c_str());
	if (!in2.is_open()) {
		cerr << "Couldn't open in2!"<<endl;
		return 0;
	}

	if (in2 >> dummy[0] /*>> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4] >> dummy[5]*/ >> r_anf /*>> dummy[6]*/) {
		out<<distance<<" "<<-log(0.35*r_anf/c)/l<<endl;
	}

	in2.close();
}

in.close();
out.close();
return 0;

}

