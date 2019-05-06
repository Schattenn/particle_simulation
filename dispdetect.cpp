#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <dirent.h>

using namespace std;

int count=100000;
int avgrange=300;

int main(int argc, const char *argv[]) {

stringstream outname, inname, dispname;
std::ifstream input;
std::ofstream out, dispout;
int t;
long double maxima[count], times[count];
long double r[3], p[3], radius, tim;
std::string filename[100];
struct dirent *pdirent;
DIR *pdir;

if (argc<2) {
	cout<<"Blah!"<<endl;
	return 0;
}

for (int n=1; n<argc; n++) {
	pdir=opendir(argv[n]);
	if (pdir==NULL) {
		cout<<"Can't open directory!"<<endl;
		return 1;
	}
	for (int f=0; f<100;) {
		if ((pdirent=readdir(pdir)) != NULL) {
			if (pdirent->d_name[0]=='c') {
				filename[f] = pdirent->d_name;
				cout<<filename[f]<<endl;
				f++;
			}
		} else {
			f+=100;
		}
	}

	dispname.str("");
	dispname << argv[n] << "/avg_disptime.txt";
	dispout.open(dispname.str().c_str());
	if (!dispout.is_open()) {
		cerr << "Couldn't open dispout!"<<endl;
		return 0;
	}

cout<<argv[n]<<endl;
for (int m=0; m<32; m++) {

//input

outname.str("");
inname.str("");

inname << argv[n] << "/"<<filename[m];
outname << argv[n] << "/avg_"<<filename[m];


input.open(inname.str().c_str());
if (!input.is_open()) {
	std::cerr << "Couldn't open input file " << filename[m] << "!";
	return 0;
}

out.open(outname.str().c_str());
if (!out.is_open()) {
	std::cerr << "Couldn't open output file " << m << "!";
	return 0;
}


t=0;
tim=0;
while (input >> r[0] >> r[1] >> r[2] >> p[0] >> p[1] >> p[2] >> radius >> tim) {
	maxima[t]=radius;
	times[t]=tim;
	t++;
	if (tim<=0) {
		cerr<<"HUH?! "<<r[0]<<" "<<r[1]<<" "<<r[2]<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<radius<<" "<<tim<<endl;
	}

}
if (t==0) {
	cerr << "File empty!" << endl;
	return 0;
}
input.close();

//average maxima

t-=avgrange;
long double rad_anf=maxima[0];
bool disp=false;
long double maxavg[t];
string outdist(filename[m].begin()+14, filename[m].end()-4);
for (int i=0; i<t; i++) {
	maxavg[i]=0;
}

for (int i=0; i<t; i++) {
	for (int j=i; j<i+avgrange; j++) { 
		maxavg[i] += maxima[j];
	}
	maxavg[i]/=avgrange;
	out << times[i+avgrange/2] << " " << maxavg[i] << endl;

	if (!disp) {
		if (maxavg[i]<rad_anf*0.35) {
			dispout << outdist << " " << times[i+avgrange/2] << endl;
			disp=true;
		}
	}

}
out.close();
}
dispout.close();
closedir(pdir);
}
return 0;

}

