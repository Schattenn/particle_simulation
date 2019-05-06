#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <pthread.h>

//#include <direct.h>

#include "moreatoms.hpp"

using namespace std;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

std::string thuruk_output(particle part)
{
	std::stringstream ss;
	ss << '[' << part.r[0] << ',' << part.r[1] << ',' << part.r[2] << ']';
	return ss.str();
}


void output(particle &part, ofstream &file) {   //output in file of properties of part
    file<<"\n"<<part.name<<"\n";
	file<<"\nCoordinates:\n";
	file<<part.r[0] << "\n";
	file<<part.r[1] << "\n";
	file<<part.r[2] << "\n";
	file<<"\nMomentum:\n";
	file<<part.p[0] << "\n";
	file<<part.p[1] << "\n";
	file<<part.p[2] << "\n";
}

void* haupt(void* qq) {

    stringstream ss_q;
    ss_q.str("");
    ss_q<<"./outputfiles/disptime.txt";
    ofstream file_q(ss_q.str().c_str());
    if (!file_q.is_open()) {
        cout<<"File_q not opened.";
        return 0;
    }
    file_q.close();

    long double qqq[20];
    qqq[0] = 0.987;
    qqq[1] = 0.9875;
    qqq[2] = 0.988;
    qqq[3] = 0.9885;
    qqq[4] = 0.989;
    qqq[5] = 0.9895;
    qqq[6] = 0.99;
    qqq[7] = 0.9905;
    qqq[8] = 0.991;
    qqq[9] = 0.9915;
    qqq[10] = 0.9925;
    qqq[11] = 0.995;
    qqq[12] = 0.9955;
    qqq[13] = 0.996;
    qqq[14] = 0.9965;
    qqq[15] = 0.997;
    qqq[16] = 0.9975;
    qqq[17] = 0.998;
    qqq[18] = 0.9985;
    qqq[19] = 0.999;


long double q = qqq[(int)qq];//(long)qq*0.005 + 0.9;
long double distanz = pow(1/((1/q)-1), 1./3) * n*n;
long double omega = sqrt(1/((1-q)*distanz*distanz*distanz));

srand(time(NULL));
long double var_time = 0;
long double postime = 0;
long double endtime;

#include "particles.hpp"


    long double cloudranf[3];
    cloudranf[0] = 0;
    cloudranf[1] = 0;
    cloudranf[2] = 0;

    for (int f=2; f<anzahl; f++) {
        cloudranf[0] += particles[f].r[0]-particles[f%2].r[0];
        cloudranf[1] += particles[f].r[1]-particles[f%2].r[1];
        cloudranf[2] += particles[f].r[2]-particles[f%2].r[2];

    }
    cloudranf[0] /= (anzahl-2);
    cloudranf[1] /= (anzahl-2);
    cloudranf[2] /= (anzahl-2);


long double r_anf = abstand(cloudranf);

std::cout<<"q = "<<q<<endl;
std::cout<<"Distance between Atoms = "<<distanz<<endl;
std::cout<<"Angular Frequency = "<<omega<<endl;
std::cout<<"Initial Radius = "<<r_anf<<endl<<endl;

bool *thuruk_firstoutput;
std::string *thuruk_part;

thuruk_firstoutput=new bool[anzahl]; for(long int i=0; i<anzahl; ++i) thuruk_firstoutput[i]=true;
thuruk_part=new std::string[anzahl];

long double m_total = 0;
for (int i=0; i<anzahl; i++) {
    m_total += particles[i].mass;
}

com[0] = 0;
com[1] = 0;
com[2] = 0;

int zahl = 0;
bool dispersed=false;
long double disperse_time;

//mkdir("./outputfiles/");

//if (fileoutput) {
    stringstream ss1;
    ss1.str("");
    ss1<<"./outputfiles/cloud_average_"<<q<<".txt";
    ofstream file3(ss1.str().c_str());
    if (!file3.is_open()) {
        cout<<"File not opened.";
        return 0;
    }
    file3.close();
//}



stringstream pos_time;
pos_time<<"var timestep=[";

while(var_time<=time_max) {

//calculte new delta_t

if (!changingtime) {
    delta_t = delta_t_max;
} else {
    long double komischesumme1=0, komischesumme2=0;
    for (int i=0; i<anzahl; i++) {
        for (int j=0; j<anzahl; j++) {
            if (((particles[i].index == particles[j].index) || (particles[j].index == 0)) && (i!=j)) {
            if (abstand(particles[i].r, particles[j].r) != 0) {
                komischesumme1 += mypow(abstand(particles[i].p) / abstand(particles[i].r, particles[j].r) / particles[i].mass, 2);
                komischesumme2 += sqrt(mypow(k * particles[i].charge * particles[j].charge / mypow(abstand(particles[i].r, particles[j].r), 3)/particles[i].mass, 2));
            }
            }
        }
    }
    delta_t = min(delta_t_max, epsilon / sqrt(komischesumme1 + epsilon / 4. * komischesumme2));


}

var_time+=delta_t;
postime+=delta_t;


//calculate new positions

rungekutta(particles);

//execute changes

execute(particles);

//rotate frame

if (rotatingframe) {
    adjustrotmatrix(rotmatrix, -omega*delta_t);

    for (int i=0; i<anzahl; i++) {
        matrixmultvector(rotmatrix, particles[i].r);
        matrixmultvector(rotmatrix, particles[i].p);
    }
}

//output


if (shouldoutput) {
if (fileoutput) {

    //Electron Cloud

    if (cloudoutput) {
        for (int muh=0; muh<anzahl-4; muh++) {
            if (!filecloud[muh].is_open()) {
                std::cout<<"Filecloud["<<muh<<"] is not open!";
                return 0;
            } else {
                filecloud[muh]<<particles[muh+4].r[0]<<" "<<particles[muh+4].r[1]<<" "<<particles[muh+4].r[2]<<" ";
                filecloud[muh]<<particles[muh+4].p[0]<<" "<<particles[muh+4].p[1]<<" "<<particles[muh+4].p[2]<<std::endl;
            }
        }
    }

    //Average
    long double cloudr[3];
    long double cloudp[3];
    cloudr[0] = 0;
    cloudr[1] = 0;
    cloudr[2] = 0;
    cloudp[0] = 0;
    cloudp[1] = 0;
    cloudp[2] = 0;
    for (int f=2; f<anzahl; f++) {
        cloudr[0] += particles[f].r[0]-particles[f%2].r[0];
        cloudr[1] += particles[f].r[1]-particles[f%2].r[1];
        cloudr[2] += particles[f].r[2]-particles[f%2].r[2];
        cloudp[0] += particles[f].p[0]-particles[f%2].p[0];
        cloudp[1] += particles[f].p[1]-particles[f%2].p[1];
        cloudp[2] += particles[f].p[2]-particles[f%2].p[2];
    }
    cloudr[0] /= (anzahl-2);
    cloudr[1] /= (anzahl-2);
    cloudr[2] /= (anzahl-2);
    cloudp[0] /= (anzahl-2);
    cloudp[1] /= (anzahl-2);
    cloudp[2] /= (anzahl-2);

    ofstream file3(ss1.str().c_str(), ios::app);
    if (!file3.is_open()) {
        cout<<"File not opened.";
        return 0;
    }
    file3<<cloudr[0]<<" "<<cloudr[1]<<" "<<cloudr[2]<<" ";
    file3<<cloudp[0]<<" "<<cloudp[1]<<" "<<cloudp[2]<<" ";
    file3<<abstand(cloudr)<<" "<<var_time/2./PI/n/n/n<<std::endl;
    file3.close();


    //Positions
}


for (long int i=0; i<anzahl; ++i)
{
	if ((calcprotons) || (particles[i].charge==-1))
	{
		std::stringstream ss;
		if(thuruk_firstoutput[i])
		{
			thuruk_firstoutput[i]=false;
			ss << thuruk_part[i] << thuruk_output(particles[i]);
		}
		else
		{
			ss << thuruk_part[i] << ',' << thuruk_output(particles[i]);
		}
		thuruk_part[i]=ss.str();
	}
}


    std::cout/*<<"\r"<<var_time<<" / "<<time_max<<" - "*/<<var_time/2./PI/n/n/n<<" / "<<time_max/2./PI/n/n/n<<" Orbits for q="<<q<<"                    "<<endl;

    zahl=0;
    shouldoutput=false;

    pos_time<<postime<<",";
    endtime=postime;
    postime=0;

}

    long double cloudr[3];
    cloudr[0] = 0;
    cloudr[1] = 0;
    cloudr[2] = 0;
    for (int f=2; f<anzahl; f++) {
        cloudr[0] += particles[f].r[0]-particles[f%2].r[0];
        cloudr[1] += particles[f].r[1]-particles[f%2].r[1];
        cloudr[2] += particles[f].r[2]-particles[f%2].r[2];
    }
    cloudr[0] /= (anzahl-2);
    cloudr[1] /= (anzahl-2);
    cloudr[2] /= (anzahl-2);



if (!dispersed) {
    if (abstand(cloudr) < r_anf/2.) {
        disperse_time = var_time/2./PI/n/n/n;
        std::cout<<"Cancelled with radius "<<abstand(cloudr)<<endl;
        var_time+=time_max;
        pthread_mutex_lock(&mutex);
        ofstream file_q(ss_q.str().c_str(), ios::app);
        if (!file_q.is_open()) {
            cout<<"File_q not opened.";
            return 0;
        }
        file_q<<q<<" "<<disperse_time<<endl;
        file_q.close();
        pthread_mutex_unlock(&mutex);
        dispersed = true;
    }
}



zahl++;
if (zahl==framerate) {
    shouldoutput=true;
    //zahl=0;
    //std::cout<<"Ping "<<q<<"     "<<endl;
}
}

if (thurukoutput) {
stringstream ss4;
ss4.str("");
ss4<<"./outputfiles/positions_"<<q<<".js";
ofstream file4(ss4.str().c_str());
if(!file4.is_open()) {
    cout<<"File4 not opened.";
    return 0;
}
file4 << "var positions=[[";

	for(long int i=0; i<anzahl; ++i)
	{
		if(!calcprotons&&particles[i].charge>0)
		{
			if(i==anzahl-1) file4 << thuruk_output(particles[i]) << "]];";
			else file4 << thuruk_output(particles[i]) << "],[";
		}
		else
		{
			if(i==anzahl-1) file4 << thuruk_part[i] << "]];";
			else file4 << thuruk_part[i] << "],[";
		}
	}
	file4 << "var colors=[";
	for(long int i=0; i<anzahl; ++i)
	{
		if(!calcprotons&&particles[i].charge>0)
		{
			if(i==anzahl-1) file4 << "\"red\"];";
			else file4 << "\"red\",";
		}
		else
		{
			if((particles[i].charge<0) && !(particles[i].dispersed))
			{
				if(i==anzahl-1) file4 << "\"blue\"];";
				else file4 << "\"blue\",";
			}
			if ((particles[i].charge<0) && (particles[i].dispersed))
            {
                if(i==anzahl-1) file4 << "\"white\"];";
				else file4 << "\"white\",";
            }
			if (particles[i].charge==0)
			{
				if(i==anzahl-1) file4 << "\"white\"];";
				else file4 << "\"white\",";
			}
		}
	}
	file4 << "var pairs=[";
    std::stringstream *thuruk_pairs; thuruk_pairs = new std::stringstream[anzahl]; long int k=0;
    for(long int i=0; i<anzahl; ++i){
        for(long int j=0; j<anzahl; ++j) if(i<j&&particles[i].index==particles[j].index&&particles[i].charge<0&&particles[j].charge<0){

            thuruk_pairs[k] << '[' << i << ',' << j << ']';
            ++k;
        }
    }
    if(k>0) for(long int i=0; i<k; ++i){
        if(i!=(k-1)) file4 << thuruk_pairs[i].str() << ',';
        else file4 << thuruk_pairs[i].str();
    }
    file4 << "];";
    pos_time << endtime <<"];";
    file4<<pos_time.str();
	file4.close();
}

	if (dispersed) {
        std::cout<<"Dispersed at "<<disperse_time<<endl;
        //file_q<<qq<<" "<<disperse_time<<endl;
	} else {
        std::cout<<"Not Dispersed!"<<endl;
    }

pthread_exit(NULL);
}

int main (int argc, char *argv[])
{
   pthread_t threads[20];
   pthread_mutex_init(&mutex, NULL);
   int rc;
   long t;
   for(t=6; t<9; t++){
      printf("In main: creating thread %ld\n", t);
      pthread_create(&threads[t], NULL, haupt, (void *)t);
      if (rc){
         printf("ERROR; return code from pthread_create() is %d\n", rc);
         exit(-1);
      }
   }


   /* Last thing that main() should do */
   pthread_exit(NULL);
}


