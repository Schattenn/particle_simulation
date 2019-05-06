#include <cmath>
#include <cstdlib>

#include "parameters.hpp"

long double abstand(long double r[3]) {
    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

long double abstand(long double r1[3], long double r2[3]) {
    return sqrt((r2[0] - r1[0])*(r2[0] - r1[0]) + (r2[1] - r1[1])*(r2[1] - r1[1]) + (r2[2] - r1[2])*(r2[2] - r1[2]));
}

long double mypow(long double basis, int exponent) {
long double ret=1;
for (int i=0; i<exponent; i++) {
    ret *= basis;
}
return ret;
}

void matrixmultvector(long double matrix[3][3], long double vec[]) { //multiply vec and matrix, result is endvector
    long double erg=0.;
    long double endvector[3];
    for (int i=0; i<3; i++) {
        for (int k=0; k<3; k++) {
            erg += matrix[i][k] * vec[k];
        }
        endvector[i] = erg;
        erg = 0;
    }
    for (int i=0; i<3; i++) {
        vec[i] = endvector[i];
    }
}

void adjustrotmatrix(long double (&rotmatrix)[3][3], long double phi) {
    //x-axis
/*    rotmatrix[0][0] = 1.;
    rotmatrix[0][1] = 0.;
    rotmatrix[0][2] = 0.;
    rotmatrix[1][0] = 0.;
    rotmatrix[1][1] = cos(phi);
    rotmatrix[1][2] = -sin(phi);
    rotmatrix[2][0] = 0.;
    rotmatrix[2][1] = sin(phi);
    rotmatrix[2][2] = cos(phi);
*/

    //z-axis
    rotmatrix[0][0] = cos(phi);
    rotmatrix[0][1] = -sin(phi);
    rotmatrix[0][2] = 0.;
    rotmatrix[1][0] = sin(phi);
    rotmatrix[1][1] = cos(phi);
    rotmatrix[1][2] = 0.;
    rotmatrix[2][0] = 0.;
    rotmatrix[2][1] = 0.;
    rotmatrix[2][2] = 1.;

}

void execute(particle particles[]) {  //execute changes for part

for (int i=0; i<anzahl; i++) {
    if ((calcprotons) || (particles[i].charge==-1)) {
        particles[i].r[0] = particles[i].dummyr[0];
        particles[i].p[0] = particles[i].dummyp[0];
        particles[i].r[1] = particles[i].dummyr[1];
        particles[i].p[1] = particles[i].dummyp[1];
        particles[i].r[2] = particles[i].dummyr[2];
        particles[i].p[2] = particles[i].dummyp[2];
    }
    if (centerofmass) {
        long double m_total = 0;
        for (int m=0; m<anzahl; m++) {
            m_total += particles[m].mass;
        }
        for (int j=0; j<anzahl; j++) {
            com[0] += particles[j].r[0]*particles[j].mass;
            com[1] += particles[j].r[1]*particles[j].mass;
            com[2] += particles[j].r[2]*particles[j].mass;
        }
        com[0] /= m_total;
        com[1] /= m_total;
        com[2] /= m_total;
        particles[i].r[0] -= com[0];
        particles[i].r[1] -= com[1];
        particles[i].r[2] -= com[2];
    }
}
}

long double force(particle &part1, particle &part2, int index) { //returns the i-th coordinate of the force of part2 working on part1
if ((index<=2) && (abstand(part1.dummyr, part2.dummyr)!=0)) {
    return (k*part1.charge*part2.charge/mypow(abstand(part1.dummyr, part2.dummyr), 3)) * (part1.dummyr[index]-part2.dummyr[index]) ;
}
}

//kinetic / potential energy of part, or part1 and part2, respectively
long double ekin(particle &part) {
return mypow(abstand(part.p), 2) / 2. / part.mass;
}
long double epot(particle &part1, particle &part2) {
return k * part1.charge * part2.charge / abstand(part1.r, part2.r);
}



//Runge-Kutta
void rungekutta(particle particles[]) {

long double k1[anzahl][6], k2[anzahl][6], k3[anzahl][6], k4[anzahl][6];

for (int i=0; i<anzahl; i++) {
    for (int j=0; j<6; j++) {
        k1[i][j] = 0;
        k2[i][j] = 0;
        k3[i][j] = 0;
        k4[i][j] = 0;
    }
}

for (int i=0; i<anzahl; i++) {
        particles[i].dummyr[0] = particles[i].r[0];
        particles[i].dummyr[1] = particles[i].r[1];
        particles[i].dummyr[2] = particles[i].r[2];
        particles[i].dummyp[0] = particles[i].p[0];
        particles[i].dummyp[1] = particles[i].p[1];
        particles[i].dummyp[2] = particles[i].p[2];
}

for (int i=0; i<anzahl; i++) {
    if ((calcprotons) || (particles[i].charge==-1)) {
        k1[i][0] = particles[i].dummyp[0] / particles[i].mass;
        k1[i][1] = particles[i].dummyp[1] / particles[i].mass;
        k1[i][2] = particles[i].dummyp[2] / particles[i].mass;
        for (int j=0; j<anzahl; j++) {
            if (((particles[i].index == particles[j].index) || (particles[j].index == 0)) && (i!=j)) {
                if (abstand(particles[i].dummyr, particles[j].dummyr) == 0) {
                    std::cout << "Collision! ONOEZ!\n";
                } else {
                    k1[i][3] += force(particles[i], particles[j], 0);
                    k1[i][4] += force(particles[i], particles[j], 1);
                    k1[i][5] += force(particles[i], particles[j], 2);
                }
            }
        }
    }
}

for (int i=0; i<anzahl; i++) {
        particles[i].dummyr[0] = particles[i].r[0] + delta_t / 2. * k1[i][0];
        particles[i].dummyr[1] = particles[i].r[1] + delta_t / 2. * k1[i][1];
        particles[i].dummyr[2] = particles[i].r[2] + delta_t / 2. * k1[i][2];
        particles[i].dummyp[0] = particles[i].p[0] + delta_t / 2. * k1[i][3];
        particles[i].dummyp[1] = particles[i].p[1] + delta_t / 2. * k1[i][4];
        particles[i].dummyp[2] = particles[i].p[2] + delta_t / 2. * k1[i][5];
}

for (int i=0; i<anzahl; i++) {
    if ((calcprotons) || (particles[i].charge==-1)) {
        k2[i][0] = particles[i].dummyp[0] / particles[i].mass;
        k2[i][1] = particles[i].dummyp[1] / particles[i].mass;
        k2[i][2] = particles[i].dummyp[2] / particles[i].mass;
        for (int j=0; j<anzahl; j++) {
            if (((particles[i].index == particles[j].index) || (particles[j].index == 0)) && (i!=j)) {
                if (abstand(particles[i].dummyr, particles[j].dummyr) == 0) {
                    std::cout << "Collision! ONOEZ!\n";
                } else {
                    k2[i][3] += force(particles[i], particles[j], 0);
                    k2[i][4] += force(particles[i], particles[j], 1);
                    k2[i][5] += force(particles[i], particles[j], 2);
                }
            }
        }
    }
}


for (int i=0; i<anzahl; i++) {
        particles[i].dummyr[0] = particles[i].r[0] + delta_t / 2. * k2[i][0];
        particles[i].dummyr[1] = particles[i].r[1] + delta_t / 2. * k2[i][1];
        particles[i].dummyr[2] = particles[i].r[2] + delta_t / 2. * k2[i][2];
        particles[i].dummyp[0] = particles[i].p[0] + delta_t / 2. * k2[i][3];
        particles[i].dummyp[1] = particles[i].p[1] + delta_t / 2. * k2[i][4];
        particles[i].dummyp[2] = particles[i].p[2] + delta_t / 2. * k2[i][5];
}

for (int i=0; i<anzahl; i++) {
    if ((calcprotons) || (particles[i].charge==-1)) {
        k3[i][0] = particles[i].dummyp[0] / particles[i].mass;
        k3[i][1] = particles[i].dummyp[1] / particles[i].mass;
        k3[i][2] = particles[i].dummyp[2] / particles[i].mass;
        for (int j=0; j<anzahl; j++) {
            if (((particles[i].index == particles[j].index) || (particles[j].index == 0)) && (i!=j)) {
                if (abstand(particles[i].dummyr, particles[j].dummyr) == 0) {
                    std::cout << "Collision! ONOEZ!\n";
                } else {
                    k3[i][3] += force(particles[i], particles[j], 0);
                    k3[i][4] += force(particles[i], particles[j], 1);
                    k3[i][5] += force(particles[i], particles[j], 2);
                }
            }
        }
    }
}


for (int i=0; i<anzahl; i++) {
        particles[i].dummyr[0] = particles[i].r[0] + delta_t * k3[i][0];
        particles[i].dummyr[1] = particles[i].r[1] + delta_t * k3[i][1];
        particles[i].dummyr[2] = particles[i].r[2] + delta_t * k3[i][2];
        particles[i].dummyp[0] = particles[i].p[0] + delta_t * k3[i][3];
        particles[i].dummyp[1] = particles[i].p[1] + delta_t * k3[i][4];
        particles[i].dummyp[2] = particles[i].p[2] + delta_t * k3[i][5];
}

for (int i=0; i<anzahl; i++) {
    if ((calcprotons) || (particles[i].charge==-1)) {
        k4[i][0] = particles[i].dummyp[0] / particles[i].mass;
        k4[i][1] = particles[i].dummyp[1] / particles[i].mass;
        k4[i][2] = particles[i].dummyp[2] / particles[i].mass;
        for (int j=0; j<anzahl; j++) {
            if (((particles[i].index == particles[j].index) || (particles[j].index == 0)) && (i!=j)) {
                if (abstand(particles[i].dummyr, particles[j].r) == 0) {
                    std::cout << "Collision! ONOEZ!\n";
                } else {
                    k4[i][3] += force(particles[i], particles[j], 0);
                    k4[i][4] += force(particles[i], particles[j], 1);
                    k4[i][5] += force(particles[i], particles[j], 2);
                }
            }
        }
    }
}


for (int i=0; i<anzahl; i++) {
    if ((calcprotons) || (particles[i].charge==-1)) {
        particles[i].dummyr[0] = particles[i].r[0] + 1./6. * ( k1[i][0] + 2. * k2[i][0] + 2. * k3[i][0] + k4[i][0] ) * delta_t;
        particles[i].dummyr[1] = particles[i].r[1] + 1./6. * ( k1[i][1] + 2. * k2[i][1] + 2. * k3[i][1] + k4[i][1] ) * delta_t;
        particles[i].dummyr[2] = particles[i].r[2] + 1./6. * ( k1[i][2] + 2. * k2[i][2] + 2. * k3[i][2] + k4[i][2] ) * delta_t;
        particles[i].dummyp[0] = particles[i].p[0] + 1./6. * ( k1[i][3] + 2. * k2[i][3] + 2. * k3[i][3] + k4[i][3] ) * delta_t;
        particles[i].dummyp[1] = particles[i].p[1] + 1./6. * ( k1[i][4] + 2. * k2[i][4] + 2. * k3[i][4] + k4[i][4] ) * delta_t;
        particles[i].dummyp[2] = particles[i].p[2] + 1./6. * ( k1[i][5] + 2. * k2[i][5] + 2. * k3[i][5] + k4[i][5] ) * delta_t;
    }
}

}

