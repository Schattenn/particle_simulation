
//initialize particles

particle particles[anzahl];

particles[0].name = "Core 1";
particles[0].r[0] = 0;//-distanz/2.;
particles[0].r[1] = 0;
particles[0].r[2] = -distanz/2.;//0;
particles[0].p[0] = 0;
particles[0].p[1] = 0;
particles[0].p[2] = 0;
particles[0].charge = 1;
particles[0].mass = 88. * 1836.1;
particles[0].index = 0;
particles[0].dispersed = false;

particles[1].name = "Core 2";
particles[1].r[0] = 0;//distanz/2.;
particles[1].r[1] = 0;
particles[1].r[2] = distanz/2.;//0;
particles[1].p[0] = 0;
particles[1].p[1] = 0;
particles[1].p[2] = 0;
particles[1].charge = 1;
particles[1].mass = 88. * 1836.1;
particles[1].index = 0;
particles[1].dispersed = false;

/*
circular states:

n=1, l=1:
r = 1, p = 1

n=300, l=300:
x = 90000, p = 1/300

in general:
r = n^2, p = 1/n



*/

//electron cloud

    stringstream ss_cloud;
    std::ofstream filecloud[anzahl-2];
    long double r[3], p[3];

//input from file
    std::ifstream ensemble("near_circular.dat");
    if (!ensemble.is_open()) {
        std::cout << "Couldn't open file!";
        return 0;
    }
    long double randcount, randanf = rand()%90000;
while ((ensemble >> r[0] >> r[1] >> r[2] >> p[0] >> p[1] >> p[2]) && randcount<randanf) {
    randcount++;
}

for (int dgang=2; (dgang<anzahl); ) {
if ((ensemble >> r[0] >> r[1] >> r[2] >> p[0] >> p[1] >> p[2]) /*&& (r[0] >= 0.75)*/) {
/*
    //gaussian distribution
        long double x1, x2, w, rx, ry;
        do {
            x1 = 2.0 * (long double)rand()/RAND_MAX - 1.0;
            x2 = 2.0 * (long double)rand()/RAND_MAX - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        rx = x1 * w;
        ry = x2 * w;

            //gaussian distribution
        long double rz, px;
        do {
            x1 = 2.0 * (long double)rand()/RAND_MAX - 1.0;
            x2 = 2.0 * (long double)rand()/RAND_MAX - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        rz = x1 * w;
        px = x2 * w;

            //gaussian distribution
        long double py, pz;
        do {
            x1 = 2.0 * (long double)rand()/RAND_MAX - 1.0;
            x2 = 2.0 * (long double)rand()/RAND_MAX - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        py = x1 * w;
        pz = x2 * w;
*/

        ss_cloud.str("");
        ss_cloud<<"Electron "<<dgang;

        particles[dgang].name = ss_cloud.str();
        particles[dgang].r[0] = r[0]/*rx * sigma/10.*/ + particles[dgang%2].r[0];
        particles[dgang].r[1] = r[1]/*ry * sigma + n*n*/ + particles[dgang%2].r[1];
        particles[dgang].r[2] = r[2]/*rz * sigma/10.*/ + particles[dgang%2].r[2];
        particles[dgang].p[0] = p[0]/*px * sigma/10.*/;
        particles[dgang].p[1] = p[1]/*py * sigma/10.*/;
        particles[dgang].p[2] = p[2]/*pz * sigma + 1/n*/;
        particles[dgang].charge = -1;
        particles[dgang].mass = 1;
        particles[dgang].index = (dgang - dgang%2)/2;
        particles[dgang].dispersed = false;
/*
        std::cout << r[0] << " / " << r[1] << " / " << r[2] << endl << p[0] << " / " << p[1] << " / " << p[2] << endl << endl;
        std::cout << particles[dgang].r[0] << " / " << particles[dgang].r[1] << " / " << particles[dgang].r[2] << endl;
        std::cout << particles[dgang].p[0] << " / " << particles[dgang].p[1] << " / " << particles[dgang].p[2] << endl << endl << endl;
*/
        if (anzahl>2) {
            ss_cloud.str("");
            ss_cloud<<"./outputfiles/cloud_"<<dgang<<".txt";
            if (cloudoutput) filecloud[dgang-2].open(ss_cloud.str().c_str());
        }
dgang++;
}
}


