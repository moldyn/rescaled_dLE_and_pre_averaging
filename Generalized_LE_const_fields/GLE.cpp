using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <list>
#include <ostream>
#include <math.h>
#include <vector>
#include <random>			//for random generator
#include <functional>		//for random generator as function
#include <Eigen/Cholesky>	//library used to do the Cholesky decomposition, library "Eigen" necessary to use the matrices/vectors for the Langevin propagation
#include <Eigen/LU>         //library used to calculate the inverse of matrices

int main(int argc, char* argv[]){
    
for (int i = 1; i < argc; ++i) {
  if(std::string(argv[i]) == "-h"){
    cout << " " <<endl;
    cout << "USAGE: ./GLE [OPTIONS]" << endl;
    cout << " " <<endl;
    cout << "INPUT: file containing negativ logarithm of the histogram representing the free energy, matrix of friction prefactors, masses of each coordinate, vector of memory decay times, temperature" << endl; 
    cout << "OUTPUT: trajectory generated according to the generalised Langevin equation, integration timestep (option -t) times output frequency (option -s) is the times difference between the points." <<endl;
    cout << " " <<endl;
    cout << "The program propagates the generalised Langevin equation based on given fields starting at a certain point (file at option -start)." << endl;
    cout << " " <<endl;
    cout << "The generalised Langevin equation is based on the free energy (drift), friction (memory included) and coloured noise. Based on F(x)=-kT*ln(P(x)), the file given to the" <<endl;
    cout << "program (option -free) needs to be (x1,x2,...,xn,-ln(P(x))) where P(x) represents the propability to be at x=(x1,x2,..,xn)^T. This means that we give some 'grid' to the" <<endl;
    cout << "program which is used to approximate the continouse free energy. The program assumes that each dimension has the same number of 'grid' points (specified by option -n), i.e., if the histogram" <<endl;
    cout << "was produced with the same number of bins per dimension everything is fine. Additionally, it is assumed that the grid points for each dimension are always seperated by the" <<endl;
    cout << "same distance. This is needed to circumvent tediouse searches for the actual free energy. To circumvent bin assignment problems, the free energy needs to be given in the following way" << endl;
    cout << " " <<endl;
    cout << "                                                                              x_0 y_0 z_0 -ln(P(x_0,y_0,z_0)) " <<endl;
    cout << "                                                                              x_0 y_0 z_1 -ln(P(x_0,y_0,z_1)) " <<endl;
    cout << "                                                                              x_0 y_0 z_2 -ln(P(x_0,y_0,z_2)) " <<endl;
    cout << "                                                                                           ........." <<endl;
    cout << "                                                                              x_0 y_1 z_0 -ln(P(x_0,y_1,z_0)) " <<endl;
    cout << "                                                                              x_0 y_1 z_1 -ln(P(x_0,y_1,z_1)) " <<endl;
    cout << "                                                                                           ........." <<endl;
    cout << "                                                                              x_1 y_0 z_0 -ln(P(x_1,y_0,z_0)) " <<endl;
    cout << "                                                                                           ........." <<endl;
    cout << "                                                                              x_n y_n z_n -ln(P(x_n,y_n,z_n)) " <<endl;
    cout << "with x_n=x_max and x_0=x_min (other coordinates similar)." <<endl;
    cout << " " <<endl;
    cout << "This means that the coordinate columns vary at first total right and at last total left. The derivative of the free energy is approximated by the difference between the local" <<endl;
    cout << "free energy and the free energy at the next grid point along the respective dimension. The differences in both directions (x_i > x_{i,actual} and x_i < x_{i,actual}) are averaged." <<endl;
    cout << "The memory kernel K(t) is a matrix and its elements are given by {K(t)}_{i,j}=gamma_{i,j}/tau_{i} * exp(-t/tau_{i}) which means that we assume an exponentially decaying memory kernel." <<endl;
    cout << "The gamma_{i,j} (no units) are given by the file at option -gamma, the tau_{i} (in ps) are given by the file at option -tau. Together with the temperature (option -T given in K) and" <<endl;
    cout << "the integration timestep (option -t) this is everything which we need. Based on normal distributed white noise vector xi (seed given by option -I), the coloured noise is propagated by" <<endl;
    cout << " " <<endl;
    cout << "                                    F_{i}(n*timestep)=exp(-timestep/tau_{i}) * F_{i}((n-1)*timestep) + sqrt(2kT)/tau_{i} * sqrt(timestep) * (sqrt(Gamma) * xi)_{i}" <<endl;
    cout << " " <<endl;
    cout << "where sqrt(Gamma) is the Cholesky decomposition of gamma, i.e., make sure that gamma is symmetric! Together with the iterative relation for the memory friction I(t)=integral_0^t(K(t-s)x_dot(s)ds)" << endl;
    cout << "which is given by" <<endl;
    cout << " " <<endl;
    cout << "                                   I_{i}(n*timestep)=exp(-timestep/tau_{i})I_{i}((n-1)*timestep) + 1/tau_{i}* exp(-timestep/tau_{i}) * timestep * (Gamma*x_dot((n-1)* timestep)_{i}" <<endl;
    cout << " " <<endl;
    cout << "and the masses along the different coordinates (option -mass) the trajectory is propagated by" <<endl;
    cout << " " <<endl;
    cout << "                                   x_dot(n*timestep)=x_dot((n-1)*timestep) + m^{-1}+(timestep * (-dF/dx((n-1*timestep)) - I((n-1)+timestep) + F((n-1)*timestep)))" <<endl;
    cout << " " <<endl;
    cout << "with m^{-1} being a diagonal matrix of the inverse masses. The position is propagated according to" <<endl;
    cout << " " <<endl;
    cout << "                                                   x(n*timestep)=x_((n-1)*timestep) + x_dot((n-1)*timestep) * timestep" <<endl;
    cout << " " <<endl;
    cout << " " <<endl;
    cout << "ATTENTION: The dimensionality (option -d) has to be given correctly since contradictions to the 'real' dimensionality (e.g., given by -gamma) are not detected." <<endl;
    cout << "ATTENTION: The free energy must have the same number of bins for all dimensions, contradictions are not detected, the format needs to fit as well (described above)" <<endl;
    cout << " " <<endl;
    cout << " " <<endl;
    cout << "OPTIONS" <<endl;
    cout << " " << " " << "-h show these lines" << endl;
    cout << " " << " " << "-start name file starting point [default startpoint]" << endl;
    cout << " " << " " << "-free name file histogram (free energy) as list (x1,...,xn,-ln(P(x1,..xn))) [default free_energy]" << endl;
    cout << " " << " " << "-gamma name file matrix prefactors of friction (Gamma11,Gamma12...) [default gamma]" << endl;
    cout << " " << " " << "-tau name file memory decay times (tau11,tau12...) [default tau]" << endl;
    cout << " " << " " << "-mass name file masses (m1,m2,...) [default mass]" << endl;
    cout << " " << " " << "-o name output trajectory (x1,..,xn) [default lang]" << endl;
    cout << " " << " " << "-t integration timestep [default 0.001 ps]" << endl;
    cout << " " << " " << "-T temperature [default 300 K]" << endl;
    cout << " " << " " << "-I seed random number generator [default: 0]" << endl;
    cout << " " << " " << "-L length of output trajectory [default: 100000 points]" << endl;
    cout << " " << " " << "-s write out every sth point [default: 1, i.e., every timestep]" << endl;
    cout << " " << " " << "-d dimensionality input [default: 1] ATTENTION: contradiction to real dimensionality of the input is not detected by program " << endl;
    cout << " " << " " << "-n number of 'grid' points per dimension of the free energy [default: 200]" << endl;
    cout << " " <<endl;
    cout << " " <<endl;
    exit(0);
  }
}

//Count Variables
int a,b,c,f;

//Input Variables
int length,dim,freeenergybins,everynth;
int seed;
double deltat;
double temperature;
std::string namefreeenergy="free_energy";
std::string namegamma="gamma";
std::string nametau="tau";
std::string nameoutput="lang";
std::string namestart="startpoint";
std::string namemass="mass";

length=100000;
dim=1;
seed=0;
deltat=0.001;
temperature=300;
freeenergybins=200;
everynth=1;

for(a=1;a<argc;a++){
    if(std::string(argv[a])=="-free"){
        if(a+1< argc){ 							// Make sure we aren't at the end of argv!
		namefreeenergy=argv[a+1]; 						// Increment 'i' so we don't get the argument as the next argv[i].
        }
    }
    else if(std::string(argv[a])=="-o"){
        if(a+1< argc){ 
		nameoutput=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-start"){
        if(a+1< argc){ 
		namestart=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-gamma"){
        if(a+1< argc){ 
		namegamma=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-tau"){
        if(a+1< argc){ 
		nametau=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-mass"){
        if(a+1< argc){ 
		namemass=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-d"){
        if(a+1<argc){
		dim=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-L"){
        if(a+1<argc){
		length=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-s"){
        if(a+1<argc){
		everynth=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-I"){
        if(a+1<argc){
		seed=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-t"){
        if(a+1<argc){
		deltat=atof(argv[a+1]);
        }
    }
    else if(std::string(argv[a])=="-T"){
        if(a+1<argc){
		temperature=atof(argv[a+1]);
        }
    }
    else if(std::string(argv[a])=="-n"){
        if(a+1<argc){
		freeenergybins=(int)(atof(argv[a+1]));
        }
    }
}

int numberfreeenergypoints=pow(freeenergybins,dim);

//Input fields

double transfer=0;

Eigen::VectorXd startpoint(dim);
Eigen::VectorXd maxcoordinates(dim);
Eigen::VectorXd mincoordinates(dim);
Eigen::VectorXd freeenergy(numberfreeenergypoints);
Eigen::VectorXd deltax(dim);
Eigen::VectorXd tauvector(dim);
Eigen::MatrixXd gammamatrix(dim,dim);
Eigen::MatrixXd choleskygamma(dim,dim);
Eigen::MatrixXd massesinvers(dim,dim);

ifstream datastart;                                                             //open input file
datastart.open(namestart,ios::in);
ifstream datagamma;                                                             //open input file
datagamma.open(namegamma,ios::in);
ifstream datatau;                                                             //open input file
datatau.open(nametau,ios::in);
ifstream datamass;                                                             //open input file
datamass.open(namemass,ios::in);

for(a=0;a<dim;a++){
    datatau >> transfer;
    tauvector(a)=transfer;
    datastart >> transfer;
    startpoint(a)=transfer;
    for(b=0;b<dim;b++){
        datagamma >> transfer;
        gammamatrix(a,b)=transfer;
        massesinvers(a,b)=0;
        if(b==a){
            datamass >> transfer;
            massesinvers(a,b)=1/transfer;
        }
    }
}

datastart.close();
datagamma.close();
datatau.close();
datamass.close();

choleskygamma=Eigen::LLT<Eigen::MatrixXd>(gammamatrix).matrixL();                       // needed for noise

ifstream datafreeenergy;                                                             //open input file
datafreeenergy.open(namefreeenergy,ios::in);

for(a=0;a<numberfreeenergypoints;a++){
    for(b=0;b<dim;b++){
        datafreeenergy >> transfer;
        if(a==0){
            mincoordinates(b)=transfer;
        }
        else if(a==numberfreeenergypoints-1){
            maxcoordinates(b)=transfer;
        }
    }
    datafreeenergy >> transfer;
    freeenergy(a)=transfer;
}

datafreeenergy.close();

for(a=0;a<dim;a++){
    if(maxcoordinates(a)<=mincoordinates(a)){
        cout << "Problems with the free energy format. Please check.";
        exit(0);
    }
}

for(a=0;a<dim;a++){
    mincoordinates(a)=mincoordinates(a)-0.000001;
    maxcoordinates(a)=maxcoordinates(a)+0.000001;
    deltax(a)=(maxcoordinates(a)-mincoordinates(a))/freeenergybins;
}

cout << "Starting point, Free energy, gamma, tau read in." << endl;
cout << "Starting point:" << endl;
for(a=0;a<dim;a++){
    cout << startpoint(a) << " ";
}
cout << endl;
cout << "Free energy range:" << endl;
cout << "Min:" << endl;
for(a=0;a<dim;a++){
    cout << mincoordinates(a) << " ";
}
cout << endl;
cout << "Max:" << endl;
for(a=0;a<dim;a++){
    cout << maxcoordinates(a) << " ";
}
cout << endl;
cout << "Delta x:" << endl;
for(a=0;a<dim;a++){
    cout << deltax(a) << " ";
}
cout << endl;
cout << "tau:" << endl;
for(a=0;a<dim;a++){
    cout << tauvector(a) << " ";
}
cout << endl;
cout << "Gamma:" << endl;
for(a=0;a<dim;a++){
    for(b=0;b<dim;b++){
    cout << gammamatrix(a,b) << " ";
    }
    cout << endl;
}
cout << endl;

//Matrices and vectors contributing to the gLE

double kT=38*temperature/300;
Eigen::VectorXd exponentialfactor(dim);
Eigen::MatrixXd gammaovertautimesdttimesexp(dim,dim);
Eigen::MatrixXd prefactorwhitenoise(dim,dim);

for(a=0;a<dim;a++){
    exponentialfactor(a)=exp((-1)*(deltat/tauvector(a)));
}
for(a=0;a<dim;a++){
    for(b=0;b<dim;b++){
        gammaovertautimesdttimesexp(a,b)=gammamatrix(a,b)/tauvector(a)*exponentialfactor(a)*deltat;
        prefactorwhitenoise(a,b)=sqrt(2*kT)*choleskygamma(a,b)/tauvector(a);
    }
}

//random number generator

if(seed<0){                                                                     //mt19337 works just with unsigned ints, I don't know exactly the influences of negative seeds
    seed=-seed;
}

std::function<float()> randomnumber;							                 //random generator as function
randomnumber=std::bind(std::normal_distribution<double>(0.0, 1.0)			//uniform distributions of real numbers between zero and one
                 , std::mt19937(seed));

cout << "Start propagation" << endl;

ofstream output;                                                    //open output file
output.open(nameoutput,ios::out | ios::trunc);

output << "#" << " ";
for(a=0;a<dim;a++){
    output << "x"<<a+1 << " ";
}
for(a=0;a<dim;a++){
    output << "drift"<<a+1 << " ";
}
for(a=0;a<dim;a++){
    output << "friction"<<a+1 << " ";
}
for(a=0;a<dim;a++){
    output << "Noise"<<a+1 << " ";
}
for(a=0;a<dim;a++){
    output << "Whitenoise"<<a+1 << " ";
}
output << endl;

double actualfreenergy=0,actualfreeenergyplusdeltax=0,actualfreenergyminusdeltax=0;
int helpdeterminebox=0;

Eigen::VectorXi actualpointbinnumbers(dim);
Eigen::VectorXd actualpoint(dim);
Eigen::VectorXd actualvelocity(dim);
Eigen::VectorXd lastvelocity(dim);
Eigen::VectorXd drift(dim);
Eigen::VectorXd friction(dim);
Eigen::VectorXd colourednoise(dim);
Eigen::VectorXd whitenoise(dim);

for(a=0;a<dim;a++){
    actualpoint(a)=startpoint(a);
    actualvelocity(a)=0;
    drift(a)=0;
    friction(a)=0;
    colourednoise(a)=0;
    whitenoise(a)=0;
}

c=0;
int alreadymentioned=0;

for(a=0;a<length;a++){
    for(f=0;f<everynth;f++){
        if(c%100000==0){
            if(alreadymentioned==0){
                cout << "gLE point" << " " << c << endl;
                alreadymentioned++;
            }
        }
        if(f==0){
            for(b=0;b<dim;b++){
                output << actualpoint(b)<< " ";
            }
            for(b=0;b<dim;b++){
                output << drift(b)<< " ";
            }
            for(b=0;b<dim;b++){
                output << friction(b)<< " ";
            }
            for(b=0;b<dim;b++){
                output << colourednoise(b)<< " ";
            }
            for(b=0;b<dim;b++){
                output << whitenoise(b)<< " ";
            }
            output << endl;
            c++;
            alreadymentioned=0;
        }
        for(b=0;b<dim;b++){                                                                                         //assign actual point to local free energy
            actualpointbinnumbers(b)=(int)((actualpoint(b)-mincoordinates(b))/(maxcoordinates(b)-mincoordinates(b))*freeenergybins);
            if(actualpointbinnumbers(b)<1){                                                                           //if the dLE point lies out of the value regions defined by the input, the last known free energy is used
                actualpointbinnumbers(b)=1;
            }
            else if(actualpointbinnumbers(b)>(freeenergybins-2)){
                actualpointbinnumbers(b)=freeenergybins-2;
            }
        }
        helpdeterminebox=0;
        for(b=0;b<dim;b++){                                                                                         //assign unique number to the actual point -> box in which the point is located
            helpdeterminebox=helpdeterminebox+(actualpointbinnumbers(b)*pow(freeenergybins,(dim-1-b)));
        }
        actualfreenergy=freeenergy(helpdeterminebox);
        for(b=0;b<dim;b++){                                                                                                                                             //calculate drift field
            helpdeterminebox=helpdeterminebox-(actualpointbinnumbers(b)*pow(freeenergybins,(dim-1-b)))+((actualpointbinnumbers(b)-1)*pow(freeenergybins,(dim-1-b)));
            actualfreenergyminusdeltax=freeenergy(helpdeterminebox);
            helpdeterminebox=helpdeterminebox-((actualpointbinnumbers(b)-1)*pow(freeenergybins,(dim-1-b)))+((actualpointbinnumbers(b)+1)*pow(freeenergybins,(dim-1-b)));
            actualfreeenergyplusdeltax=freeenergy(helpdeterminebox);
            helpdeterminebox=helpdeterminebox-((actualpointbinnumbers(b)+1)*pow(freeenergybins,(dim-1-b)))+(actualpointbinnumbers(b)*pow(freeenergybins,(dim-1-b)));
            drift(b)=-((actualfreeenergyplusdeltax-actualfreenergy)/deltax(b)+(actualfreenergy-actualfreenergyminusdeltax)/deltax(b))*0.5*kT;
        }
        for(b=0;b<dim;b++){
            whitenoise(b)=randomnumber();
        }
        actualpoint=actualpoint+actualvelocity*deltat;
        lastvelocity=actualvelocity;
        actualvelocity=actualvelocity+(massesinvers*((drift*deltat)-(friction*deltat)+(colourednoise*deltat)));
        for(b=0;b<dim;b++){
            friction(b)=friction(b)*exponentialfactor(b);
            colourednoise(b)=colourednoise(b)*exponentialfactor(b);
        }
        friction=friction+(gammaovertautimesdttimesexp*lastvelocity);
        colourednoise=colourednoise+((prefactorwhitenoise*whitenoise)*sqrt(deltat));
    }
}

output.close();
return EXIT_SUCCESS;
}
