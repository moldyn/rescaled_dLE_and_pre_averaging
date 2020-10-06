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
    cout << "USAGE: ./dLE_rescaled_testmodel [OPTIONS]" << endl;
    cout << " " <<endl;
    cout << "INPUT: coordinate set x, column indicating endings of part trajecties by 0 and continouse parts by 1 (x1,...,xn,continuity)" << endl; 
    cout << "ADDITIONAL TESTMODEL INPUT: trajectory y, column indicating endings of part trajecties by 0 and continouse parts by 1 (y1,...,yn,continuity)." << endl; 
    cout << "ADDITIONAL INPUT:  matrix to rescale friction and noise to slow down/speed up system dynamics" << endl;
    cout << "OUTPUT: trajectory as given as ADDITIONAL TESTMODEL INPUT, field estimates, noise" <<endl;
    cout << " " <<endl;
    cout << "The program constructs a testmodel according to the data-driven Langevin equation (second order, i.e., Markovian) based on a given input file (option -i)." << endl;
    cout << "The input file contains the trajectory and a column which indicates the end of part trajectories by 0 (and continuity by 1)." << endl;
    cout << "The trajectory which should be tested is given as well (option -itest) so that the basis for the field estimates do not need to be the same trajectory." <<endl;
    cout << " " <<endl;
    cout << "For a certain number of ADDITIONAL INPUT frames (option -L, default full file), the dLE looks for the next neighbours (given in the input file -i) " <<endl;
    cout << "and takes a certain number of them (option -k) to estimate the fields. Having estimated the fields, the matrix M (given as file via -alterfric)" <<endl;
    cout << "is used to change friction and noise so that the Fluctuation-Dissipation theorem stays intact, i.e., K'K'^T=sqrt(M)*KK^T*sqrt(M)^T=2kT(sqrt(M)*Gamma*sqrt(M)^T)" <<endl;
    cout << "is calculated. M is meant to be a diagonal matrix where the diagonal entries quantify the 'increase of viskosity' for this coordinate." <<endl;
    cout << " " <<endl;
    cout << "The input -itest needs to have an additional column which indicates ends of trajectory pieces by 0. So, 1 marks concatenated pieces of the trajectory." <<endl;
    cout << " " <<endl;
    cout << "ATTENTION: closeness of points is defined by the Eucledian distance between them, so all input coordinates need to have the same 'scale'" <<endl;
    cout << "for distances to prevent that closeness along one dimension is favored to closeness along another one (e.g., x1 out of [0,1], x2 out of [-180,180]" <<endl;
    cout << "would lead to neighbourhoods which massively favour closeness along x1). For the sake of transparence, the dLE does not apply any rescaling." <<endl;
    cout << "If it is wished, the write out of the field estimates could be surpressed (option -nofields set to everything but 0)." <<endl;
    cout << " " <<endl;
    cout << "To speed up the neighbourhood determination, the points used for the field estimates (-i) are assigned to boxes so that not all points have to be" <<endl;
    cout << "scanned in each dLE step (box assisted search). The total number of boxes can be specified by option -b (default 10000). By default, the full coordinate" <<endl;
    cout << "space is divided into these boxes but it is possible to do the boxing only in the first two dimensions of the coordinate space (option -subspace). This leads to" <<endl;
    cout << "errors in the field estimation but has the advantage that the speed up of the neighbourhood determination due to the box assistance could be preserved. For high-" <<endl;
    cout << "dimensional data this speed up could fade away if the boxing is done in the full dimensionality." << endl;
    cout << "ATTENTION: The two-dimensional subspace is just used to assign the points to boxes. To determine the distance between actual dLE point and potential neighbour," <<endl;
    cout << "the full dimensional Eucledian distance is used." <<endl;
    cout << " " <<endl;
    cout << "ATTENTION: The dimensionality (option -d) has to be given correctly since contradictions to the 'real' dimensionality (input file) are not detected." <<endl;
    cout << "ATTENTION: the input files need to be complete and the columns have to be ordered in the right way, missing columns or wrong orders are not detected." <<endl;
    cout << " " <<endl;
    cout << " " <<endl;
    cout << "OPTIONS" <<endl;
    cout << " " << " " << "-h show these lines" << endl;
    cout << " " << " " << "-i name input trajectory (x1,...,xn,continuity) [default input]" << endl;
    cout << " " << " " << "-itest name input trajectory (y1,...,yn,future) on which the testmodel should be applied [default input.test]" << endl;
    cout << " " << " " << "-o name output trajectory (x1,..,xn,f1,...,fn,g11,g12,...,K11,K21,...,xi1,...,xin,max_distance) [default input.lang]" << endl;
    cout << " " << " " << "-k number of considered next neighbours [default 100]" << endl;
    cout << " " << " " << "-L length of tested trajectory piece [default: full file -itest]" << endl;
    cout << " " << " " << "-b number of boxes for box assisted search [default: 10000]" << endl;
    cout << " " << " " << "-nofields fields should not be written out is indicated by everything but 0 [default: 0]" << endl;
    cout << " " << " " << "-subspace box assignment of input -i should only be done in the first two dimensions is indicated by everything but 0 [default: 0]" << endl;
    cout << " " << " " << "-d dimensionality input [default: 1] ATTENTION: contradiction to real dimensionality of the input is not detected by program " << endl;
    cout << " " << " " << "-alterfric name matrix file used to alter friction and noise [default: modifyfriction] " << endl;
    cout << " " <<endl;
    cout << " " <<endl;
    exit(0);
  }
}

//Static parameters
int boxesperdim;                                                        //Number of boxes per dimension, are calculated later on as pow(numberboxes,(1/dim)).
int consideredpast=2;

//Input Variables
int nrneighbours,length,fieldsno,dim,subspace;
std::string nameinput="input";
std::string nametest="input.test";
std::string nameoutput="input.lang";
std::string namealterfric="modifyfriction";

int numberboxes=10000;                                                  //This is the number of boxes to which input points are assigend later on before the dLE propagation (i.e., the neighbourhood serchings) starts.
nrneighbours=100;
length=-1;
fieldsno=0;
dim=1;
subspace=0;

//Count/Helping variables
int a,b,c,d;
std::vector<int> considerpoint;                                                         //use vector to check whether input point can be used for field estimate
std::vector<double> transferposition;                                                   //vectors to store positions, velocities when read from input or later on from storing vectors
std::vector<double> actualposition;
std::vector<double> lastposition;
std::vector<double> transferforwarddisplacement;
std::vector<double> transferbackwarddisplacement;
int transferstatus=0,checkstatus=0;
int lengthinput=0;
std::vector<double> mininput;                                                           //Here, I store the value range which is spanned by the input data
std::vector<double> maxinput;
std::vector<int> boxindexactualpoint;                                                   //Store index of actual point box here
std::vector<int> boxindexprevpoint;                                                     //store index of previouse dLE point here to skip useless determination of neighbouring boxes if box was not left
std::vector<int> boxindexboxcomparetopoint;                                             //Store index of box which should be tested whether it is a neighbour of the actual point box
int helpdeterminebox=0;
int summedweightneighbourboxes=0;                                                       //use this variable to determine how many "boxlayers" around the box of the actual dLE point are nedded to find enough neighbourweigth
int numberlayers=0;                                                                     //use this variable to store the number of "boxlayers" which are needed for the acctual dLE point
int bindistance=0;
int checkboxhaschanged=0;
int sumweight=0;
double maxdistance=0;
double transferalterfric=0;
int statusprevpoint=0,statusactualpoint=0,statusnextpoint=0;
int lengthtesttrajectory=0;
int waswrittenout=0;

int checkedpoints=0, sortedpoints=0;

//variables for box assisted neighbourhood determination
std::vector<std::vector<std::vector<double>>> positionsassigntoboxes;                               //"box assisted" means that I assign all input points to boxes which form the full coordinate space, so that the nearest neighbour search only needs to scan the neighbouring boxes of the acctual bin
std::vector<std::vector<std::vector<double>>> forwardvelocitiesassignedtoboxes;
std::vector<std::vector<std::vector<double>>> backwardvelocityassigendtoboxes;
std::vector<int> entriesassigendtoboxes;
std::vector<int> boxshouldbeconsidered;

//variables to construct/store the neighbourhoodlist
std::vector<std::vector<double>> neighbourhoodpositions;
std::vector<std::vector<double>> neighbourhoodforwarddisplacements;
std::vector<std::vector<double>> neighbourhoodbackwarddisplacements;
std::vector<double> neighbourhooddistances;
std::vector<double> replaceposition;
std::vector<double> replaceforwarddisplacements;
std::vector<double> replacebackwarddisplacements;
double transferdistance=0,replacedistance=0;

std::vector<std::vector<double>> neighbourhoodordering;
int replacebox=0,replacepointnumber=0,transferbox=0,transferpointnumber=0;

for(a=1;a<argc;a++){
    if(std::string(argv[a])=="-i"){
        if(a+1< argc){ 							// Make sure we aren't at the end of argv!
		nameinput=argv[a+1]; 						// Increment 'i' so we don't get the argument as the next argv[i].
        }
    }
    else if(std::string(argv[a])=="-itest"){
        if(a+1< argc){
		nametest=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-o"){
        if(a+1< argc){ 
		nameoutput=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-d"){
        if(a+1<argc){
		dim=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-k"){
        if(a+1<argc){
		nrneighbours=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-b"){
        if(a+1<argc){
		numberboxes=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-L"){
        if(a+1<argc){
		length=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-nofields"){
        if(a+1<argc){
		fieldsno=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-subspace"){
        if(a+1<argc){
		subspace=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-alterfric"){
        if(a+1< argc){ 
		namealterfric=argv[a+1];
        }
    }
}

for(a=0;a<numberboxes;a++){
    positionsassigntoboxes.push_back(std::vector<std::vector<double>>());
    forwardvelocitiesassignedtoboxes.push_back(std::vector<std::vector<double>>());
    backwardvelocityassigendtoboxes.push_back(std::vector<std::vector<double>>());
}

for(a=0;a<numberboxes;a++){
    entriesassigendtoboxes.push_back(0);
    boxshouldbeconsidered.push_back(0);
}

for(a=0;a<consideredpast;a++){
    considerpoint.push_back(0);
}

cout << "Options detected" << endl;

//dLE vectors/fields and propagation
Eigen::VectorXd position(dim);
Eigen::VectorXd prev_position(dim);
Eigen::VectorXd next_position(dim);
Eigen::VectorXd forwardvelocity(dim);
Eigen::VectorXd backwardvelocity(dim);
Eigen::MatrixXd covarianceforwardforward(dim,dim);
Eigen::MatrixXd covariancebackwardbackward(dim,dim);
Eigen::MatrixXd covariancebackwardbackwardinverse(dim,dim);
Eigen::MatrixXd covarianceforwardbackward(dim,dim);
Eigen::VectorXd drift(dim);
Eigen::MatrixXd friction(dim,dim);
Eigen::MatrixXd frictiontranspose(dim,dim);
Eigen::MatrixXd diffusiontransposediffusion(dim,dim);
Eigen::MatrixXd diffusion(dim,dim);
Eigen::MatrixXd diffusiontranspose(dim,dim);
Eigen::MatrixXd diffusioninverse(dim,dim);
Eigen::MatrixXd diffusiontransposeinverse(dim,dim);
Eigen::MatrixXd helpmatrix(dim,dim);
Eigen::MatrixXd mass(dim,dim);
Eigen::MatrixXd masstranspose(dim,dim);
Eigen::MatrixXd massinverse(dim,dim);
Eigen::VectorXd noise(dim);
Eigen::MatrixXd alterfrictionmatrix(dim,dim);
Eigen::MatrixXd alterfrictionmatrixtranspose(dim,dim);

for(a=0;a<dim;a++){
  transferposition.push_back(0);                                                        //initialise with zeros
  lastposition.push_back(0);
  actualposition.push_back(0);
  transferbackwarddisplacement.push_back(0);
  transferforwarddisplacement.push_back(0);
  replaceposition.push_back(0);
  replaceforwarddisplacements.push_back(0);
  replacebackwarddisplacements.push_back(0);
  boxindexactualpoint.push_back(0);
  boxindexprevpoint.push_back(-1);
  boxindexboxcomparetopoint.push_back(0);
}

for(a=0;a<nrneighbours;a++){
    neighbourhoodpositions.push_back(std::vector<double>());
    neighbourhoodforwarddisplacements.push_back(std::vector<double>());
    neighbourhoodbackwarddisplacements.push_back(std::vector<double>());
    for(b=0;b<dim;b++){
        neighbourhoodpositions[a].push_back(0);
        neighbourhoodforwarddisplacements[a].push_back(0);
        neighbourhoodbackwarddisplacements[a].push_back(0);
    }
    neighbourhooddistances.push_back(-1);
}

for(a=0;a<nrneighbours;a++){
    neighbourhoodordering.push_back(std::vector<double>());
    for(b=0;b<3;b++){
        neighbourhoodordering[a].push_back(0);
    }
}

cout << "Vectors initialised" << endl;

if(dim<3){                                                                  //if dim<3, the restriction to a two-dimensional subspace makes no sense
    subspace=0;
}

if(subspace==0){
    boxesperdim=(int)(pow(numberboxes,(1.0/dim)));                                            //calculate the number of boxes needed later on if the box assistance should be applied to the full dimensional space
    if(boxesperdim<1){
        boxesperdim=1;
    }
}
else{
    boxesperdim=(int)(sqrt(numberboxes));                                              //calculate the number of boxes needed later on if the box assiatance should be applied only on a two-dimensional subspace
}

cout << "Boxes per dimension" << " " << boxesperdim <<endl;

b=0;

ifstream datafieldestimate;                                                             //open input file
datafieldestimate.open(nameinput,ios::in);

while(datafieldestimate>>transferposition[0]){                                          //read in input, determine value range and dLE startposition
    for(a=1;a<dim;a++){
        datafieldestimate>>transferposition[a];
    }
    datafieldestimate>>transferstatus;
    lengthinput++;
    if(b==0){                                                                              //determine value range of input and starting point of dLE (last input point)
        for(a=0;a<dim;a++){
            mininput.push_back(transferposition[a]);
            maxinput.push_back(transferposition[a]);
        }
    }
    else{
        for(a=0;a<dim;a++){
            if(transferposition[a]<mininput[a]){
                mininput[a]=transferposition[a];
            }
            else if(transferposition[a]>maxinput[a]){
                maxinput[a]=transferposition[a];
            }
        }
    }
    b++;
}


for(a=0;a<dim;a++){                                                                         //prevent possible problems at the borders of the value ranges
    mininput[a]=mininput[a]-0.0001;
    maxinput[a]=maxinput[a]+0.0001;
}

cout << "Value ranges Input determined." <<endl;

datafieldestimate.clear();                                                                  //Return to beginning of input file
datafieldestimate.seekg(0, ios::beg);

for(a=0;a<lengthinput;a++){                                                                 //assign input to the various boxes
    for(b=0;b<dim;b++){
        lastposition[b]=actualposition[b];
        actualposition[b]=transferposition[b];
        datafieldestimate>>transferposition[b];
    }
    datafieldestimate>>transferstatus;
    checkstatus=0;
    for(b=0;b<consideredpast;b++){
        if(considerpoint[b]==0){
            checkstatus++;
        }
    }
    if(checkstatus==0){
        if(subspace==0){
            for(b=0;b<dim;b++){                                                                                         //determine boxindex of actual point if the box assistance is used in the full dimensional space
                boxindexactualpoint[b]=(int)((maxinput[b]-actualposition[b])/(maxinput[b]-mininput[b])*boxesperdim);
            }
            helpdeterminebox=0;
            for(b=0;b<dim;b++){                                                                                         //assign unique number to the actual point -> box in which the point is located
                helpdeterminebox=helpdeterminebox+(boxindexactualpoint[b]*pow(boxesperdim,b));
            }
        }
        else{
            for(b=0;b<2;b++){                                                                                         //determine boxindex of actual point if only a two-dimensional subspace is used
                boxindexactualpoint[b]=(int)((maxinput[b]-actualposition[b])/(maxinput[b]-mininput[b])*boxesperdim);
            }
            helpdeterminebox=0;
            for(b=0;b<2;b++){                                                                                         //assign unique number to the actual point -> box in which the point is located
                helpdeterminebox=helpdeterminebox+(boxindexactualpoint[b]*pow(boxesperdim,b));
            }
        }
        positionsassigntoboxes[helpdeterminebox].push_back(actualposition);
        for(b=0;b<dim;b++){
            transferforwarddisplacement[b]=transferposition[b]-actualposition[b];
            transferbackwarddisplacement[b]=actualposition[b]-lastposition[b];
        }
        forwardvelocitiesassignedtoboxes[helpdeterminebox].push_back(transferforwarddisplacement);
        backwardvelocityassigendtoboxes[helpdeterminebox].push_back(transferbackwarddisplacement);
        entriesassigendtoboxes[helpdeterminebox]++;
    }
    for(b=0;b<consideredpast-1;b++){
        considerpoint[b]=considerpoint[b+1];
    }
    considerpoint[consideredpast-1]=transferstatus;
}


cout << "Box assignment performed" << endl;

datafieldestimate.close();

ifstream dataalterfriction;                                                             //open input file
dataalterfriction.open(namealterfric,ios::in);

cout << "Alter friction by following matrix" << endl;

for(a=0;a<dim;a++){
    for(b=0;b<dim;b++){
        dataalterfriction >> transferalterfric;
        alterfrictionmatrix(a,b)=sqrt(transferalterfric);
        cout << alterfrictionmatrix(a,b)<< " ";
    }
    cout << endl;
}

dataalterfriction.close();
alterfrictionmatrixtranspose=alterfrictionmatrix.transpose();

ifstream datatest;                                                             //open input file
datatest.open(nametest,ios::in);

while(datatest>>transferposition[0]){                                          //determine length of testmodel trajectory
    if(dim>1){
        for(a=1;a<dim;a++){
            datatest>>transferposition[a];
        }
    }
    datatest>>statusactualpoint;
    lengthtesttrajectory++;
}

if(length==-1){
    length=lengthtesttrajectory-2;                                              // -2 because first and last point of trajectory are neglected in every case.
}
else if(length>(lengthtesttrajectory-2)){                                       // make sure length parameter is not larger than possible considered points.
    length=lengthtesttrajectory-2;
}

datatest.clear();                                                                  //Return to beginning of input file
datatest.seekg(0, ios::beg);

ofstream output;                                                    //open output file
output.open(nameoutput,ios::out | ios::trunc);

output << "#" << " ";                                               //insert comment line
for(c=0;c<dim;c++){
    output << "x"<<c+1 << " ";
}
if(fieldsno==0){
    for(c=0;c<dim;c++){
        output << "f"<<c+1 << " ";
    }
    for(c=0;c<dim;c++){
        for(d=0;d<dim;d++){
            output << "g_("<<c+1<<","<<d+1<<")" << " ";
        }
    }
    for(c=0;c<dim;c++){
        for(d=0;d<dim;d++){
            if(d<=c){
                output << "K_("<<c+1<<","<<d+1<<")" << " ";
            }
        }
    }
}
for(c=0;c<dim;c++){
    output << "xi"<<c+1 << " ";
}
output << "future";
if(fieldsno==0){
    output << " " << "distance";
}
output<<endl;

cout << "Start testmodel calculation" << endl;

for(a=0;a<length;a++){
    if(a%100000==0){
        cout << "dLE point" << " " << a << endl;
    }
    if(a==0){
        for(b=0;b<dim;b++){                                                                                     //initialise points of testmodel trajectory
            datatest>>prev_position(b);
        }
        datatest>>statusprevpoint;
        for(b=0;b<dim;b++){
            datatest>>position(b);
        }
        datatest>>statusactualpoint;
        for(b=0;b<dim;b++){
            datatest>>next_position(b);
        }
        datatest>>statusnextpoint;
    }
    else{                                                                                                           //Later on: goo through trajectory
        if(waswrittenout==1){
            for(b=0;b<dim;b++){
                boxindexprevpoint[b]=boxindexactualpoint[b];
            }
        }
        prev_position=position;
        position=next_position;
        statusprevpoint=statusactualpoint;
        statusactualpoint=statusnextpoint;
        for(b=0;b<dim;b++){
            datatest>>next_position(b);
        }
        datatest>>statusnextpoint;
    }
    waswrittenout=0;
    if(statusactualpoint!=0){                                                                                               //only consider points which have a follower
        if(statusprevpoint!=0){  
            if(subspace==0){
                for(b=0;b<dim;b++){                                                                                         //assign actual point to one of the boxes, box assistance on full dimensional space
                    boxindexactualpoint[b]=(int)((maxinput[b]-position(b))/(maxinput[b]-mininput[b])*boxesperdim);
                    if(boxindexactualpoint[b]<0){                                                                           //if the dLE point lies out of the value regions defined by the input, it is assigned to the nearest box
                        boxindexactualpoint[b]=0;
                    }
                    else if(boxindexactualpoint[b]>(boxesperdim-1)){
                        boxindexactualpoint[b]=boxesperdim-1;
                    }
                }
                checkboxhaschanged=0;
                for(b=0;b<dim;b++){                                                                                            //check if the dLE left the box in which it was located in the prevoiuse step
                    if(boxindexactualpoint[b]!=boxindexprevpoint[b]){
                        checkboxhaschanged++;
                    }
                }
            }
            else{
                for(b=0;b<2;b++){                                                                                         //assign actual point to one of the boxes, box assistance on subspace
                    boxindexactualpoint[b]=(int)((maxinput[b]-position(b))/(maxinput[b]-mininput[b])*boxesperdim);
                    if(boxindexactualpoint[b]<0){                                                                           //if the dLE point lies out of the value regions defined by the input, it is assigned to the nearest box
                        boxindexactualpoint[b]=0;
                            }
                    else if(boxindexactualpoint[b]>(boxesperdim-1)){
                        boxindexactualpoint[b]=boxesperdim-1;
                    }
                }
                checkboxhaschanged=0;
                for(b=0;b<2;b++){                                                                                            //check if the dLE left the box in which it was located in the prevoiuse step
                    if(boxindexactualpoint[b]!=boxindexprevpoint[b]){
                        checkboxhaschanged++;
                    }
                }
            }
            //cout << "After index determination" << endl;
            if(checkboxhaschanged!=0){
                summedweightneighbourboxes=0;
                numberlayers=0;
                for(c=0;c<numberboxes;c++){
                    boxshouldbeconsidered[c]=0;
                    helpdeterminebox=c;
                    if(subspace==0){
                        boxindexboxcomparetopoint[0]=helpdeterminebox%boxesperdim;                                          //determine indices of actual point, box assitance on full dimensionality
                        for(d=1;d<dim;d++){
                            helpdeterminebox=helpdeterminebox/boxesperdim;
                            boxindexboxcomparetopoint[d]=helpdeterminebox%boxesperdim;
                        }
                        for(d=0;d<dim;d++){                                                                                 //check if box is neighbour of box which contains the actual point
                            bindistance=boxindexactualpoint[d]-boxindexboxcomparetopoint[d];
                            if(bindistance<0){
                                bindistance=-bindistance;
                            }
                            if(bindistance>boxshouldbeconsidered[c]){
                                boxshouldbeconsidered[c]=bindistance;
                            }
                        }
                    }
                    else{
                        boxindexboxcomparetopoint[0]=helpdeterminebox%boxesperdim;                                          //determine indices of actual point, box assitance in two-dimensional subspace
                        for(d=1;d<2;d++){
                            helpdeterminebox=helpdeterminebox/boxesperdim;
                            boxindexboxcomparetopoint[d]=helpdeterminebox%boxesperdim;
                        }
                        for(d=0;d<2;d++){                                                                                 //check if box is neighbour of box which contains the actual point
                            bindistance=boxindexactualpoint[d]-boxindexboxcomparetopoint[d];
                            if(bindistance<0){
                                bindistance=-bindistance;
                            }
                            if(bindistance>boxshouldbeconsidered[c]){
                                boxshouldbeconsidered[c]=bindistance;
                            }
                        }
                    }
                }
                while(summedweightneighbourboxes<nrneighbours){ 
                    numberlayers++;
                    summedweightneighbourboxes=0;
                    for(d=0;d<numberboxes;d++){
                        if(boxshouldbeconsidered[d]<=numberlayers){
                            summedweightneighbourboxes=summedweightneighbourboxes+entriesassigendtoboxes[d];
                            boxshouldbeconsidered[d]=1;                                                                     //set to one for consistancy with script parts below
                        }
                    }
                }
            }
            for(b=0;b<nrneighbours;b++){
                for(c=0;c<dim;c++){
                    neighbourhoodpositions[b][c]=0;
                    neighbourhoodforwarddisplacements[b][c]=0;
                    neighbourhoodbackwarddisplacements[b][c]=0;
                }
                neighbourhooddistances[b]=-1;
            }
            for(b=0;b<nrneighbours;b++){
                neighbourhoodordering[b][0]=-1;
                neighbourhoodordering[b][1]=0;
                neighbourhoodordering[b][2]=0;
            }
            for(b=0;b<numberboxes;b++){                                                                                     //determine actual k nearest neighbours
                if(boxshouldbeconsidered[b]==1){                                                                            //use box assistance, i.e., the result of the binning above
                    for(c=0;c<entriesassigendtoboxes[b];c++){
                        transferposition=positionsassigntoboxes[b][c];
                        transferdistance=0;
                        transferbox=b;
                        transferpointnumber=c;
                        for(d=0;d<dim;d++){                                                                                 //determine the distance between dLE point and actual input point
                            transferdistance=transferdistance+pow((position(d)-transferposition[d]),2);
                        }
                        checkedpoints++;
                        if(neighbourhoodordering[nrneighbours-1][0]==-1){
                            sortedpoints++;
                            for(d=0;d<nrneighbours;d++){                                                                          //determine whether the actual point is part of the numweight next neighbours
                                if(neighbourhoodordering[d][0]>transferdistance){
                                    replacebox=neighbourhoodordering[d][1];
                                    replacepointnumber=neighbourhoodordering[d][2];
                                    replacedistance=neighbourhoodordering[d][0];
                                    neighbourhoodordering[d][0]=transferdistance;
                                    neighbourhoodordering[d][1]=transferbox;
                                    neighbourhoodordering[d][2]=transferpointnumber;
                                    transferdistance=replacedistance;
                                    transferbox=replacebox;
                                    transferpointnumber=replacepointnumber;
                                }
                                else if(neighbourhoodordering[d][0]==-1){
                                    replacebox=neighbourhoodordering[d][1];
                                    replacepointnumber=neighbourhoodordering[d][2];
                                    replacedistance=neighbourhoodordering[d][0];
                                    neighbourhoodordering[d][0]=transferdistance;
                                    neighbourhoodordering[d][1]=transferbox;
                                    neighbourhoodordering[d][2]=transferpointnumber;
                                    transferdistance=replacedistance;
                                    transferbox=replacebox;
                                    transferpointnumber=replacepointnumber;
                                    d=nrneighbours;
                                }
                            }
                        }
                        else if(transferdistance<neighbourhoodordering[nrneighbours-1][0]){
                            sortedpoints++;
                            for(d=0;d<nrneighbours;d++){                                                                          //determine whether the actual point is part of the numweight next neighbours
                                if(neighbourhoodordering[d][0]>transferdistance){
                                    replacebox=neighbourhoodordering[d][1];
                                    replacepointnumber=neighbourhoodordering[d][2];
                                    replacedistance=neighbourhoodordering[d][0];
                                    neighbourhoodordering[d][0]=transferdistance;
                                    neighbourhoodordering[d][1]=transferbox;
                                    neighbourhoodordering[d][2]=transferpointnumber;
                                    transferdistance=replacedistance;
                                    transferbox=replacebox;
                                    transferpointnumber=replacepointnumber;
                                }
                                else if(neighbourhoodordering[d][0]==-1){
                                    replacebox=neighbourhoodordering[d][1];
                                    replacepointnumber=neighbourhoodordering[d][2];
                                    replacedistance=neighbourhoodordering[d][0];
                                    neighbourhoodordering[d][0]=transferdistance;
                                    neighbourhoodordering[d][1]=transferbox;
                                    neighbourhoodordering[d][2]=transferpointnumber;
                                    transferdistance=replacedistance;
                                    transferbox=replacebox;
                                    transferpointnumber=replacepointnumber;
                                    d=nrneighbours;
                                }
                            }
                        }
                    }
                }
            }
            //cout << "Before collecting neighbourhoodlist" << endl;
            for(b=0;b<nrneighbours;b++){
                if(neighbourhoodordering[b][0]!=-1){
                transferbox=neighbourhoodordering[b][1];
                transferpointnumber=neighbourhoodordering[b][2];
                neighbourhoodpositions[b]=positionsassigntoboxes[transferbox][transferpointnumber];
                neighbourhoodbackwarddisplacements[b]=backwardvelocityassigendtoboxes[transferbox][transferpointnumber];
                neighbourhoodforwarddisplacements[b]=forwardvelocitiesassignedtoboxes[transferbox][transferpointnumber];
                neighbourhooddistances[b]=neighbourhoodordering[b][0];
                }
            }
            for(b=0;b<dim;b++){
                forwardvelocity(b)=0;
                backwardvelocity(b)=0;
                for(c=0;c<dim;c++){
                    covarianceforwardforward(b,c)=0;
                    covariancebackwardbackward(b,c)=0;
                    covarianceforwardbackward(b,c)=0;
                }
            }
            sumweight=0;
            maxdistance=0;
            b=0;
            while(sumweight<nrneighbours){                                                                                                                //calculate primary sums which lead to the averages/covariances from the neighbourhood
                transferposition=neighbourhoodpositions[b];
                transferforwarddisplacement=neighbourhoodforwarddisplacements[b];
                transferbackwarddisplacement=neighbourhoodbackwarddisplacements[b];
                maxdistance=neighbourhooddistances[b];
                for(c=0;c<dim;c++){
                    forwardvelocity(c)=forwardvelocity(c)+(transferforwarddisplacement[c]);
                    backwardvelocity(c)=backwardvelocity(c)+(transferbackwarddisplacement[c]);
                    for(d=0;d<dim;d++){
                        covarianceforwardforward(c,d)=covarianceforwardforward(c,d)+(transferforwarddisplacement[c]*transferforwarddisplacement[d]);
                        covariancebackwardbackward(c,d)=covariancebackwardbackward(c,d)+(transferbackwarddisplacement[c]*transferbackwarddisplacement[d]);
                        covarianceforwardbackward(c,d)=covarianceforwardbackward(c,d)+(transferforwarddisplacement[c]*transferbackwarddisplacement[d]);
                    }
                }
                sumweight++;
                b++;
            }
            for(c=0;c<dim;c++){                                                                                                                                         //Finish calculation of averages and covariances
                forwardvelocity(c)=forwardvelocity(c)/sumweight;
                backwardvelocity(c)=backwardvelocity(c)/sumweight;
            }
            for(c=0;c<dim;c++){ 
                for(d=0;d<dim;d++){
                    covarianceforwardforward(c,d)=(covarianceforwardforward(c,d)/(sumweight-1))-((sumweight/(sumweight-1))*forwardvelocity(c)*forwardvelocity(d));
                    covariancebackwardbackward(c,d)=(covariancebackwardbackward(c,d)/(sumweight-1))-((sumweight/(sumweight-1))*backwardvelocity(c)*backwardvelocity(d));
                    covarianceforwardbackward(c,d)=(covarianceforwardbackward(c,d)/(sumweight-1))-((sumweight/(sumweight-1))*forwardvelocity(c)*backwardvelocity(d));
                }
            }
            covariancebackwardbackwardinverse=covariancebackwardbackward.inverse();
            friction=((-1)*covarianceforwardbackward)*covariancebackwardbackwardinverse;                                                                                 //Calculate field estimates
            frictiontranspose=friction.transpose();
            drift=forwardvelocity+(friction*backwardvelocity);
            diffusiontransposediffusion=covarianceforwardforward-(friction*covariancebackwardbackward*frictiontranspose);
            diffusion=Eigen::LLT<Eigen::MatrixXd>(diffusiontransposediffusion).matrixL();                                                                                //perform Cholesky decomposition to get diffusion matrix
    
            // Calculate actual mass and resscale friction
                                                                                                    // Calculate mass based on 1/(2kTdelta_t^2)*m^T=K_tilde_transpose^(-1)*K_tilde^(-1)*(Gamma_tilde+1)
            diffusiontranspose=diffusion.transpose();
            diffusioninverse=diffusion.inverse();
            diffusiontransposeinverse=diffusiontranspose.inverse();

            for(c=0;c<dim;c++){ 
                for(d=0;d<dim;d++){
                    helpmatrix(c,d)=friction(c,d);
                    if(c==d){
                        helpmatrix(c,d)=helpmatrix(c,d)+1;
                    }
                }
            }

            masstranspose=(diffusiontransposeinverse*diffusioninverse)*helpmatrix;
            mass=masstranspose.transpose();
            massinverse=mass.inverse();
            helpmatrix=massinverse*alterfrictionmatrix*((mass*helpmatrix)*alterfrictionmatrixtranspose);            //Alter friction matrix and rescale explicite timestep dependence
            for(c=0;c<dim;c++){ 
                for(d=0;d<dim;d++){
                    friction(c,d)=helpmatrix(c,d);
                    if(c==d){
                        friction(c,d)=friction(c,d)-1;
                    }
                }
            }
            diffusion=massinverse*(alterfrictionmatrix*(mass*diffusion));                             //Alter noise amplitude and rescale explicit timestep dependence
            diffusioninverse=diffusion.inverse();
            noise=diffusioninverse*(next_position-position-drift+(friction*(position-prev_position)));                                                                           //perform propagation
            for(c=0;c<dim;c++){                                                                                                                                           //generate output
            output << position(c) << " ";
            }
            if(fieldsno==0){
                for(c=0;c<dim;c++){
                    output << drift(c) << " ";
                }
                for(c=0;c<dim;c++){
                    for(d=0;d<dim;d++){
                        output << friction(c,d) << " ";
                    }
                }
                for(c=0;c<dim;c++){
                    for(d=0;d<dim;d++){
                        if(d<=c){
                            output << diffusion(c,d) << " ";
                        }
                    }
                }
            }
            for(c=0;c<dim;c++){
                output << noise(c) << " ";
            }
            if(statusnextpoint==0){
                output << 0;
            }
            else{
                output << 1;
            }
            if(fieldsno==0){
                output << " " << sqrt(maxdistance);
            }
            output<<endl;
            waswrittenout=1;
        }
    }
}


datatest.close();
output.close();
return EXIT_SUCCESS;
}
