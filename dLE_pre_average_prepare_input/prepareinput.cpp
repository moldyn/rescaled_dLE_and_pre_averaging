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

int main(int argc, char* argv[]){
    
for (int i = 1; i < argc; ++i) {
  if(std::string(argv[i]) == "-h"){
    cout << " " <<endl;
    cout<< "USAGE: ./prepareinput [OPTIONS]" << endl;
    cout << " " <<endl;
    cout << "INPUT: trajectory x, forward displacements v, backward displacements w (x1 x2, ..., xn, v1, ..., vn, w1, ..., wn), file specifying the binning" <<endl;
    cout << "OUTPUT: averages of vn, wn, v_n*v_n, w_n*w_n, v_n,w_n for the binning defined by the options" <<endl;
    cout << " " <<endl;
    cout<< "The program performs a binning of the observables which are relevant for the dLE." << endl;
    cout<< "As input, the program need the trajectory (option -i), the forward displacements vi" << endl;
    cout<< "and the backward displacements wi as (x1 x2, ..., xn, v1, ..., vn, w1, ..., wn)" << endl;
    cout<< "and a file with two lines representing the lower and the upper borders" << endl;
    cout<< "for the binning and n colums where each column describes one dimension (option -f)." << endl;
    cout<< "The line ordering is unimportant, upper and lower borders are detected, but" << endl;
    cout<< "the first column must be the borders for x1, the second for x2 and so on." << endl;
    cout<< "Additionally, the dimensionality of the trajectory must be given (option -d)." << endl;
    cout<< "The output is the result of the binning of the trajectory and the calculations" << endl;
    cout<< "of the obserables (option -o). It is given in the following order " << endl;
    cout<< "(x1_av,bin, x2_av,bin ..., xn_bin, v1_av, ..., w1_av, (v1*v1)_av, (v1*v2)_av, ..., (w1*w1)_av,..., (v1*w1)_av..., N)." << endl;
    cout << " " <<endl;
    cout<< "The option -s should prevent memory problems. If one coarse bin contains more points than it is wished (option -b)" << endl;
    cout<< "this bin is subdivided in smallbins=(N_bin/N_max)^{1/dim} per dimension. In case the width of this new bins are below a certain threshold" << endl;
    cout<< "(option -wmin), smallbins is decresed again so that the threshold is not undercut. The option -wmax defines the minimal resolution of the" << endl;
    cout << "binning in the same way. In case the width of the bins defined by option -b are above -wmax, smallbins is increased." <<endl;
    cout << " " <<endl;
    cout << " " <<endl;
    cout << "OPTIONS" <<endl;
    cout << " " << " " << "-h show these lines" << endl;
    cout << " " << " " << "-i name input trajectory [default input]" << endl;
    cout << " " << " " << "-o name output trajectory [default input.binned]" << endl;
    cout << " " << " " << "-f name file which contains the binning borders [default input.borders]" << endl;
    cout << " " << " " << "-s number of coarse bins per dimension [default: 400]" << endl;
    cout << " " << " " << "-b number of maximal points per bin  [default: 100]" << endl;
    cout << " " << " " << "-wmin minimal width of bin after subdivison defined by -b [default:0.00001]" << endl;
    cout << " " << " " << "-wmax maximal width of bin after subdivison defined by -b [default:0.0001]" << endl;
    cout << " " << " " << "-d dimensionality input [default: 1] ATTENTION: contradiction to real dimensionality (trajectory) is not detected by program " << endl;
    cout << " " <<endl;
    cout << " " <<endl;
    exit(0);
  }
}

//Input variables
int dim,maxpointsbin,maxbins,minbins,coarsebins,bins,numbercoarsebinsatall,numbercomponentsaverageviwi,numbercomponentsnumberentries,numbercomponentscovariances;
double minimalwidth,maximalwidth,smallestrange;

//Count variables
int a,b,c,d,e;

//local variables
std::vector<int> binindex;								    //determine by this vector to which bin a trajectory point belongs 
std::vector<int> binindexcoarse;                             //this vector stores the "coarse bin"
std::vector<double> borders;                                //store binning borders in vector
std::vector<double> actualcoarsebinborders;                  //store binning borders of actual coarse grained bin in vector
std::vector<double> actualpoint;                            //store actual trajectory point in vector
std::vector<double> averagexi;                              //store average point position in bin
std::vector<double> averagevi;                              //store average forward displacement in vector
std::vector<double> averagewi;                              //store average backward displacement in vector
std::vector<double> numberentries;                          //store number points per bin in vector
std::vector<double> averagevivi;                          //store Cov(vi,vi) in vector
std::vector<double> averagewiwi;                          //store Cov(wi,wi) in vector
std::vector<double> averageviwi;                          //store Cov(vi,wi) in vector
double borderactual;
int helpdeterminecoarsebin,helpdeterminebin;
double checkborders,checkwidth,checknumbereintries;
int checkinbin;
int actualbin,actualcoarsebin;
int usedpoints=0;

//preordering of points which should be binned
std::vector<std::vector<std::vector<double>>> pointsassigenedtocoarsebins;
std::vector<int> numberentriespointsassigendtocoarsebins;

maxpointsbin=100;
maxbins=1;
smallestrange=0;
dim=1;
coarsebins=400;
bins=1;
minimalwidth=0.00001;
maximalwidth=0.0001;
std::string name="input";
std::string nameoutput="input.binned";
std::string nameborders="input.borders";


    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-i") {
            if (i + 1 < argc) { 							// Make sure we aren't at the end of argv!
		name=argv[i+1]; 						// Increment 'i' so we don't get the argument as the next argv[i].
            }
        }
        else if (std::string(argv[i]) == "-o") {
            if (i + 1 < argc) {
		nameoutput=argv[i+1];
            }
        }
        else if (std::string(argv[i]) == "-f") {
            if (i + 1 < argc) {
		nameborders=argv[i+1];
            }
        }
        else if (std::string(argv[i]) == "-d") {
            if (i + 1 < argc) {
		dim=(int)(atof(argv[i+1]));
            }
        }
        else if (std::string(argv[i]) == "-b") {
            if (i + 1 < argc) {
		maxpointsbin=(int)(atof(argv[i+1]));
            }
        }
        else if (std::string(argv[i]) == "-wmin") {
            if (i + 1 < argc) {
		minimalwidth=atof(argv[i+1]);
            }
        }
        else if (std::string(argv[i]) == "-wmax") {
            if (i + 1 < argc) {
		maximalwidth=atof(argv[i+1]);
            }
        }
        else if (std::string(argv[i]) == "-s") {
            if (i + 1 < argc) {
		coarsebins=(int)(atof(argv[i+1]));
            }
        }
    }
    
if(minimalwidth>maximalwidth){
    checkwidth=minimalwidth;
    minimalwidth=maximalwidth;
    maximalwidth=checkwidth;
}
    
numbercoarsebinsatall=pow(coarsebins,dim);

cout << "dim" << " " << dim << " " << "nr.coarsebins per dim" << " " << coarsebins << " " << "nr.coarsebins" << " " << numbercoarsebinsatall << endl;

for(a=0;a<numbercoarsebinsatall;a++){
    pointsassigenedtocoarsebins.push_back(std::vector<std::vector<double>>());
    numberentriespointsassigendtocoarsebins.push_back(0);
}
    
for(a=0;a<dim;a++){
  binindex.push_back(0);                                        //initialise it with zeros
}
for(a=0;a<dim;a++){
  binindexcoarse.push_back(0);
}
for(a=0;a<(2*dim);a++){
  borders.push_back(0);  
}
for(a=0;a<(2*dim);a++){
  actualcoarsebinborders.push_back(0);  
}
for(a=0;a<(3*dim);a++){
  actualpoint.push_back(0);  
}
    
ifstream inputborders;
inputborders.open(nameborders,ios::in);

a=0;

while(inputborders>>borderactual){                                  //Store borders
 borders[a]=borderactual;
 a++;
}

inputborders.close();

for(a=0;a<dim;a++){                                                 //Make sure that lower borders are lower borders
    if(borders[a]>borders[a+dim]){
        checkborders=borders[a];
        borders[a]=borders[a+dim];
        borders[a+dim]=checkborders;
    }
}

smallestrange=borders[dim]-borders[0];                              //Determine maximal number of small bins if the large bin is subdivided
for(a=1;a<dim;a++){
    if(borders[a+dim]-borders[a]<smallestrange){
        smallestrange=borders[a+dim]-borders[a];
    }
}
maxbins=(int)(smallestrange/(minimalwidth*coarsebins));
minbins=(int)(smallestrange/(maximalwidth*coarsebins));

cout << "maxbins" << " " << maxbins << " " << "minbins" << " " << minbins << endl;

ifstream data;                                                      //open input file
data.open(name,ios::in);

cout << "open" << endl;

while(data>>actualpoint[0]){                                                    //assign points to coarse bins
    for(b=1;b<3*dim;b++){
        data >> actualpoint[b];
    }
    checkinbin=0;
    for(b=0;b<dim;b++){
        if(actualpoint[b]>borders[b+dim]){
            checkinbin++;
        }
        else if(actualpoint[b]<borders[b]){
            checkinbin++;
        }
    }
    if(checkinbin==0){
        for(b=0;b<dim;b++){
            binindexcoarse[b]=(int)((actualpoint[b]-borders[b])/(borders[b+dim]-borders[b])*coarsebins);
        }
        actualcoarsebin=0;
        for(b=0;b<dim;b++){                                                                                                                 //assign vector binindex to unambigouse single integer
            actualcoarsebin=actualcoarsebin+(int)((binindexcoarse[b]*pow(coarsebins,b)));
        }
        pointsassigenedtocoarsebins[actualcoarsebin].push_back(actualpoint);
        numberentriespointsassigendtocoarsebins[actualcoarsebin]=numberentriespointsassigendtocoarsebins[actualcoarsebin]+1;
        usedpoints++;
    }
}

cout << "Assigned points" << " " << usedpoints << endl;
usedpoints=0;

data.close();

ofstream output;                                                    //open output file
output.open(nameoutput,ios::out | ios::trunc);

output << "#" << " " << "./prepareinput" << " " << "-i" << " " << name << " " << "-s" << " " << coarsebins << " " << "-wmin" << " " << minimalwidth << " " << "-wmax" << " " << maximalwidth << " " << "-b" << " " << maxpointsbin << " " << "-d" << " " << dim <<endl;
output << "#" << " ";                                               //insert comment line
for(c=0;c<dim;c++){
    output << "x"<<c << " ";
}
for(c=0;c<dim;c++){
    output << "v"<<c << " ";
}
for(c=0;c<dim;c++){
    output << "w"<<c << " ";
}
for(c=0;c<dim;c++){
    for(d=0;d<dim;d++){
        output << "v"<<c<<"*v"<<d<< " ";
    }
}
for(c=0;c<dim;c++){
    for(d=0;d<dim;d++){
        output << "w"<<c<<"*w"<<d<< " ";
    }
}
for(c=0;c<dim;c++){
    for(d=0;d<dim;d++){
        output << "v"<<c<<"*w"<<d << " ";
    }
}
output << "entries" << endl;

for(a=0;a<numbercoarsebinsatall;a++){
    if(numberentriespointsassigendtocoarsebins[a]>0){
        helpdeterminecoarsebin=a;
        binindexcoarse[0]=a%coarsebins;                                     //Determine actual coarse bin;
        for(b=1;b<dim;b++){
            helpdeterminecoarsebin=helpdeterminecoarsebin/coarsebins;
            binindexcoarse[b]=helpdeterminecoarsebin%coarsebins;
        }
        for(b=0;b<dim;b++){
            actualcoarsebinborders[b]=borders[b]+(binindexcoarse[b]*((borders[b+dim]-borders[b])/coarsebins));
            actualcoarsebinborders[b+dim]=borders[b]+((binindexcoarse[b]+1)*((borders[b+dim]-borders[b])/coarsebins));
        }
        if(numberentriespointsassigendtocoarsebins[a]>maxpointsbin){                                                                                 //determine partitioning of coarse bin
            checknumbereintries=double(numberentriespointsassigendtocoarsebins[a])/double(maxpointsbin);
            checknumbereintries=pow(checknumbereintries,double(1.0/dim));
            bins=int(checknumbereintries);
        }
        if(bins>maxbins){                                                                                                               // make sure minimal binwidth is preserved 
            bins=maxbins;
        }
        if(bins<minbins){
            bins=minbins;
        }
        cout << "coarse bin" << " " << a << " " << "bins" << " " << bins << " " << "points" << " " << numberentriespointsassigendtocoarsebins[a] << endl;
        numbercomponentsnumberentries=pow(bins,dim);                                                                                    //initialise calculated quantities
        numbercomponentsaverageviwi=pow(bins,dim)*dim;
        numbercomponentscovariances=pow(bins,dim)*pow(dim,dim);
        for(b=0;b<numbercomponentsaverageviwi;b++){
            averagexi.push_back(0);  
            }
        for(b=0;b<numbercomponentsaverageviwi;b++){
            averagevi.push_back(0);  
            }
        for(b=0;b<numbercomponentsaverageviwi;b++){
            averagewi.push_back(0);  
            }
        for(b=0;b<numbercomponentsnumberentries;b++){
            numberentries.push_back(0);  
            }
        for(b=0;b<numbercomponentscovariances;b++){
            averagevivi.push_back(0);  
            }
        for(b=0;b<numbercomponentscovariances;b++){
            averagewiwi.push_back(0);  
            }
        for(b=0;b<numbercomponentscovariances;b++){
            averageviwi.push_back(0);  
            }
        for(e=0;e<numberentriespointsassigendtocoarsebins[a];e++){
            actualpoint=pointsassigenedtocoarsebins[a][e];
            for(b=0;b<dim;b++){                                                                                                                 //transfer spacial information to integers -> which bin?
                binindex[b]=(int)((actualpoint[b]-actualcoarsebinborders[b])/(actualcoarsebinborders[b+dim]-actualcoarsebinborders[b])*bins);
                if(binindex[b]<0){
                    binindex[b]=0;
                }
                else if(binindex[b]>bins-1){
                    binindex[b]=bins-1;
                }
                }
                actualbin=0;
                for(b=0;b<dim;b++){                                                                                                                 //assign vector binindex to unambigouse single integer
                    actualbin=actualbin+(int)((binindex[b]*pow(bins,b)));
                }   
                numberentries[actualbin]=numberentries[actualbin]+1;                                                                                //count how many points per bin
                for(b=0;b<dim;b++){
                    averagexi[actualbin+(b*pow(bins,dim))]+=actualpoint[b];
                    averagevi[actualbin+(b*pow(bins,dim))]+=actualpoint[dim+b];                                                                     //first bins^dim vector elements: v_1, next bins^dim elements v_2 ...
                    averagewi[actualbin+(b*pow(bins,dim))]+=actualpoint[(2*dim)+b];                                                                 //first bins^dim vector elements: w_1, next bins^dim elements w_2 ...
                }
                for(b=0;b<dim;b++){
                    for(c=0;c<dim;c++){                                                                                                     //First bins^dim vector elements: Cov(v1,v1), next bins^dim elements Cov(v1,v2)...  
                        averagevivi[actualbin+(c*pow(bins,dim))+(b*dim*pow(bins,dim))]+=(actualpoint[dim+b]*actualpoint[dim+c]);
                        averagewiwi[actualbin+(c*pow(bins,dim))+(b*dim*pow(bins,dim))]+=(actualpoint[(2*dim)+b]*actualpoint[(2*dim)+c]);
                        averageviwi[actualbin+(c*pow(bins,dim))+(b*dim*pow(bins,dim))]+=(actualpoint[dim+b]*actualpoint[(2*dim)+c]);
                    }
                }
            } 
        for(b=0;b<numbercomponentsnumberentries;b++){                                                                                           //dividing by number entries
            if(numberentries[b]>1){
                for(c=0;c<dim;c++){
                    averagexi[b+(c*pow(bins,dim))]=averagexi[b+(c*pow(bins,dim))]/numberentries[b];
                    averagevi[b+(c*pow(bins,dim))]=averagevi[b+(c*pow(bins,dim))]/numberentries[b];
                    averagewi[b+(c*pow(bins,dim))]=averagewi[b+(c*pow(bins,dim))]/numberentries[b];
                }
                for(c=0;c<dim;c++){                                                                                                             //Calculate covariances
                    for(d=0;d<dim;d++){
                        averagevivi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))]=averagevivi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))]/numberentries[b];
                        averagewiwi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))]=averagewiwi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))]/numberentries[b];
                        averageviwi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))]=averageviwi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))]/numberentries[b];
                    }
                }
            }
        }
        for(b=0;b<numbercomponentsnumberentries;b++){                                     // Write out
            if(numberentries[b]>0){
                helpdeterminebin=b;
                binindex[0]=b%bins;                                                          //Determine actual bin;
                for(c=1;c<dim;c++){
                    helpdeterminebin=helpdeterminebin/bins;
                    binindex[c]=helpdeterminebin%bins;
                }
                for(c=0;c<dim;c++){
                    output << averagexi[b+(c*pow(bins,dim))] << " ";
                }
                for(c=0;c<dim;c++){
                    output << averagevi[b+(c*pow(bins,dim))] << " ";
                }
                for(c=0;c<dim;c++){
                    output << averagewi[b+(c*pow(bins,dim))] << " ";
                }
                for(c=0;c<dim;c++){
                    for(d=0;d<dim;d++){
                        output << averagevivi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))] << " ";
                    }
                }
                for(c=0;c<dim;c++){
                    for(d=0;d<dim;d++){
                        output << averagewiwi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))] << " ";
                    }
                }
                for(c=0;c<dim;c++){
                    for(d=0;d<dim;d++){
                        output << averageviwi[b+(d*pow(bins,dim))+(c*dim*pow(bins,dim))] << " ";
                    }
                }
                output << numberentries[b] << endl;
                usedpoints=usedpoints+numberentries[b];
            }
        }
        averagexi.clear();
        averagevi.clear();
        averagewi.clear();
        numberentries.clear();
        averagevivi.clear();
        averagewiwi.clear();
        averageviwi.clear();
    }   
}

cout << "Used points" << " " << usedpoints << endl;

output.close();
return EXIT_SUCCESS;

}
