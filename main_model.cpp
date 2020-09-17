#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <ctime> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



const int xdim = 220; //x length of simulation in microns
const int ydim = 220; //y length of simulation in microns
const int zdim = 220; //z length of simulation in microns

//long int cell_max = int(xdim * ydim * zdim); //max nummber ofcells allowed assuming near denset packing




const int param_N = 4; // number of cell state parameters -> center x coordinate,center y coordinate,center z coordinate, radius


//parameters:
float R = 4.125;
float a = 2.5;


const int cell_num  = 1000;
int xcells = 10;
int zcells = 10;
int ycells= 10;



int main(){
	using namespace std;
	//cout<<cell_num<<endl;

	//cout<<cell_num<<"\n";


	/*if (cell_num > cell_max){

		std::cout << "maximum initial density exceeded" <<"\n";
		exit(EXIT_FAILURE);



	}*/
	long double cells[cell_num][param_N] = {0}; //array containing centers, radius
	long double forces[cell_num][3] = {0}; 
	long double x_space =  double(xdim)/double(xcells+1);
	long double y_space =  double(ydim)/double(ycells+1);
	long double z_space = double(zdim)/double(zcells+1);
	//cout<<x_space<<endl;

	int cell_count =0;

	for(int z =0; z<zcells;z++){
		for(int y =0; y<ycells;y++){
			for(int x =0; x<xcells;x++){
				if(cell_count<cell_num){
					cells[cell_count][0]= double(.5*x_space+x*x_space);
					cells[cell_count][1] = double(.5*y_space+y*y_space);
					cells[cell_count][2]= double(.5*z_space+z*z_space);

					cell_count+=1;

				}

			}
		}
	}
	ofstream fprof;
	time_t time_start;
	clock_t c_init = clock();
	struct tm * timeinfo;
	char buffer [80];
	float ht = 0.5;
    time (&time_start);
	timeinfo = localtime (&time_start);




	strftime (buffer,80,"%F-%H-%M-%S",timeinfo);

	ostringstream date_time;
	date_time << buffer;

	fprof.open("prof_init_" + date_time.str()+"_.txt");

 

	for(int n=0;n<10;n++){
		int cell_ind[cell_num] ={0};
		for(int i = 0; i <cell_num; i++){
			cell_ind[i] =i;
		}
		random_shuffle(cell_ind, cell_ind + cell_num);

	    for(int i = 0; i <cell_num; i++){
	    	int ind = cell_ind[i];
	    	double min_dist =500;
	    	
	    	
	    	for (int j =0; j <cell_num;j++){
	    		double dist=0;

	    		for(int c=0; c< 3;c++){
	    			dist+=pow(cells[ind][c]-cells[j][c],2);
	    		}
	    		//cout<<pow(dist,.5)<<endl;
	    		if ((pow(dist,.5)<min_dist)&&(dist>0)){
	    			min_dist=double(pow(dist,.5));
	    		}

	    	}

	    	cells[i][3] = (double(min_dist)/2)/double(a);

		}
	}

	for(int i = 0; i < cell_num; i++){
		rsum+=cells[i][3];


		fprof <<setprecision(8)<< fixed << i << ", " << cells[i][0] << ", " << cells[i][1] << ", " << cells[i][2] <<", " << cells[i][3]<<endl;

	}
	//cout<<setprecision(8)<< fixed <<double(rsum)/double(cell_num)<<endl;


	//fprof <<setprecision(8)<< fixed << i << ", " << cells[i][0] << ", " << cells[i][1] << ", " << cells[i][2] <<", " << cells[i][3]endl;



	return 0;
}




