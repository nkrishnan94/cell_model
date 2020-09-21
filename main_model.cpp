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



const int xdim = 90; //x length of simulation in microns
const int ydim = 1000; //y length of simulation in microns
const int zdim = 150; //z length of simulation in microns
unsigned long CPD_flag = 1;

//long int cell_max = int(xdim * ydim * zdim); //max nummber ofcells allowed assuming near denset packing

const int param_N = 4; // number of cell state parameters -> center x coordinate,center y coordinate,center z coordinate, radius, y force

unsigned int tStop = 2000;
unsigned int CPD_time = 0;
float dt =2;



//parameters:
float R = 4.125; //um
float a = 2.5; //um
float K= .088;//*pow(10,0); N/um
float s_b = 1;

float fric = 0.4;//*pow(10,-6); N sec. um
float Dc = 0.01;
float delc = 1.4*a*R; //um
float deld =1.4*2*a*R; //um 
float delca = 1*a*R; //um
float delda=1.8*a*R; //um


const int cell_num  = 13500;
int xcells = 9;
int ycells = 100;
int zcells= 15;
int record_time =int(10/dt);



int main(){
	using namespace std;
	//cout<<cell_num<<endl;

	//cout<<cell_num<<"\n";


	/*if (cell_num > cell_max){

		std::cout << "maximum initial density exceeded" <<"\n";
		exit(EXIT_FAILURE);



	}*/
	default_random_engine generator;
	normal_distribution<double> distribution(0,1);

	long double cells[cell_num][param_N] = {0}; //array containing centers, radius
	
	long double x_space =  (double(xdim)-double(2*R))/double(xcells);
	long double y_space =  (double(ydim)-double(2*R))/double(ycells);
	long double z_space = (double(zdim)-double(2*R))/double(zcells);
	//cout<<x_space<<endl;

	int cell_count =0;

	for(int z =0; z<zcells;z++){
		for(int y =0; y<ycells;y++){
			for(int x =0; x<xcells;x++){
				if(cell_count<cell_num){
					cells[cell_count][0]= double(R+x*x_space);
					cells[cell_count][1] = double(R+y*y_space);
					cells[cell_count][2]= double(R+z*z_space);
					cells[cell_count][3]= double(R);

					cell_count+=1;

				}

			}
		}
	}

    for(int i = 0; i <cell_num; i++){
    	if( (cells[i][1]>int(ydim/2)-300)&&(cells[i][1]<int(ydim/2)+300) &&(cells[i][2] <abs(int(ydim/2)-i)/8 ) ) {

    		cells[i][2] =-30;
    		cells[i][3] = 0;
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

	fprof.open("CPD/prof_init_" + date_time.str()+"_.txt");

 
	/*establish radii
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
	}*/

	for(int i = 0; i < cell_num; i++){



		fprof <<setprecision(8)<< fixed << i << ", " << cells[i][0] << ", " << cells[i][1] << ", " << cells[i][2] <<", " << cells[i][3]<<endl;

	}



	for(int t =0;t<int(tStop/dt);t++){
		cout <<t<<endl;

		//refresh force vector
		long double forces[cell_num][3] = {0}; 

	    for(int i = 0; i <cell_num; i++){
	    	//int ind = cell_ind[i];
	    	for (int j =0; j <cell_num; j++){
	    		if(i!=j){
	    			float dist=0;
	    			for(int c=0; c<3;c++){
	    				//cout <<cells[ind][c] <<", " <<cells[j][c]<<endl;
	    				dist+=pow(cells[i][c]-cells[j][c],2);
	    			}
	    			if(pow(dist,.5)<delc){

		    			float xij = a*(cells[i][3]+cells[j][3]) - pow(dist,.5);
		    			int sign = xij/abs(xij);

		    			float fij = K*xij*tanh(s_b*abs(xij));

		    			for(int c=0; c<3;c++){
		    				forces[i][c] += fij * sign*(( cells[i][c]-cells[j][c]) /pow(dist,.5));
		    				
		    			}	


	    			}

	    		}
	    		if(i==j){
	    			if(cells[i][2]<delca){
	    				float xij = a*(cells[i][3]) - cells[i][2];
	    				float fij = K*xij*tanh(s_b*abs(xij));
	    				forces[i][2] += fij;
	    			} 


	    			

	    		}	    		

    		}

    	}
    	//cout<<forces[0][0]<<endl;
    	//cout<< dt*forces[0][0]/fric + pow(2*Dc*dt,.5)*distribution(generator)<<endl;
    	//cout<< cells[0][0] <<endl;

    	
		for(int i = 0; i <cell_num; i++){

			//cout<<i<<endl;
			if(cells[i][0] + dt*forces[i][0]/fric + pow(2*Dc*dt,.5)*distribution(generator) < 0 ){
				cells[i][0] = xdim + cells[i][0] + dt*forces[i][0]/fric + pow(2*Dc*dt,.5)*distribution(generator);

			}
			if(cells[i][0] + dt*forces[i][0]/fric + pow(2*Dc*dt,.5)*distribution(generator) >xdim){
				cells[i][0] = cells[i][0] + dt*forces[i][0]/fric + pow(2*Dc*dt,.5)*distribution(generator)-xdim;

			}

			if(cells[i][1] + dt*forces[i][1]/fric + pow(2*Dc*dt,.5)*distribution(generator) <(R)){
				cells[i][1] = R+abs(cells[i][1] + dt*forces[i][1]/fric + pow(2*Dc*dt,.5)*distribution(generator) );

			}

			if(cells[i][1] + dt*forces[i][1]/fric + pow(2*Dc*dt,.5)*distribution(generator)+R >ydim){
				cells[i][1] = cells[i][1] + dt*forces[i][1]/fric + pow(2*Dc*dt,.5)*distribution(generator)-(a*R);

			}
				
    		cells[i][2] = cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator);

		}

		

		


	    if ((t % record_time)==0){
			ostringstream strT;
			strT << int(t*dt);

			string proftName = "prof_T_" + strT.str() + "_"+ date_time.str() + ".txt";
			//string yforceName = "force_T_" + strT.str() + "_"+ date_time.str() + ".txt";
			ofstream fproft;
		    fproft.open("CPD/"+proftName);
		    for(int i = 0; i <cell_num; i++){

	    		fproft<< setprecision(11)<< fixed << i << ", " << cells[i][0] << ", " << cells[i][1] << ", " << cells[i][2] <<", " << cells[i][3]<<", " << forces[i][0]<< ", "<<forces[i][1]<< ", "<<forces[i][2]<< endl;
		    	
			}

	    }

	    /*if (t == int(CPD_time/dt) ){
	    	if (CPD_flag==1){
			    for(int i = 0; i <cell_num; i++){
			    	if( (cells[i][1]>int(ydim/2)-300)&&(cells[i][1]<int(ydim/2)+300) &&(cells[i][2] <(int(a*R)+42+abs(int(ydim/2)-i)/4 ) )){

			    		cells[i][2] =-30;
			    		cells[i][3] = 0;
			    	}
			    }
	    	}

	    }*/

    }

    clock_t c_fin = clock();
    double run_time = double(c_fin - c_init)/CLOCKS_PER_SEC;
    cout << "Finished!" << "\n";
    cout << "Finished in " << run_time << " seconds \n";
	puts (buffer);
	//cout<<setprecision(8)<< fixed <<double(rsum)/double(cell_num)<<endl;


	//fprof <<setprecision(8)<< fixed << i << ", " << cells[i][0] << ", " << cells[i][1] << ", " << cells[i][2] <<", " << cells[i][3]endl;



	return 0;
}




