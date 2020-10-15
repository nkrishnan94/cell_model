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



const double xdim = 100 *pow(10,-6); //x length of simulation in micron
const double ydim = 1400*pow(10,-6); //y length of simulation in microns
const double zdim = 100*pow(10,-6); //z length of simulation in microns
unsigned long CPD_flag = 1;

//long int cell_max = int(xdim * ydim * zdim); //max nummber ofcells allowed assuming near denset packing

const int param_N = 4; // number of cell state parameters -> center x coordinate,center y coordinate,center z coordinate, radius, y force

unsigned int tStop = 40000;
unsigned int CPD_time = 0;
float dt =2;



//parameters:
float R = 4.5*pow(10,-6); //um
float a = 2.5; //um
float Kc= 2.2*pow(10,-8);//*pow(10,0); N/um
float Ka= 2.2*pow(10,-8);//*pow(10,0); N/um
//float K= 2.2*pow(10,-10);//*pow(10,0); N/um

float s_c = 1*pow(10,6);

float s_a = .05*pow(10,6);

float fric = 0.4*pow(10,-6); //N sec. um
//float Dc = 0.01*pow(10,-6);
float Dc=0.0001;
float delc = 1*2*a*R; //um
float deld = 2.0*2*a*R; //um
float delca = 1.4*a*R; //um
float delda=2.0*a*R; //um
double vel_thresh = .5*pow(10,-6);



const int xcells = 1;
const int ycells = 80;
const int zcells=6;
const int cell_num  = int(xcells*ycells*zcells);


//const int cell_num = 1;
int record_time =int(50/dt);


int bond_check(long double cells_arr[cell_num][param_N], int x, int y,int del_c, int del_d, int del_ca, int del_da ){


    float dist=0;
    int return_int;
    for(int c=0; c<3;c++){
        //cout <<cells[ind][c] <<", " <<cells[j][c]<<endl;
        dist+=pow(cells_arr[x][c]-cells_arr[y][c],2);
        
    }
    
    if((x!=y)&&(pow(dist,.5)<delc)){
        return_int=1;
        
        
    }
       
   	if((x!=y)&&(pow(dist,.5)>deld)){
        return_int=0;
       
       
   	}
      

	if((x==y)&&(pow(dist,.5)>delda)){
	    return_int=0;
	  
	  
	}
    
	if((x==y)&&(pow(dist,.5)<delca)){
	    return_int=1;
	 
	 
	}
	return return_int;
		
	
} 



int main(){
    using namespace std;
    //cout<<cell_num<<endl;

    //cout<<cell_num<<"\n";


    /*if (cell_num > cell_max){

        std::cout << "maximum initial density exceeded" <<"\n";
        exit(EXIT_FAILURE);



    }*/
    
    default_random_engine generator;
    normal_distribution<double> distribution(0,1*pow(10,-6));

    long double cells[cell_num][param_N] = {0}; //array containing centers, radius
    long int bonds[cell_num][cell_num] = {0};
    
    long double x_space =  (double(xdim)-double(a*R))/double(xcells);
    long double y_space =  (double(ydim)-double(a*R))/double(ycells);
    long double z_space = (double(zdim)-double(a*R))/double(zcells);
    //long double y_space= 2*1.3*a*R;
    cout<< x_space<<endl;
    cout<< y_space<<endl;
    cout<< z_space<<endl;
    
    //cout<<x_space<<endl;

    int cell_count =0;
    //cells[0][0]= xdim/2;
    //cells[0][1] = ydim/2;
    //cells[0][2]= 1.3*a*R;
    //cells[0][3] = R;
    

    

    
    for(int z =0; z<zcells;z++){
        for(int y =0; y<ycells;y++){
            for(int x =0; x<xcells;x++){
                if(cell_count<cell_num){
                    cells[cell_count][0]= double(xdim/2+x*x_space);
                    cells[cell_count][1] = double(a*R+y*y_space);
                    cells[cell_count][2]= double(a*R+z*z_space);
                    cells[cell_count][3]= double(R);
                    
        //cout << cells[cell_count][2]<<endl;

                    cell_count+=1;

                }

            }
        }
    }








  

    ofstream fprof,fvels;
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
    fvels.open("CPD/vel_"+date_time.str()+"_.txt");
    //vector <double> fvels;

 
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
        
        /*if(t<int(CPD_time/dt)){
            
            Dc = .1;
            K= 2.2*pow(10,-9);
            delc = 1.05*2*a*R; //um
            
        } else{
            
            Dc = .01;
            K= 2.2*pow(10,-8);
            delc = 1.4*2*a*R; //um
            
        }*/
        
        cout <<t<<endl;

        //refresh force vector
        long double forces[cell_num][3] = {0};
        long int bonds[cell_num][3] = {0};
        
        long int cell_num_arri[cell_num] = {0};
        long int cell_num_arrj[cell_num] = {0};
        for(int i=0; i<cell_num;i++){
            
            cell_num_arri[i]=i;
            cell_num_arrj[i]=i;
        }
        random_shuffle(begin(cell_num_arri), end(cell_num_arrj));
        random_shuffle(begin(cell_num_arri), end(cell_num_arrj));

        for(int i = 0; i <cell_num; i++){
            //int ind = cell_ind[i];
            for (int j =0; j <cell_num; j++){
                int ai=cell_num_arri[i];
                int aj=cell_num_arrj[j];
                bonds[ai][aj] = bond_check(cells,ai,aj,delc,deld,delca,delda);


                if((ai!=aj)&&(bonds[ai][aj]==1)){
					float dist=0;
					for(int c=0; c<3;c++){
					//cout <<cells[ind][c] <<", " <<cells[j][c]<<endl;
						dist+=pow(cells[ai][c]-cells[aj][c],2);

					}



                    if(pow(dist,.5)<delc){

                        float xij = a*(cells[ai][3]+cells[aj][3]) - pow(dist,.5);
                        //int sign = -xij/abs(xij);

                        float fij =Kc*xij*tanh(s_c*abs(xij));

                        for(int c=0; c<3;c++){
                            forces[ai][c] += fij *(( cells[ai][c]-cells[aj][c]) /pow(dist,.5));
                            
                        }
                        //cout<<setprecision(11)<<forces[i][2]<<endl;
                        //cout<<setprecision(11)<<cells[i][2]<<endl;
                        //cout<<setprecision(11)<<cells[i][3]<<endl;


                    }

                }
                if( (ai==aj) && (t>int(CPD_time/dt)) && (bonds[ai][aj]==1) ){

                    float xij = a*(cells[ai][3]) - cells[ai][2];
                    float fij = Ka*xij*tanh(s_a*abs(xij));
                    forces[ai][2] +=fij ;
                    

                }

            }

        }

        //forces[cell_num][3] = {0};
        //cout<<forces[0][0]<<endl;
        //cout<< dt*forces[0][0]/fric + pow(2*Dc*dt,.5)*distribution(generator)<<endl;
        //cout<< cells[0][0] <<endl;

        double x;
        double y;
        double z;
        double dx=0;
        int cell_thresh =0;
        int new_cell_num=0;

        long int cell_num_arrk[cell_num] = {0};

        long double noise_terms[cell_num][3] = {0};
        for(int k=0; k<cell_num;k++){
            
            cell_num_arrk[k]=k;

            for(int c = 0; c<3;c++){

            	noise_terms[k][c] = distribution(generator);



            }
 
        }


        random_shuffle(begin(cell_num_arrk), end(cell_num_arrk));
        random_shuffle(begin(noise_terms), end(noise_terms));


        
        for(int k = 0; k <cell_num; k++){
        	int i = cell_num_arrk[k];

            if (cells[i][2]>-1){
                //cout<<pow(2*Dc*dt,.5)*distribution(generator)<<endl;

                //cout<<i<<endl;
                x = cells[i][0];
                y = cells[i][1];
                z = cells[i][2];
                //what about if theyre not meeting these conditions. you fucking idiot.....
                /*if((cells[i][0] + dt*forces[i][0]/fric + pow(2*Dc*dt,.5)*distribution(generator) < R )||(cells[i][0] + dt*forces[i][0]/fric + pow(2*Dc*dt,.5)*distribution(generator) + R>xdim )){
                    forces[i][0] = -forces[i][0];
                    forces[i][0] = 0;
                    

                }
                    

                
                if((cells[i][1] + dt*forces[i][1]/fric + pow(2*Dc*dt,.5)*distribution(generator) < R )||(cells[i][1] + dt*forces[i][1]/fric + pow(2*Dc*dt,.5)*distribution(generator) + R>ydim )){
                    forces[i][1] = -forces[i][1];
                    forces[i][1] = 0;
                    

                }
                    

                
                if(((cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator) < R )||(cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator) + R>zdim )) ){
                    forces[i][2] = -forces[i][2];
                    //cells[i][2] = cells[i][2] - dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator);
                    //forces[i][2] = 0
                    

                }*/
                
                
                if(((cells[i][1] + dt*forces[i][1]/fric + pow(2*0*dt,.5)*distribution(generator) < a*R )||(cells[i][1] + dt*forces[i][1]/fric + pow(2*0*dt,.5)*distribution(generator) >(ydim - a*R) ) )){
                    forces[i][1] = -2*forces[i][1];
                    //forces[i][1] =0;

                }
                if(((cells[i][2] + dt*forces[i][1]/fric + pow(2*0*dt,.5)*noise_terms[i][2] < a*R ) )){
                    cells[i][2] = a*R+abs( dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*noise_terms[i][2] );
                    //cells[i][2] = a*R;
                    //forces[i][1] =0;

                } else{
                    cells[i][2] = cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*noise_terms[i][2];
                    
                }

                
                cells[i][0] = cells[i][0] + dt*forces[i][0]/fric + pow(2*0*dt,.5)*noise_terms[i][0];
                cells[i][1] = cells[i][1] + dt*forces[i][1]/fric + pow(2*Dc*dt,.5)*noise_terms[i][1];
                //cells[i][2] = cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator);
                //cells[i][2] = cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator);
                
                
                

               /* if((cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator) < R) ){
                    
                    cells[i][2] =R;
                                    
                    
                    
                } else{
                    cells[i][2] = cells[i][2] + dt*forces[i][2]/fric + pow(2*Dc*dt,.5)*distribution(generator);
                    
                }*/
                    
                

                new_cell_num+=1;
                dx= pow(pow(x-cells[i][0],2)+pow(y-cells[i][1],2)+pow(z-cells[i][2],2),.5);
                if (dx>vel_thresh){
                    cell_thresh+=1;
                
                }
            }

        }


        

        


        if ((t % record_time)==0){
            ostringstream strT;
        
            strT << int(t*dt);

            string proftName = "prof_T_" + strT.str() + "_"+ date_time.str() + ".txt";
            //string yforceName = "force_T_" + strT.str() + "_"+ date_time.str() + ".txt";
            ofstream fproft;
            fproft.open("CPD/"+proftName);
            for(int i = 0; i <cell_num; i++){

                fproft<< setprecision(11)<< fixed << i << ", " << cells[i][0]*pow(10,6) << ", " << cells[i][1]*pow(10,6) << ", " << cells[i][2]*pow(10,6)<<", " << cells[i][3]*pow(10,6)<<", " << forces[i][0]*pow(10,6)<< ", "<<forces[i][1]*pow(10,6)<< ", "<<forces[i][2]*pow(10,6)<< endl;
                
                
            }
            fvels << setprecision(11)<< t <<", "<< cell_thresh << endl;

        }

        /*if (t == int(CPD_time/dt) ){
            if (CPD_flag==1){

                for(int i = 0; i <cell_num; i++){
                    if( (cells[i][1]>ydim/2-400*pow(10,-6)) &&(cells[i][1]<ydim/2+400*pow(10,-6))  &&(cells[i][2] < 40*pow(10,-6)) ){

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



