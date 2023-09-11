#include  <iostream>
#include  <fstream>
#include  <string>
#include <vector>
#include <array>
#include <cmath>   // <-- include cmath here
#include <chrono>

using namespace std;

int main(int argc, const char* argv[]) {

	int nx = std::atoi(argv[1]); //484
	int ny = std::atoi(argv[2]); //187

	int startfl=std::atoi(argv[3]);
	int endfl=std::atoi(argv[4]);
	int regionname= std::atoi(argv[5]);
	string path= argv[6];
	
	bool filepresent=true;
	//  creating X Y grids

	double dx = 1.0/(nx-1);
	double dy = 1.0/(ny-1);
	double dxdy = dx * dy;
	//=============================================
	// define the X Y grids
	double Xgrid[nx];
	for (int i=0; i<nx;i++){
		Xgrid[i] =  double(i)*dx;
	//    cout <<Xgrid[i]<<"\n";
	}

	double Ygrid[ny];
	for (int j=0; j<ny;j++){
		Ygrid[j] = double(j)*dy;
	//    cout << Ygrid[j]<<"\n";
}
//============================================
// create 2D array of Bx and By

double* Bx[nx];
for (int i=0; i<nx; i++){
    Bx[i] = new double[ny];
}
double* By[nx];
for (int i=0; i<nx; i++){
    By[i] = new double[ny];
}
//============================================
// fill the 2D array

const float a = 1.0;
double X, Y, X_1, Y_1, r, r_X_1, r_Y_1;
float pwr = 3.0/2.0;

int ry = ny - 1;
int rx = nx - 1;

// Start measuring time

auto begin = std::chrono::high_resolution_clock::now();
    
string filename;   // Name of the file
//cout<<"Enter filename:";
//cin>>filename;

string line;   // To read each line from code

//cout << "startfl" << startfl<<endl;
//cout << "endfl" << endfl<<endl;

for(int fN=startfl; fN<endfl;fN++){
//============================================
//  read bz files  
	cout << "path = " << path <<endl;

    int i=0;    // Variable to keep count of each line
    string arr[(nx*ny)+1];  // array to store each line
	
	filename= "bz_" + to_string(regionname) + "_"+to_string(fN) + ".txt";
	cout << filename<<endl;
    ifstream mFile (path+"/" + filename);   
    if(mFile.is_open()) 
    {
        while(!mFile.eof())
        {
            getline(mFile, line);
            arr[i]= line;
            i++;
        }
        mFile.close();
    }
    else{
		filepresent=false;
        cout<<"Couldn't open the file\n"; }

	double* bz[nx];
 
    for(int i=0;i<nx;i++)
    { bz[i] = new double[ny]; }

    // creat 2D bz 
    int indx=0;
    for(int j=0;j<ny;j++){
		for (int i=0; i<nx; i++){
			if(arr[i].size()!=0){
				bz[i][j] = stod(arr[indx]);
				indx++; }
		else{continue;}
        
    }}
//===============================================================
// calculate potential field interpolation using Green's function
//cout << "calculating" <<endl;
//int cj = 0;
if(filepresent==true){
cout << "calculating.." <<endl;
for (int i = 0; i<nx/2; i++){
//	cout << "i = "<< i <<endl;
//	cout << "cj = "<< cj <<endl;
//	cj = 0;
    for (int j=0; j<ny/2; j++){
//		cj++;
        for (int i1=0; i1<nx;i1++){
            
            X = Xgrid[i] - Xgrid[i1];
            X_1 = Xgrid[i+1] - Xgrid[i1];
            
            for (int j1=0; j1<ny;j1++){
                Y = Ygrid[j] - Ygrid[j1];
                Y_1 = Ygrid[j+1] - Ygrid[j1];

                r = pow(pow(X,2) + pow(Y,2) + pow(a,2),pwr);
                r_X_1 = pow(pow(X_1,2) + pow(Y,2) + pow(a,2),pwr);
                r_Y_1 = pow(pow(X,2) + pow(Y_1,2) + pow(a,2),pwr);

                // fill top and bottom of x-axis (rows) in the matric                    
 
                Bx[i][j] += bz[i1][j1] * X * dxdy/r;
                Bx[rx-i][ry-j] += ( bz[rx-i1][ry-j1] * -1*X * dxdy) / r;
 
//                cout << "Bx[i][j] = "<< j1 <<endl;				
 
                // if lx is even number
                if (nx%2 == 0){             
//                  ##  fill rows from right and left sids of the matrix 

                    if ((i+1) != nx/2){         //## i is not at the middle of lx
                        Bx[i][ry-j] += ( bz[i1][ry-j1] * X * dxdy ) / r ;         // # fill rows from right
                        Bx[rx-i][j] += ( bz[rx-i1][ry-j1] * -1*X * dxdy ) / r;   }     //# fill rows from left
                    else{
                        Bx[i][ry-j] += ( bz[i1][ry-j1] * X * dxdy ) / r;
                        Bx[i+1][j] +=  ( bz[i1][j1] * X_1 * dxdy ) / r_X_1;   }
                        
                }
                else{                           // # if lx is odd number
                    Bx[i][ry-j] += ( bz[i1][ry-j1] * X * dxdy  ) / r ;
                    Bx[rx-i][j] += ( bz[rx-i1][j1] * -1*X * dxdy ) / r ; }
                    
                    
                // # fill the middle of x-axis in the matrix if lx is odd number'''
                if (((i+1)==nx/2) && (nx%2 != 0) ){
                    Bx[i+1][j] +=  ( bz[i1][j1] * X_1 * dxdy ) / r_X_1 ;
                    Bx[i+1][ry-j] +=  ( bz[i1][ry-j1] * X_1 * dxdy ) / r_X_1 ; }
                    
                    
                // # fill the middle of y-axis in the matrix if ly is odd number'''
                if (ny%2 !=0){
                    if ((j+1) == ny/2){
                        Bx[i][j+1] +=  ( bz[i1][j1] * X * dxdy ) / r_Y_1 ;
                        Bx[rx-i][j+1] +=  (bz[rx-i1][j1] * -1*X * dxdy ) / r_Y_1 ;    }}

                // ## fill the center of the matrix if ly and lx are odd numbers                    
                if ((ny%2 !=0) && (nx%2 !=0) && ((i+1)==nx/2)  && ((j+1)==ny/2)){
                    
                        Bx[rx/2][ry/2] +=  ( bz[rx/2][ry/2] * X_1 *dxdy ) / pow(pow(X_1,2) + pow(Y_1,2) + pow(a,2),pwr) ; }

// ============================================ By =========================================================

                // fill top and bottom of x-axis (rows) in the matric                    
 
                By[i][j] += (bz[i1][j1] * Y * dxdy)/r;
                By[rx-i][ry-j] += ( bz[rx-i1][ry-j1] * -1*Y * dxdy) / r;

                // if lx is even number
                if (nx%2 == 0){             
//                  ##  fill rows from right and left sids of the matrix 

                    if ((i+1) != nx/2){         //## i is not at the middle of lx
                        By[i][ry-j] += ( bz[i1][ry-j1] * -1*Y * dxdy ) / r ;         // # fill rows from right
                        By[rx-i][j] += ( bz[rx-i1][ry-j1] * Y * dxdy ) / r;   }     //# fill rows from left
                    else{
                        By[i][ry-j] += ( bz[i1][ry-j1] * -1*Y * dxdy ) / r;
                        By[i+1][j] +=  ( bz[i1][j1] * Y * dxdy ) / r_X_1;   }
                }
                else{                           // # if lx is odd number
                    By[i][ry-j] += ( bz[i1][ry-j1] * -1*Y * dxdy  ) / r ;
                    By[rx-i][j] += ( bz[rx-i1][j1] * Y * dxdy ) / r ; }
                    
                    
                // # fill the middle of x-axis in the matrix if lx is odd number'''
                if (((i+1)==nx/2) && (nx%2 != 0) ){
                    By[i+1][j] +=  ( bz[i1][j1] * Y * dxdy ) / r_X_1 ;
                    By[i+1][ry-j] +=  ( bz[i1][ry-j1] * -1*Y * dxdy ) / r_X_1 ; }
                    
                    
                // # fill the middle of y-axis in the matrix if ly is odd number'''
                if (ny%2 !=0){
                    if ((j+1) == ny/2){
                        By[i][j+1] +=  ( bz[i1][j1] * Y_1 * dxdy ) / r_Y_1 ;
                        By[rx-i][j+1] +=  (bz[rx-i1][j1] * Y_1 * dxdy ) / r_Y_1 ;    }}

                // ## fill the center of the matrix if ly and lx are odd numbers                    
                if ((ny%2 !=0) && (nx%2 !=0) && ((i+1)==nx/2)  && ((j+1)==ny/2)){
                    
                        By[rx/2][ry/2] +=  ( bz[rx/2][ry/2] * Y_1 *dxdy ) / pow(pow(X_1,2) + pow(Y_1,2) + pow(a,2),pwr) ;


                }
//                cout <<  Bx[i][j]<< "\n";
            }
        }
        
    }
}

}
// Stop measuring time and calculate the elapsed time
auto end = std::chrono::high_resolution_clock::now();
auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);


printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);

//int rows = sizeof(Bx);///sizeof(Bx[0]);
//int cols = sizeof(Bx[0])/sizeof(Bx[0][0]);

//cout << "rows:" << std::extent_v<decltype(Bx), 0> << '\n';

/*
cout <<  rx<< '\t'<<ry<<"\n";
// cout <<  Bx[i][j]<< "\n";x[0][0]<< "\n";
cout <<  Bx[0][0]<< "\n";
cout <<  Bx[2][1]<< "\n";
cout <<  By[0][0]<< "\n";
cout <<  By[2][3]<< "\n";
cout <<  Bx[rx][ry-1]<< "\n";
cout <<  Bx[rx][ry]<< "\n";
//cout <<  Bx[0][6]<< "\n";
//cout <<  Bx[0][7]<< "\n";
//cout <<  Bx[0][8]<< "\n";

*/

    
// ==== print matix in a file =====

	string namex=path+"/Bxp_" + to_string(regionname) + "_"+to_string(fN)+".txt";
	ofstream Bxp(namex);

	for (int j=0;j<ny;j++)
	{
		for (int i=0;i<nx;i++)
		{   
			Bxp << to_string(Bx[i][j]) << endl; // behaves like cout - cout is also a stream
		}
//		Bxp << "\n";
}  
    Bxp.close();


	string namey=path+"/Byp_" + to_string(regionname) + "_"+to_string(fN)+".txt";	
	ofstream Byp(namey);

	for (int j=0;j<ny;j++)
	{
		for (int i=0;i<nx;i++)
		{   
			Byp << to_string(By[i][j]) << endl; // behaves like cout - cout is also a stream
		}
//		Byp << "\n";
}  
    Byp.close();
  
}




}
