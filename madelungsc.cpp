//#include "montecarlo.hpp"

#include <ctime>
#include <vector>
#include <cmath>
#include <fstream> 
#include <complex>
#include <string>
#include <iostream>
#include <iterator>
#include <iomanip>

using namespace std;

//typedef vector<int> Row;
//typedef vector<Row> Matrix;
//typedef vector<double> Roww;
//typedef vector<Roww> Matrixx;
//typedef vector<Matrix> Matrix3d;
typedef complex<double> dcomp;

//typedef struct {
//    int x,y,z;
//} xyz;

typedef struct {
    double x,y,z;
} xyz;

double dotproduct(vector<double> v1, vector<double> v2);
double Vsr12(int i1, int i2, vector<xyz>& h, int Nmax,  vector<xyz>& s, double s0);
double Vlr12(int i1, int i2, vector<xyz>& m, int Lmax,  vector<xyz>& s, double s0);
double zeromode( vector<xyz>& m, int Lmax, double s0);
void FormFactor( vector<xyz>& m, int Lmax, vector<xyz> & coord,int z,vector<dcomp> & Sk);
double VSelfEnergy(double s0,double Z);
dcomp VLongRange(vector<xyz>& m, int Lmax, double s0, vector<dcomp> & Sk);
double VelectronProton(double s0,double Z);
int GenerateNmodes(vector<xyz> & h, double a0, int iNtot);
int Generatekmodes(vector<xyz> & m, double a0 ,int  Ltot,double L);
void generate_coord(vector<xyz> & v2 , int z, int choice, double aLatt);



int main(void){

dcomp VLr ;
int iLmax, iNmax;
int Ltot = 500000 , iNtot = 50000;
vector<dcomp> Sk(Ltot);
double L =  5;
double Z =  pow(L,3);
int z   = 5;
int V   = pow(z,3);
int fcc = 0;

double aEff = 1/L;
double aLatt = 1/L;
double madelung,VSE,VEP,Vsr,kzeromode,Vlr,en12,d12;
double a0 = .12*L;
double da = .002;
double s0;

//Matrixx coord;
//m = Matrixx(3,Roww(Ltot));
//h = Matrixx(3,Roww(iNtot));

vector<xyz> coord(V);
vector<xyz> m(Ltot);
vector<xyz> h(iNtot);
    

generate_coord(coord,z,fcc,aLatt);

cout << "hello \n ";


if (fcc == 1){
aEff = 1/(pow(4,.333333)*L);
}

double sc = 2*aEff*(1/Z);
cout <<"sc" <<sc << "\n ";
cout<< "V" << V << "\t"<< "1/L = \t"<< 1/L <<"\t aEff:" << aEff<< "\n";

string ss = to_string(z);

ofstream ewalde,zx;

ewalde.open("madelung_friends"+ss+"sc"+".txt");

/*q.open("long_range"+ss+"sc"+".txt");
w.open("short_range"+ss+"sc"+".txt");
x.open("madelung"+ss+"sc"+".txt");*/
zx.open("energy_distance"+ss+"sc"+".txt");


for(unsigned int p = 0; p<V; p++ ){
cout<<p<<"::" <<L*coord[p].x<<":"<< L*coord[p].y<<":" << L*coord[p].z  <<"\n";
}

for(unsigned int y = 0; y < 100; y++){

a0 += da;

s0 = a0/L;
cout<< "s0 = " << s0 << "\n";

Vsr = 0;
Vlr = 0;

iLmax = Generatekmodes(m, s0, Ltot, L);
iNmax = GenerateNmodes(h, s0, iNtot);

cout<< "the size of m is"<< m.size()<<"\n";

cout<< " VSR " << Vsr<< "\n";

for(unsigned int i2 =1; i2 < V; i2++ ){
for(unsigned int i1 =0; i1 < i2;i1++){

	Vsr += Vsr12(i1,i2,h,iNmax,coord,s0);

}
}
    
//cout<< "Vlr12 = "<<  Vlr12(1,2,m,iLmax,coord,s0)<< "\n" ;
//cout<< "zero mode is " << (Z/2)*zeromode(m,iLmax,s0)<<"\n";
//cout<< "ilmax = "<<  iLmax <<" m[iLmax].x= "<< m[iLmax].x<< endl;

for(unsigned int i2 =1; i2 < V; i2++ ){
for(unsigned int i1 =0; i1 < i2;i1++){

	Vlr +=  Vlr12(i1,i2,m,iLmax,coord,s0);

}
}

cout<< "\n Vlr +Vsr " << Vsr12(0,1,h,iNmax,coord,s0) +  Vlr12(0,1,m,iLmax,coord,s0) << "\n";

//vector<double> s(3);
//double less,more,total_energy;
//
//if(y == 1){
//
//total_energy = 0;
//int i1 = 0;
//cout<< " s0 = " << s0 << "\n";
//for(unsigned int i2 =1; i2 < V; i2++ ){
//	en12 = 0;
//	s[0] = L*(coord[i1].x - coord[i2].x);
//	s[1] = L*(coord[i1].y - coord[i2].y);
//	s[2] = L*(coord[i1].z - coord[i2].z);
//
//	en12 =  Vlr12(i1,i2,m,iLmax,coord,s0) + Vsr12(i1,i2,h,iNmax,coord,s0);
//	total_energy += en12;
//	zx << abs(s[0])+ L*abs(s[1])+ L*L*abs(s[2]) << "\t" << en12 << "\n"; //<< "\t"<< less<< "\t"  << more <<"\n";
//
//}
//}


//FormFactor(m, iLmax,coord,V,Sk);
//VLr = VLongRange(m,iLmax, s0, Sk);
VSE = VSelfEnergy(s0,Z);
VEP = VelectronProton(s0,Z);
kzeromode = (Z/2)*zeromode(m,iLmax,s0);
Vlr += kzeromode;
    
//cout<<"\n"<< y <<" total_energy = " <<  sc*(total_energy +VSE + VEP + kzeromode)  << "\n";

madelung = Vlr + VSE + VEP + Vsr;  //VLr.real()

cout << "U0" << sc(kzeromode+VSE+VEP);

cout<<  "Lr= " << Vlr << " VSE = " << VSE<< " VEP:  " << VEP<< " Vsr:  " << Vsr<< "  Vsr+VEP+VSE  : "<<VEP+VSE+Vsr <<"\n";

cout <<"V : "<< sc*madelung<< "\n";

ewalde<< s0 << "\t" << sc*Vlr << "\t" << sc*(VEP+VSE+Vsr) <<"\t" <<sc*madelung <<"\t" <<endl;

}
    
//q.close();
//w.close();
//x.close();
zx.close();
ewalde.close();

return(0);

}


int Generatekmodes(vector<xyz> & m, double a0 ,int  Ltot,double L){

vector<double> k(3);
int kmax = (int) (5*L/(a0*M_PI));
int Lmax = 0;
cout<< "kmax = " << kmax <<"\n";
int ksquared;
for(int kx = -kmax+1; kx< kmax; kx++){
for(int ky = -kmax+1; ky< kmax; ky++){
for(int kz = -kmax+1; kz< kmax; kz++){

    ksquared = kx*kx + ky*ky + kz*kz;
	if(ksquared == 0 ||  ksquared > kmax*kmax)
        continue;
        
        if(Lmax <Ltot){

		m[Lmax].x = kx;
		m[Lmax].y = ky;
		m[Lmax].z = kz;
		Lmax++;
        }
}
}
}

cout<< "Lmax = " << Lmax << "m.x = "<< m[2].x  << endl;

	return Lmax;
}


int GenerateNmodes(vector<xyz> & h, double a0 ,int  iNtot){

int nmax = (int) (5 * a0 + 1);
int iNmax = 0;

vector<double> n(3);

cout<< "nmax = " << nmax <<"\n";

for(int nx  = 0; nx < 5; nx++){
for(int ny  = 0; ny < 5; ny++){
for(int nz  = 0; nz < 5; nz++){

    if( nx*nx + ny*ny + nz*nz < pow(10,2)&& iNmax< iNtot) {

		h[iNmax].x =  nx;
		h[iNmax].y =  ny;
		h[iNmax].z =  nz;
		iNmax++;
}

}
}
}

cout<< "Nmax = " << iNmax << "\n";

	return iNmax;
}


double Vsr12(int i1, int i2, vector<xyz> & h, int Nmax,  vector<xyz>& s, double s0){

	double d12;

	int nx, ny, nz, zeros = 0;

//	cout<< "sx1: " << s[0][i1]<< " sx2 "<< s[0][i2] << " sy1 "<<s[1][i1]<< " sy2 "<< s[1][i2]<<" or "<< sy2  <<"\n";
	double Vs12 = 0;

/*	for(int i  = 0; i < Nmax; i++){
	cout<< "h[0][iNmax] = " << h[0][i]<< "h[1][iNmax] = " << h[1][i]  << "\n"  ;
}  */

	vector<double> d(3);
	//cout<< "Nmax =" << Nmax<< "\n";
	for(int nc = 0; nc< Nmax; nc++  ){

        d[0] = s[i1].x - s[i2].x + h[nc].x - 2;
        d[1] = s[i1].y - s[i2].y + h[nc].y - 2;
        d[2] = s[i1].z - s[i2].z + h[nc].z - 2;
		d12 = sqrt( d[0]*d[0]+ d[1]*d[1] + d[2]*d[2] );

		if(d12 == 0){
			zeros++;
		cout << " hnc ="<< nc  <<" d0 = "<< d[0] <<" d1 = "<< d[1] <<" d2 = "<< d[2] << " sx1 " << s[i1].x <<" sx2 "<<  s[i2].x << " sy1 "<< s[i1].y <<   " sy2 "<< s[i2].y  << "\n";
}
		if(d12 != 0){
		Vs12 +=    erfc(d12/s0)/d12;

}
}

	return Vs12;
}

double Vlr12(int i1, int i2, vector<xyz>& m, int Lmax,  vector<xyz>& s, double s0){

    double dl2,ld, angle,cosi;
    double pi = M_PI;

    int zeros = 0;

    //	cout<< "sx1: " << s[0][i1]<< " sx2 "<< s[0][i2] << " sy1 "<<s[1][i1]<< " sy2 "<< s[1][i2]<<" or "<< sy2  <<"\n";
    dcomp Vl12 = 0;
    dcomp im = dcomp(0,1);

    /*	for(int i  = 0; i < Nmax; i++){
     cout<< "h[0][iNmax] = " << h[0][i]<< "h[1][iNmax] = " << h[1][i]  << "\n"  ;
     }  */
    vector<double> d(3),l(3);

    //cout<< "Nmax =" << Nmax<< "\n";
    for(int nc = 0; nc< Lmax; nc++  ){
        

        d[0] = s[i1].x - s[i2].x ;
        d[1] = s[i1].y - s[i2].y ;
        d[2] = s[i1].z - s[i2].z ;


        dl2 = m[nc].x * m[nc].x  + m[nc].y * m[nc].y + m[nc].z * m[nc].z ;
        

        angle = 2*pi*(m[nc].x*d[0] + m[nc].y*d[1] + m[nc].z*d[2]);
        cosi =   cos(angle);
        
        
//        if(dl2 == 0){
//            zeros++;
//            cout  <<" mx and nc"<< m[nc].x<<"  "<< nc  << " my "<<  m[nc].y  << endl;
//        }
//
//        if(dl2!= 0){
        Vl12 +=   exp(-pi*pi*s0*s0*dl2) *cos(angle)*(1/(pi*dl2));
         
//
//        }
    }

    return Vl12.real();
}


dcomp VLongRange(vector<xyz> & m, int Lmax, double s0, vector<dcomp> & Sk){
	vector<double> l(3);
    double pi = M_PI;
	//cout<< "pi ="<< pi<<"\n";
	dcomp Vlr = 0;
	dcomp Sk2;
	double dl2;
	double Expo;
	for( int il =0; il < Lmax; il++){
		l[0] = m[il].x;
		l[1] = m[il].y;
		l[2] = m[il].z;
		dl2 = dotproduct(l,l);
		Expo = exp(-pow(pi,2)*pow(s0,2)*dl2);
		//Sk2 = Sk[il]*conj(Sk[il]);
		Sk2 = norm(Sk[il]);
		Vlr += Sk2*(Expo/dl2 );

	}
	Vlr = Vlr/(2*pi);
	return Vlr;

}

void FormFactor(vector<xyz> & m, int Lmax, vector<xyz> & coord,int z, vector<dcomp> & Sk){
	vector<double> l(3) , s(3);
	dcomp expo;
	dcomp im;
	im = dcomp(0,1);
	double pi = M_PI;
    double cosi, angle, sini ;
	for(unsigned int i= 0; i< Lmax; i++){

		Sk[i] = 0;
	for(unsigned int j= 0; j< z; j++){
        angle = 2*pi*(m[i].x*coord[j].x + m[i].y*coord[j].x + m[i].z*coord[j].z);
        cosi =   cos(angle);
        //sini =   sin(angle);
        Sk[i] += cosi ;
}
}
}


double VSelfEnergy(double s0, double Z){

	double Vcr=0 ;
	double pi = 2 * asin(1);


	Vcr -= Z/(s0*sqrt(pi));
	cout << "Vcr"<<Vcr<< "\n";
	return Vcr;
}


double VelectronProton(double s0, double Z){

	double Vep=0;
	double pi = 2 * asin(1);
	Vep = - Z*Z*M_PI*pow(s0,2)/2;
	cout << "Vep"<<Vep<< "\n";
	return Vep;
}


double dotproduct(vector<double> v1, vector<double> v2){
	double dott=0;
	for(unsigned int i=0; i<3; i++){
		dott += v1[i]*v2[i];
	}

	return dott;
	}

void generate_coord(vector<xyz> & v2, int z, int choice, double aLatt){

	int j = 0;

if(choice==0){

	for(unsigned int ix=0; ix<z; ix++){
	for(unsigned int iy=0; iy<z; iy++){
	for(unsigned int iz=0; iz<z; iz++){
		v2[j].x =  ix*aLatt;
		v2[j].y =  iy*aLatt;
		v2[j].z =  iz*aLatt;

		j += 1;
		
	}
	}
	}

}
    
//    for(unsigned int p = 0; p<V; p++ ){
//        cout<<p<<"::" <<v2[p].x<<":"<< v2[p].y<<":" << v2[p].z  <<"\n";
//    }

if(choice==1){
	for(unsigned int ix=0; ix<z; ix++){
	for(unsigned int iy=0; iy<z; iy++){
	for(unsigned int iz=0; iz<z; iz++){

		v2[j+0].x =  ix*aLatt;
		v2[j+0].y =  iy*aLatt;
		v2[j+0].z =  iz*aLatt;

		v2[j+1].x  =  (ix+ 0)*aLatt;
		v2[j+1].y  =  (iy+.5)*aLatt;
		v2[j+1].z  =  (iz+.5)*aLatt;

		v2[j+2].x  =  (ix+.5)*aLatt;
		v2[j+2].y  =  (iy+ 0)*aLatt;
		v2[j+2].z  =  (iz+.5)*aLatt;

		v2[j+3].x  =  (ix+.5)*aLatt;
		v2[j+2].y  =  (iy+.5)*aLatt;
		v2[j+2].z  =  (iz+ 0)*aLatt;



		j += 4;
		//cout<< i <<"\t"<<ii<<"\t"<< iii<<"\n";
	}
	}
	}
}

	cout<< "j is:  "<< j<<"\n";
}


double zeromode(vector<xyz> & m, int Lmax, double s0){

    double dl2, zeros = 0;
    double pi = 2 * asin(1);
    double Vl12 = 0;
//    vector<double> d(3),l(3);
    //cout<< "Nmax =" << Nmax<< "\n";
    for(unsigned int nc = 0; nc< Lmax; nc++  ){


        dl2 = m[nc].x * m[nc].x  + m[nc].y * m[nc].y + m[nc].z * m[nc].z ;
        //	cout << "d0 ="<<d[0]<<"d1 "<<d[1]<<"d2"<<d[2]  << "\n";

        if(dl2 == 0){
            zeros++;
        }

        if(dl2!= 0){
            Vl12 +=   exp(-pow(pi,2)*pow(s0,2)*dl2)*(1/(pi*dl2));
        }
    }

    return Vl12;
}
