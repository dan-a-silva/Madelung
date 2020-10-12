//#include "montecarlo.hpp"

#include <ctime>
#include <vector>
#include <cmath>
#include <fstream>
#include <complex>
#include <string>
#include <iostream>
#include <chrono>


using namespace std;

typedef vector<int> Row;
typedef vector<Row> Matrix;
typedef vector<double> Roww;
typedef vector<Roww> Matrixx;
typedef vector<Matrix> Matrix3d;
typedef complex<double> dcomp;
typedef struct{
	int x,y,z;
} xyz;
double dotproduct(vector<double> v1, vector<double> v2);
double Vsr12(int i1, int i2, vector<xyz> & h, int Nmax,  Matrixx & s, double s0);
void FormFactor(vector<xyz> & m, int Lmax,Matrixx & coord,int z,vector<dcomp> & Sk);
double VSelfEnergy(double s0,double Z);
dcomp VLongRange(vector<xyz> & m, int Lmax, double s0, vector<dcomp> & Sk);
double VelectronProton(double s0,double Z);
int GenerateNmodes(vector<xyz> & h, double a0, int iNtot);
int Generatekmodes(vector<xyz> & m, double a0 ,int  Ltot,double L);
void generate_coord(Matrixx& v2, int z, int choice, double aLatt);




//void generate_coord(Matrixx & v2,Matrixx & prim, int z, int choice);


int main(void){

dcomp VLr ;
int iLmax, iNmax;
int Ltot = 50000, iNtot = 50000;
ofstream q,w,x,zx ;
vector<dcomp> Sk(Ltot);
double L =  2;
double Z =  4*pow(L,3);
int z   = 2;
int V   = 4*pow(z,3);
int fcc = 1;

double aEff = L;
double aLatt = 1/L;
double madelung,VSE,VEP,Vsr;
double a0 = .225;
double da = .002;
double s0;

Matrixx coord;
coord = Matrixx(3,Roww(V));
vector<xyz> m(Ltot);
vector<xyz> h(iNtot);
// m = Matrixx(3,Roww(Ltot));
//h = Matrixx(3,Roww(iNtot));

generate_coord(coord,z,fcc,aLatt);


if (fcc == 1){
aEff = 1/(pow(4,.333333)*L);
}

double sc = 2*aEff*(1/Z);
cout <<"sc" <<sc << "\n ";
cout<< "V" << V << "\t"<< "1/L = \t"<< 1/L <<"\t aEff:" << aEff<< "\n";

string ss = to_string(z);

q.open("long_range"+ss+"fcc"+".txt");
w.open("short_range"+ss+"fcc"+".txt");
x.open("madelung"+ss+"fcc"+".txt");
zx.open("energy_distance"+ss+"fcc"+".txt");


// for(unsigned int p = 0; p<V; p++ ){
//
// cout<<p<<"::" <<L*coord[0][p]<<":"<< L*coord[1][p]<<":" << L*coord[2][p]  <<"\n";
// }

std::chrono::time_point<std::chrono::system_clock> start, end;
std::chrono::duration<double> elapse;
for(unsigned int y = 0; y < 150; y++){
start = std::chrono::system_clock::now();
a0 += da;

s0 = a0/L;
cout<< "s0 = " << s0 << "\n";

Vsr = 0;

iLmax = Generatekmodes(m, s0, Ltot, L);
end = std::chrono::system_clock::now();
elapse = end- start;
cout<< "elapsed =" << elapse.count()<<endl;
start = std::chrono::system_clock::now();

iNmax = GenerateNmodes(h, s0, iNtot);
end = std::chrono::system_clock::now();
elapse = end- start;
cout<< "elapsed =" << elapse.count()<<endl;

cout<< "the size of m is"<< m.size()<<"\n";


/*for(int i  = 0; i < iNmax; i++){
	cout<< i << " h[0][i] = " << h[0][i]<< "h[1][iNmax] = " << h[1][i] << "h[2][iNmax] = " << h[2][i]  << "\n"  ;
} */

cout<< " VSR " << Vsr<< "\n";

start = std::chrono::system_clock::now();

for(unsigned int i2 =1; i2 < V; i2++ ){
for(unsigned int i1 =0; i1 < i2;i1++){

	Vsr += Vsr12(i1,i2,h,iNmax,coord,s0);

}
}
    
end = std::chrono::system_clock::now();
elapse = end- start;
cout<< "elapsed =" << elapse.count()<<endl;

start = std::chrono::system_clock::now();
FormFactor(m, iLmax,coord,V,Sk);
end = std::chrono::system_clock::now();
elapse = end- start;
cout<< "elapsed =" << elapse.count()<<endl;
    
VLr = VLongRange(m,iLmax, s0, Sk);
//cout<< "hello \n";
VSE = VSelfEnergy(s0,Z);
VEP = VelectronProton(s0,Z);

madelung = VLr.real()+ VSE + VEP + Vsr;

cout<<  "Lr= " <<VLr.real()<< " VSE = " << VSE<< " VEP:  " << VEP<< " Vsr:  " << Vsr<< "  Vsr+VEP+VSE  : "<<VEP+VSE+Vsr <<"\n";

cout <<"V : "<< sc*madelung<< "\n";

q << s0 << "\t" << sc*VLr.real() <<"\n";
w << s0 << "\t" << sc*(VEP+VSE+Vsr) <<"\n";
x << s0 << "\t" << sc*madelung <<"\n";

}
q.close();
w.close();
x.close();
zx.close();
return(0);

}


int Generatekmodes(vector<xyz> & m, double a0 ,int  Ltot,double L){

int kmax = (int) (5*L/(a0*M_PI));
int Lmax = 0;
cout<< "kmax = " << kmax <<"\n";
int ksquared;
for(int kx = -kmax+1; kx< kmax; kx++){
for(int ky = -kmax+1 ; ky< kmax; ky++){
for(int kz = -kmax+1 ; kz< kmax; kz++){

	ksquared = kx*kx + ky*ky + kz*kz;
	if(ksquared == 0 || ksquared > kmax*kmax)
	continue;
	//cout<< "k dot d = " << dotproduct(k,k)<<"\n";

	if(  Lmax < Ltot){

		m[Lmax].x =  kx;
		m[Lmax].y =  ky;
		m[Lmax].z =  kz;
		Lmax++;

}
}
}
}

//cout<< "Lmax = " << Lmax << "\n";
	return Lmax;
}

int GenerateNmodes(vector<xyz> & h, double a0 ,int  iNtot){

int nmax = (int) (5 * a0 + 1);
int iNmax = 0;

vector<double> n(3);

cout<< "nmax = " << nmax <<"\n";

for(int nx  = 0; nx < 5; nx++){
for(int ny  = 0; ny < 5; ny++){
for(int nz =  0;  nz <5; nz++){
	//cout<< "the answer" <<  t << "\n";
	//cout<< "n dot n = " << dotproduct(n,n)<<"\n";
	if( nx*nx + ny*ny + nz*nz < pow(10,2)&& iNmax< iNtot) {

		h[iNmax].x = nx;
		h[iNmax].y = ny;
		h[iNmax].z = nz;
		iNmax++;
}

}
}
}

cout<< "Nmax = " << iNmax << "\n";

	return iNmax;
}


double Vsr12(int i1, int i2, vector<xyz> & h, int Nmax,  Matrixx & s, double s0){

	double sx1,sy1,sz1, sx2,sy2,sz2,d12;

	sx1 = s[0][i1];
	sy1 = s[1][i1];
	sz1 = s[2][i1];

	sx2 = s[0][i2];
	sy2 = s[1][i2];
	sz2 = s[2][i2];

	int nx, ny, nz, zeros = 0;

//	cout<< "sx1: " << s[0][i1]<< " sx2 "<< s[0][i2] << " sy1 "<<s[1][i1]<< " sy2 "<< s[1][i2]<<" or "<< sy2  <<"\n";
	double Vs12 = 0;

	vector<double> d(3);
	//cout<< "Nmax =" << Nmax<< "\n";
	for(unsigned int nc = 0; nc< Nmax; nc++  ){

		// nx = h[0][nc];
		// ny = h[1][nc];
		// nz = h[2][nc];

		d[0] = sx1 - sx2 + h[nc].x -2;
		d[1] = sy1 - sy2 + h[nc].y -2;
		d[2] = sz1 - sz2 + h[nc].z -2;
		d12 = sqrt(dotproduct(d,d));
	//	cout << "d0 ="<<d[0]<<"d1 "<<d[1]<<"d2"<<d[2]  << "\n";

		if(d12 == 0){
			zeros++;
		cout << " nc ="<< nc  <<" d0 = "<<d[0]<<" d1 = "<<d[1]<<" d2 = "<<d[2] << " sx1 " << sx1 <<" sx2 "<<sx2<< " sy1 "<<sy1<<   " sy2 "<< sy2  << "\n";
}

		if(d12 != 0){
		Vs12 +=    erfc(d12/s0)/d12;

}
}



	return Vs12;
}


dcomp VLongRange(vector<xyz>& m, int Lmax, double s0, vector<dcomp> & Sk){

	double pi = M_PI;
	//cout<< "pi ="<< pi<<"\n";
	dcomp Vlr = 0;
	dcomp Sk2;
	double dl2;
	double Expo;
	for( int il =0; il < Lmax; il++){
		dl2 = m[il].x * m[il].x + m[il].y *m[il].y + m[il].z *m[il].z;
		Expo = exp(-pi*pi*s0*s0*dl2);
		//Sk2 = Sk[il]*conj(Sk[il]);
		Sk2 = norm(Sk[il]);
		Vlr += Sk2*(Expo/dl2 );

	}
	Vlr = Vlr/(2*pi);
	return Vlr;

}

void FormFactor(vector<xyz>& m, int Lmax,Matrixx & coord,int z, vector<dcomp> & Sk){
	vector<double> l(3) , s(3);
	dcomp expo;
	dcomp im;
	im = dcomp(0,1);
	double pi = M_PI;
	double cosi, angle, sini ;


	for(unsigned int i= 0; i< Lmax; i++){

		Sk[i] = 0;
	for(unsigned int j= 0; j< z; j++){

		//expo = exp(-2*pi*im*(m[0][i]*coord[0][j]+m[1][i]*coord[1][j]+m[2][i]*coord[2][j]));
		angle = 2*pi*(m[i].x*coord[0][j]+m[i].y*coord[1][j]+m[i].z*coord[2][j]);
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


double VelectronProton(double s0,double Z){

	double Vep=0;
	double pi = 2 * asin(1);
	Vep = - pow(Z,2)*pi*pow(s0,2)/2;
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

void generate_coord(Matrixx& v2, int z, int choice, double aLatt){

	int j = 0;

if(choice==0){

	for(unsigned int ix=0; ix<z; ix++){
	for(unsigned int iy=0; iy<z; iy++){
	for(unsigned int iz=0; iz<z; iz++){
		v2[0][j] =  ix*aLatt;
		v2[1][j] =  iy*aLatt;
		v2[2][j] =  iz*aLatt;

		j += 1;
		//cout<< i <<"\t"<<ii<<"\t"<< iii<<"\n";
	}
	}
	}

}

if(choice==1){
	for(unsigned int ix=0; ix<z; ix++){
	for(unsigned int iy=0; iy<z; iy++){
	for(unsigned int iz=0; iz<z; iz++){

		v2[0][j+0] =  ix*aLatt;
		v2[1][j+0] =  iy*aLatt;
		v2[2][j+0] =  iz*aLatt;

		v2[0][j+1]  =  (ix+ 0)*aLatt;
		v2[1][j+1]  =  (iy+.5)*aLatt;
		v2[2][j+1]  =  (iz+.5)*aLatt;

		v2[0][j+2]  =  (ix+.5)*aLatt;
		v2[1][j+2]  =  (iy+ 0)*aLatt;
		v2[2][j+2]  =  (iz+.5)*aLatt;

		v2[0][j+3]  =  (ix+.5)*aLatt;
		v2[1][j+3]  =  (iy+.5)*aLatt;
		v2[2][j+3]  =  (iz+ 0)*aLatt;



		j += 4;
		//cout<< i <<"\t"<<ii<<"\t"<< iii<<"\n";
	}
	}
	}
}

	cout<< "j is:  "<< j<<"\n";
}
