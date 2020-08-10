#ifndef _LBM_3D_LATTICE_H_
#define _LBM_3D_LATTICE_H_

#include "../SimpleArray.h"
#include <iostream>
using namespace std;

class LBM3dLattice{
public:
	LBM3dLattice(int nx,int ny,int nz,double dx);
	virtual ~LBM3dLattice();
	int nx,ny,nz;
	double lx,ly,lz;
	double dx;
	Array3d<double> u,v,w;
	Array3d<double> f[15];//distribution function
	Array3d<double> f0[15];//distribution function previous time
	Array3d<double> bfx,bfy,bfz;//body force
	Array3d<double> rho;//density
	Array3d<double> vis;//viscosity field
	Array3d<double> I;//maker field
private:
};

#endif // _LBM_3D_LATTICE_H_
