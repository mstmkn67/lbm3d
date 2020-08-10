#include "LBM3dLattice.h"

LBM3dLattice::LBM3dLattice(int _nx,int _ny,int _nz,double _dx)
:nx(_nx),ny(_ny),nz(_nz),lx(_nx*_dx),ly(_ny*_dx),lz(_nz*dx),dx(_dx),
	u(0,nx,0,ny,0,nz),v(0,nx,0,ny,0,nz),w(0,nx,0,ny,0,nz),
	bfx(0,nx,0,ny,0,nz),bfy(0,nx,0,ny,0,nz),bfz(0,nx,0,ny,0,nz),
	rho(0,nx,0,ny,0,nz),vis(0,nx,0,ny,0,nz),I(0,nx,0,ny,0,nz)
{
	for(int i=0;i<15;i++){
		f[i].setBounds(0,nx,0,ny,0,nz);
		f0[i].setBounds(0,nx,0,ny,0,nz);
	}
}

LBM3dLattice::~LBM3dLattice(){
}
