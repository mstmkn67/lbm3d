#include "LBM3dPlates.h"

LBM3dPlates::LBM3dPlates(
	LBM3dLattice* lattice,
	double viscosity,double density,double dt,
	double delta_p,
	double u_down,double v_down,double w_down,
	double u_up,double v_up,double w_up)
	:LBM3dSimulator(lattice,viscosity,density,dt,delta_p),
	 u0(u_down),v0(v_down),w0(w_down),u2(u_up),v2(v_up),w2(w_up){
}

LBM3dPlates::~LBM3dPlates(){}

void LBM3dPlates::collide(){
	for(int k=0;k<=lattice->nz;k++){
		for(int j=1;j<lattice->ny;j++){
			for(int i=0;i<=lattice->nx;i++){
				for(int n=0;n<15;n++){
					calc_collide(n,i,j,k);
				}
			}
		}
	}
	for(int k=0;k<=lattice->nz;k++){
		for(int i=0;i<=lattice->nx;i++){
			collide_down_wall(i,k);//down
			collide_up_wall(i,k);//up
		}
	}
}

void LBM3dPlates::collide_down_wall(int i,int k){
	int j=0;
	//0,1,3,4,5,6,9,11,12,14
	calc_collide(0,i,j,k);
	calc_collide(1,i,j,k);
	calc_collide(3,i,j,k);
	calc_collide(4,i,j,k);
	calc_collide(5,i,j,k);
	calc_collide(6,i,j,k);
	calc_collide(9,i,j,k);
	calc_collide(11,i,j,k);
	calc_collide(12,i,j,k);
	calc_collide(14,i,j,k);
	//2,7,8,10,13
	nonslip_down_wall(i,k);
}

void LBM3dPlates::collide_up_wall(int i,int k){
	int j=lattice->ny;
	//0,1,2,3,4,6,7,8,10,13
	calc_collide(0,i,j,k);
	calc_collide(1,i,j,k);
	calc_collide(2,i,j,k);
	calc_collide(3,i,j,k);
	calc_collide(4,i,j,k);
	calc_collide(6,i,j,k);
	calc_collide(7,i,j,k);
	calc_collide(8,i,j,k);
	calc_collide(10,i,j,k);
	calc_collide(13,i,j,k);
	//5,9,11,12,14
	nonslip_up_wall(i,k);
}

void LBM3dPlates::nonslip_down_wall(int i,int k){
	//2,7,8,10,13
	int j=0;
	double r=lattice->f[0][i][j][k]+lattice->f[1][i][j][k]+lattice->f[3][i][j][k]
	        +lattice->f[4][i][j][k]+lattice->f[6][i][j][k]
	        +2.*(lattice->f[5][i][j][k]+lattice->f[9][i][j][k]+lattice->f[11][i][j][k]
	            +lattice->f[12][i][j][k]+lattice->f[14][i][j][k]);
	r/=1.-v0/c;
	lattice->f[2][i][j][k]=lattice->f[5][i][j][k]+2./3.*r*v0/c;
	lattice->f[7][i][j][k]=lattice->f[11][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k])
	                       +0.25*r*(u0+w0)/c+r*v0/c/12.;
	lattice->f[8][i][j][k]=lattice->f[12][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k])
	                       +0.25*r*(-u0+w0)/c+r*v0/c/12.;
	lattice->f[10][i][j][k]=lattice->f[14][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k])
	                       +0.25*r*(u0-w0)/c+r*v0/c/12.;
	lattice->f[13][i][j][k]=lattice->f[9][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k])
	                       +0.25*r*(-u0-w0)/c+r*v0/c/12.;
}

void LBM3dPlates::nonslip_up_wall(int i,int k){
	//5,9,11,12,14
	int j=lattice->ny;
	double r=lattice->f[0][i][j][k]+lattice->f[1][i][j][k]+lattice->f[3][i][j][k]
	        +lattice->f[4][i][j][k]+lattice->f[6][i][j][k]
	        +2.*(lattice->f[2][i][j][k]+lattice->f[13][i][j][k]+lattice->f[7][i][j][k]
	            +lattice->f[8][i][j][k]+lattice->f[10][i][j][k]);
	r/=1.+v2/c;
	lattice->f[5][i][j][k]=lattice->f[2][i][j][k]-2./3.*r*v2/c;
	lattice->f[9][i][j][k]=lattice->f[13][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k])
	                       +0.25*r*(u2+w2)/c-r*v2/c/12.;
	lattice->f[11][i][j][k]=lattice->f[7][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k])
	                       +0.25*r*(-u2-w2)/c-r*v2/c/12.;
	lattice->f[12][i][j][k]=lattice->f[8][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k])
	                       +0.25*r*(u2-w2)/c-r*v2/c/12.;
	lattice->f[14][i][j][k]=lattice->f[10][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k])
	                       +0.25*r*(-u2+w2)/c-r*v2/c/12.;
}

void LBM3dPlates::set_boundary_force(){
	/*
	for(int j=0;j<=lattice->ny;j++){
		lattice->bfx[lattice->nx][j]=lattice->bfx[0][j];
		lattice->bfy[lattice->nx][j]=lattice->bfy[0][j];
	}*/
}

void LBM3dPlates::set_boundary_viscosity(){
	/*
	for(int j=0;j<=lattice->ny;j++){
		lattice->I[lattice->nx][j]=lattice->I[0][j];
		lattice->vis[lattice->nx][j]=lattice->vis[0][j];
	}*/
}
