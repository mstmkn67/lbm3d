#include "LBM3dRect.h"

LBM3dRect::LBM3dRect(
	LBM3dLattice* lattice,
	double viscosity,double density,double dt,
	double delta_p)
	:LBM3dSimulator(lattice,viscosity,density,dt,delta_p){
}

LBM3dRect::~LBM3dRect(){}

void LBM3dRect::collide(){
	for(int k=1;k<lattice->nz;k++){
		for(int j=1;j<lattice->ny;j++){
			for(int i=0;i<=lattice->nx;i++){
				for(int n=0;n<15;n++){
					calc_collide(n,i,j,k);
				}
			}
		}
	}
	for(int k=1;k<lattice->nz;k++){
		for(int i=0;i<=lattice->nx;i++){
			collide_down_wall(i,k);//down
			collide_up_wall(i,k);//up
		}
	}
	for(int j=1;j<lattice->ny;j++){
		for(int i=0;i<=lattice->nx;i++){
			collide_back_wall(i,j);
			collide_front_wall(i,j);
		}
	}
	for(int i=0;i<=lattice->nx;i++){
		collide_down_back_edge(i);
		collide_down_front_edge(i);
		collide_up_back_edge(i);
		collide_up_front_edge(i);
	}
}

void LBM3dRect::collide_down_wall(int i,int k){
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

void LBM3dRect::collide_up_wall(int i,int k){
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

void LBM3dRect::collide_back_wall(int i,int j){
	int k=0;
	//0,1,2,4,5,6,10,11,12,13
	calc_collide(0,i,j,k);
	calc_collide(1,i,j,k);
	calc_collide(2,i,j,k);
	calc_collide(4,i,j,k);
	calc_collide(5,i,j,k);
	calc_collide(6,i,j,k);
	calc_collide(10,i,j,k);
	calc_collide(11,i,j,k);
	calc_collide(12,i,j,k);
	calc_collide(13,i,j,k);
	//3,7,8,9,14
	nonslip_back_wall(i,j);
}

void LBM3dRect::collide_front_wall(int i,int j){
	int k=lattice->nz;
	//0,1,2,3,4,5,7,8,9,14
	calc_collide(0,i,j,k);
	calc_collide(1,i,j,k);
	calc_collide(2,i,j,k);
	calc_collide(3,i,j,k);
	calc_collide(4,i,j,k);
	calc_collide(5,i,j,k);
	calc_collide(7,i,j,k);
	calc_collide(8,i,j,k);
	calc_collide(9,i,j,k);
	calc_collide(14,i,j,k);
	//6,10,11,12,13
	nonslip_front_wall(i,j);
}

void LBM3dRect::nonslip_down_wall(int i,int k){
	//2,7,8,10,13
	int j=0;
	lattice->f[2][i][j][k]=lattice->f[5][i][j][k];
	lattice->f[7][i][j][k]=lattice->f[11][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k]);
	lattice->f[8][i][j][k]=lattice->f[12][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k]);
	lattice->f[10][i][j][k]=lattice->f[14][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k]);
	lattice->f[13][i][j][k]=lattice->f[9][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k]);
}

void LBM3dRect::nonslip_up_wall(int i,int k){
	//5,9,11,12,14
	int j=lattice->ny;
	lattice->f[5][i][j][k]=lattice->f[2][i][j][k];
	lattice->f[9][i][j][k]=lattice->f[13][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k]);
	lattice->f[11][i][j][k]=lattice->f[7][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k]);
	lattice->f[12][i][j][k]=lattice->f[8][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	                              +lattice->f[3][i][j][k]-lattice->f[6][i][j][k]);
	lattice->f[14][i][j][k]=lattice->f[10][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]-lattice->f[4][i][j][k]
	                              -lattice->f[3][i][j][k]+lattice->f[6][i][j][k]);
}

void LBM3dRect::nonslip_back_wall(int i,int j){
	//3,7,8,9,14
	int k=0;
	lattice->f[3][i][j][k]=lattice->f[6][i][j][k];
	lattice->f[7][i][j][k]=lattice->f[11][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]-lattice->f[2][i][j][k]
	                              +lattice->f[4][i][j][k]+lattice->f[5][i][j][k]);
	lattice->f[8][i][j][k]=lattice->f[12][i][j][k]
	                       +0.25*(-lattice->f[2][i][j][k]-lattice->f[4][i][j][k]
	                              +lattice->f[1][i][j][k]+lattice->f[5][i][j][k]);
	lattice->f[9][i][j][k]=lattice->f[13][i][j][k]
	                       +0.25*(-lattice->f[1][i][j][k]-lattice->f[5][i][j][k]
	                              +lattice->f[2][i][j][k]+lattice->f[4][i][j][k]);
	lattice->f[14][i][j][k]=lattice->f[10][i][j][k]
	                       +0.25*(-lattice->f[4][i][j][k]-lattice->f[5][i][j][k]
	                              +lattice->f[1][i][j][k]+lattice->f[2][i][j][k]);
}

void LBM3dRect::nonslip_front_wall(int i,int j){
	//6,10,11,12,13
	int k=lattice->nz;
	lattice->f[6][i][j][k]=lattice->f[3][i][j][k];
	lattice->f[11][i][j][k]=lattice->f[7][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]+lattice->f[2][i][j][k]
	                              -lattice->f[4][i][j][k]-lattice->f[5][i][j][k]);
	lattice->f[12][i][j][k]=lattice->f[8][i][j][k]
	                       +0.25*( lattice->f[2][i][j][k]+lattice->f[4][i][j][k]
	                              -lattice->f[1][i][j][k]-lattice->f[5][i][j][k]);
	lattice->f[13][i][j][k]=lattice->f[9][i][j][k]
	                       +0.25*( lattice->f[1][i][j][k]+lattice->f[5][i][j][k]
	                              -lattice->f[2][i][j][k]-lattice->f[4][i][j][k]);
	lattice->f[10][i][j][k]=lattice->f[14][i][j][k]
	                       +0.25*( lattice->f[4][i][j][k]+lattice->f[5][i][j][k]
	                              -lattice->f[1][i][j][k]-lattice->f[2][i][j][k]);
}

void LBM3dRect::collide_down_back_edge(int i){
	int j=0,k=0;
	//0,1,4,5,6,11,12
	calc_collide(0,i,j,k); calc_collide(1,i,j,k);
	calc_collide(4,i,j,k); calc_collide(5,i,j,k);
	calc_collide(6,i,j,k); calc_collide(11,i,j,k);
	calc_collide(12,i,j,k);
	//2,3,7,8, 9,10,13,14
	nonslip_down_back_edge(i);
}

void LBM3dRect::collide_down_front_edge(int i){
	int j=0,k=lattice->nz;
	//0,1,3,4,5,9,14
	calc_collide(0,i,j,k); calc_collide(1,i,j,k);
	calc_collide(3,i,j,k); calc_collide(4,i,j,k);
	calc_collide(5,i,j,k); calc_collide(9,i,j,k);
	calc_collide(14,i,j,k);
	//2,6,10,13, 7,8,11,12
	nonslip_down_front_edge(i);
}

void LBM3dRect::collide_up_back_edge(int i){
	int j=lattice->ny,k=0;
	//0,1,2,4,6,10,13
	calc_collide(0,i,j,k); calc_collide(1,i,j,k);
	calc_collide(2,i,j,k); calc_collide(4,i,j,k);
	calc_collide(6,i,j,k); calc_collide(10,i,j,k);
	calc_collide(13,i,j,k);
	//3,5,9,14, 7,8,11,12
	nonslip_up_back_edge(i);
}

void LBM3dRect::collide_up_front_edge(int i){
	int j=lattice->ny,k=lattice->nz;
	//0,1,2,3,4,7,8
	calc_collide(0,i,j,k); calc_collide(1,i,j,k);
	calc_collide(2,i,j,k); calc_collide(3,i,j,k);
	calc_collide(4,i,j,k); calc_collide(7,i,j,k);
	calc_collide(8,i,j,k);
	//5,6,11,12, 9,10,13,14
	nonslip_up_front_edge(i);
}

void LBM3dRect::nonslip_down_back_edge(int i){
	int j=0,k=0;
	lattice->f[2][i][j][k]=lattice->f[5][i][j][k];
	lattice->f[3][i][j][k]=lattice->f[6][i][j][k];
	lattice->f[7][i][j][k]=lattice->f[11][i][j][k]
	                      -0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	lattice->f[8][i][j][k]=lattice->f[12][i][j][k]
	                      +0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	double r=18./17.*(
	         lattice->f[0][i][j][k]+lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	        +2.*(lattice->f[5][i][j][k]+lattice->f[6][i][j][k]
	            +lattice->f[11][i][j][k]+lattice->f[12][i][j][k])
	         )/72;
	lattice->f[9][i][j][k]=r; lattice->f[10][i][j][k]=r;
	lattice->f[13][i][j][k]=r;lattice->f[14][i][j][k]=r;
}

void LBM3dRect::nonslip_down_front_edge(int i){
	int j=0,k=lattice->nz;
	lattice->f[2][i][j][k]=lattice->f[5][i][j][k];
	lattice->f[6][i][j][k]=lattice->f[3][i][j][k];
	lattice->f[10][i][j][k]=lattice->f[14][i][j][k]
	                      -0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	lattice->f[13][i][j][k]=lattice->f[9][i][j][k]
	                      +0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	double r=18./17.*(
	         lattice->f[0][i][j][k]+lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	        +2.*(lattice->f[3][i][j][k]+lattice->f[5][i][j][k]
	            +lattice->f[9][i][j][k]+lattice->f[14][i][j][k])
	         )/72;
	lattice->f[7][i][j][k]=r; lattice->f[8][i][j][k]=r;
	lattice->f[11][i][j][k]=r;lattice->f[12][i][j][k]=r;
}

void LBM3dRect::nonslip_up_back_edge(int i){
	int j=lattice->ny,k=0;
	lattice->f[3][i][j][k]=lattice->f[6][i][j][k];
	lattice->f[5][i][j][k]=lattice->f[2][i][j][k];
	lattice->f[9][i][j][k]=lattice->f[13][i][j][k]
	                      -0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	lattice->f[14][i][j][k]=lattice->f[10][i][j][k]
	                      +0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	double r=18./17.*(
	         lattice->f[0][i][j][k]+lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	        +2.*(lattice->f[2][i][j][k]+lattice->f[6][i][j][k]
	            +lattice->f[10][i][j][k]+lattice->f[13][i][j][k])
	         )/72;
	lattice->f[7][i][j][k]=r; lattice->f[8][i][j][k]=r;
	lattice->f[11][i][j][k]=r;lattice->f[12][i][j][k]=r;
}

void LBM3dRect::nonslip_up_front_edge(int i){
	int j=lattice->ny,k=lattice->nz;
	lattice->f[5][i][j][k]=lattice->f[2][i][j][k];
	lattice->f[6][i][j][k]=lattice->f[3][i][j][k];
	lattice->f[12][i][j][k]=lattice->f[8][i][j][k]
	                      -0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	lattice->f[11][i][j][k]=lattice->f[7][i][j][k]
	                      +0.5*(lattice->f[1][i][j][k]-lattice->f[4][i][j][k]);
	double r=18./17.*(
	         lattice->f[0][i][j][k]+lattice->f[1][i][j][k]+lattice->f[4][i][j][k]
	        +2.*(lattice->f[2][i][j][k]+lattice->f[3][i][j][k]
	            +lattice->f[7][i][j][k]+lattice->f[8][i][j][k])
	         )/72;
	lattice->f[9][i][j][k]=r; lattice->f[10][i][j][k]=r;
	lattice->f[13][i][j][k]=r;lattice->f[14][i][j][k]=r;
}

void LBM3dRect::set_boundary_force(){
	/*
	for(int j=0;j<=lattice->ny;j++){
		lattice->bfx[lattice->nx][j]=lattice->bfx[0][j];
		lattice->bfy[lattice->nx][j]=lattice->bfy[0][j];
	}*/
}

void LBM3dRect::set_boundary_viscosity(){
	/*
	for(int j=0;j<=lattice->ny;j++){
		lattice->I[lattice->nx][j]=lattice->I[0][j];
		lattice->vis[lattice->nx][j]=lattice->vis[0][j];
	}*/
}
