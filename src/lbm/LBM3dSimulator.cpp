#include "LBM3dSimulator.h"

#define PI 3.1415926535

LBM3dSimulator::LBM3dSimulator(LBM3dLattice* lat,double vis,double den,double t,double dp)
:lattice(lat),delp(dp),viscosity(vis),density(den),dt(t){
	c=lattice->dx/dt;
	cs2=c*c/3.0;
	cx[ 0]= 0;cy[ 0]= 0;cz[ 0]= 0;w[ 0]=2.0/9.0;
	cx[ 1]= 1;cy[ 1]= 0;cz[ 1]= 0;w[ 1]=1.0/9.0;
	cx[ 2]= 0;cy[ 2]= 1;cz[ 2]= 0;w[ 2]=1.0/9.0;
	cx[ 3]= 0;cy[ 3]= 0;cz[ 3]= 1;w[ 3]=1.0/9.0;
	cx[ 4]=-1;cy[ 4]= 0;cz[ 4]= 0;w[ 4]=1.0/9.0;
	cx[ 5]= 0;cy[ 5]=-1;cz[ 5]= 0;w[ 5]=1.0/9.0;
	cx[ 6]= 0;cy[ 6]= 0;cz[ 6]=-1;w[ 6]=1.0/9.0;
	cx[ 7]= 1;cy[ 7]= 1;cz[ 7]= 1;w[ 7]=1.0/72.0;
	cx[ 8]=-1;cy[ 8]= 1;cz[ 8]= 1;w[ 8]=1.0/72.0;
	cx[ 9]= 1;cy[ 9]=-1;cz[ 9]= 1;w[ 9]=1.0/72.0;
	cx[10]= 1;cy[10]= 1;cz[10]=-1;w[10]=1.0/72.0;
	cx[11]=-1;cy[11]=-1;cz[11]=-1;w[11]=1.0/72.0;
	cx[12]= 1;cy[12]=-1;cz[12]=-1;w[12]=1.0/72.0;
	cx[13]=-1;cy[13]= 1;cz[13]=-1;w[13]=1.0/72.0;
	cx[14]=-1;cy[14]=-1;cz[14]= 1;w[14]=1.0/72.0;
}

LBM3dSimulator::~LBM3dSimulator(){
}

void LBM3dSimulator::initial(){
	for(int k=0;k<=lattice->nz;k++){
		for(int j=0;j<=lattice->ny;j++){
			for(int i=0;i<=lattice->nx;i++){
				for(int n=0;n<15;n++){
					lattice->f[n][i][j][k]=density*w[n];
					lattice->f0[n][i][j][k]=density*w[n];
				}
				lattice->vis[i][j][k]=viscosity;
			}
		}
	}
	calc_density_velocity();
}

void LBM3dSimulator::update(){
	set_boundary_force();
	set_boundary_viscosity();
	collide();
	calc_density_velocity();
	time_update();
}

void LBM3dSimulator::reset_bforce(){
	for(int k=0;k<=lattice->nz;k++){
		for(int j=0;j<=lattice->ny;j++){
			for(int i=0;i<=lattice->nx;i++){
				lattice->bfx[i][j][k]=delp/lattice->lx;
				lattice->bfy[i][j][k]=0.0;
				lattice->bfz[i][j][k]=0.0;
			}
		}
	}
}

//void LBM2dSimulator::reset_body_force(){
//}

void LBM3dSimulator::calc_density_velocity(){
	for(int k=0;k<=lattice->nz;k++){
		for(int j=0;j<=lattice->ny;j++){
			for(int i=0;i<=lattice->nx;i++){
				lattice->rho[i][j][k]=lattice->f[0][i][j][k];
				for(int n=1;n<15;n++){
					lattice->rho[i][j][k]+=lattice->f[n][i][j][k];
				}
			}
		}
	}
	for(int k=0;k<=lattice->nz;k++){
		for(int j=0;j<=lattice->ny;j++){
			for(int i=0;i<=lattice->nx;i++){
				lattice->u[i][j][k]=lattice->f[0][i][j][k]*(cx[0]*c);
				lattice->v[i][j][k]=lattice->f[0][i][j][k]*(cy[0]*c);
				lattice->w[i][j][k]=lattice->f[0][i][j][k]*(cz[0]*c);
				for(int n=1;n<15;n++){
					lattice->u[i][j][k]+=lattice->f[n][i][j][k]*(cx[n]*c);
					lattice->v[i][j][k]+=lattice->f[n][i][j][k]*(cy[n]*c);
					lattice->w[i][j][k]+=lattice->f[n][i][j][k]*(cz[n]*c);
				}
				lattice->u[i][j][k]/=lattice->rho[i][j][k];
				lattice->v[i][j][k]/=lattice->rho[i][j][k];
				lattice->w[i][j][k]/=lattice->rho[i][j][k];
			}
		}
	}
}

void LBM3dSimulator::calc_collide(int n,int i,int j,int k){
	int ii=i-cx[n],jj=j-cy[n],kk=k-cz[n];
	if(ii==-1)ii=lattice->nx;else if(ii==lattice->nx+1)ii=0;
	if(jj==-1)jj=lattice->ny;else if(jj==lattice->ny+1)jj=0;
	if(kk==-1)kk=lattice->nz;else if(kk==lattice->nz+1)kk=0;
	double tau=lattice->vis[ii][jj][kk]/lattice->rho[ii][jj][kk]/cs2/dt+0.5;
	double ueq=lattice->u[ii][jj][kk]+tau*lattice->bfx[ii][jj][kk]*dt/lattice->rho[ii][jj][kk];
	double veq=lattice->v[ii][jj][kk]+tau*lattice->bfy[ii][jj][kk]*dt/lattice->rho[ii][jj][kk];
	double weq=lattice->w[ii][jj][kk]+tau*lattice->bfz[ii][jj][kk]*dt/lattice->rho[ii][jj][kk];
	double uc=ueq*(cx[n]*c)+veq*(cy[n]*c)+weq*(cz[n]*c);
	double uu=ueq*ueq+veq*veq+weq*weq;
	double feq=lattice->rho[ii][jj][kk]*w[n]*(1.0+uc/cs2+0.5*uc*uc/cs2/cs2-0.5*uu/cs2);
	double a=lattice->f0[n][ii][jj][kk]-(lattice->f0[n][ii][jj][kk]-feq)/tau;
	//a+=w[n]*dt/cs2*(lattice->bfx[ii][jj]*cx[n]*c+lattice->bfy[ii][jj]*cy[n]*c);
	lattice->f[n][i][j][k]=a;
}

void LBM3dSimulator::time_update(){
	for(int k=0;k<=lattice->nz;k++){
		for(int j=0;j<=lattice->ny;j++){
			for(int i=0;i<=lattice->nx;i++){
				for(int n=0;n<15;n++){
					lattice->f0[n][i][j][k]=lattice->f[n][i][j][k];
				}
			}
		}
	}
}

void LBM3dSimulator::set_bulk_viscosity(){
}
	
#undef PI
