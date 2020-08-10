#include "LBM3d.h"
//
#include "lbm/LBM3dPlates.h"
#include "lbm/LBM3dRect.h"
#include "lbm/LBM3dLeesEdwards.h"

LBM3d::LBM3d(UDFManager* _in,UDFManager* _out)
:in(_in),out(_out),iteration(0),
 fou(false),fop(false),forho(false),fof(false){
	timer=new Timer;
}

LBM3d::~LBM3d(){
	delete lbm;
	delete lattice;
	delete timer;
}

void LBM3d::update(){
	input_lbm();
	input_other();
	lbm->initial();
	lbm->set_bulk_viscosity();
	lbm->reset_bforce();
	output();
	//
	int total=in->i("simulation.time.simulation_steps");
	int report=in->i("simulation.time.record_steps");
	for(iteration=1;iteration<=total;iteration++){
		lbm->reset_bforce();
		lbm->update(); //
		if(iteration%report==0){
			cout << "report step at " << iteration << endl;
			output();
		}
	}
}

void LBM3d::input_lbm(){
	double dt=in->d("simulation.time.dt");
	double viscosity=in->d("fluid.viscosity");
	double density=in->d("fluid.mass_density");
	//
	int nx=in->i("system_size.nx");
	int ny=in->i("system_size.ny");
	int nz=in->i("system_size.nz");
	double dx=in->d("system_size.dx");
	lattice=new LBM3dLattice(nx,ny,nz,dx);
	if(in->s("simulation.boundary_condition.type")=="periodic"){
		lbm=new LBM3dLeesEdwards(lattice,viscosity,density,dt,
			                       &iteration,0.0);
	}else if(in->s("simulation.boundary_condition.type")=="Couette"){
		double vxu=in->d("simulation.boundary_condition.Couette.velocity_x_up");
		double vxd=in->d("simulation.boundary_condition.Couette.velocity_x_down");
		lbm=new LBM3dPlates(lattice,viscosity,density,dt,
		                    0.0, vxd,0.0,0.0, vxu,0.0,0.0);
	}else if(in->s("simulation.boundary_condition.type")=="PlatesPoiseuille"){
		double delp=in->d("simulation.boundary_condition.PlatesPoiseuille.delta_pressure");
		lbm=new LBM3dPlates(lattice,viscosity,density,dt,
		                    delp, 0.0,0.0,0.0, 0.0,0.0,0.0);
	}else if(in->s("simulation.boundary_condition.type")=="RectPoiseuille"){
		double delp=in->d("simulation.boundary_condition.RectPoiseuille.delta_pressure");
		lbm=new LBM3dRect(lattice,viscosity,density,dt,delp);
	}else if(in->s("simulation.boundary_condition.type")=="LeesEdwards"){
		double gd=in->d("simulation.boundary_condition.LeesEdwards.shear_rate");
		lbm=new LBM3dLeesEdwards(lattice,viscosity,density,dt,&iteration, gd);
	}
}

void LBM3d::input_other(){
	if(in->s("simulation.field_output.u")=="true")fou=true;
	if(in->s("simulation.field_output.p")=="true")fop=true;
	if(in->s("simulation.field_output.rho")=="true")forho=true;
	if(in->s("simulation.field_output.f")=="true")fof=true;
}

void LBM3d::output(){
	out->newRecord();
	for(int i=0;i<=lattice->nx;i++){
		for(int j=0;j<=lattice->ny;j++){
			for(int k=0;k<=lattice->nz;k++){
				if(fou){
					double u=lattice->u[i][j][k];
					nan_check(u,'u');//////////////////////////////////
					out->put("simulation_result.lattice.u[][][].u",u,INDEX(i,j,k));
					double v=lattice->v[i][j][k];
					nan_check(v,'v');/////////////////////////
					out->put("simulation_result.lattice.u[][][].v",v,INDEX(i,j,k));
					double w=lattice->w[i][j][k];
					nan_check(w,'w');/////////////////////////
					out->put("simulation_result.lattice.u[][][].w",w,INDEX(i,j,k));
				}
				if(fop){
					double gp=-lbm->delp/lattice->lx;
					double x=lattice->dx*i;
					nan_check(lattice->rho[i][j][k],'p');/////////////////////////
					out->put("simulation_result.lattice.p[][][]",lattice->rho[i][j][k]/3.0+gp*x,INDEX(i,j,k));
				}
				if(forho){
					out->put("simulation_result.lattice.rho[][][]",lattice->rho[i][j][k],INDEX(i,j,k));
				}
				if(fof){
					for(int n=0;n<15;n++){
						double f=lattice->f[n][i][j][k];
						nan_check(f,'f');//////////////////////////////////
						out->put("simulation_result.lattice.f[][][][]",f,INDEX(n,i,j,k));
					}
				}
			}
		}
	}
	//
	out->put("simulation_result.cpu_time",timer->get());
	out->write();
}

void LBM3d::nan_check(double x,char c){
  if( isnan(x) == true ){ 
    cout  << "Not a Number is found: " << c << endl;
    exit(1);
  }
}
