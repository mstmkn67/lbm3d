#ifndef _LBM_3D_SIMULATOR_H_
#define _LBM_3D_SIMULATOR_H_

#include "LBM3dLattice.h"
#include <cmath>
#include <stdlib.h>
using namespace std;
#ifdef _OPENMP
#include <omp.h>
#endif


class LBM3dSimulator{
public:
	LBM3dSimulator(LBM3dLattice* lattice,
	               double viscosity,double density,double dt,double delp=0.0);
	virtual ~LBM3dSimulator();
	virtual void initial();
	virtual void update();
	//
	virtual void set_bulk_viscosity();
	virtual void reset_bforce();
	//
	LBM3dLattice* lattice;
	//
	double delp;
protected:
	virtual void collide()=0;
	virtual void calc_density_velocity();
	virtual void calc_collide(int n,int i,int j,int k);
	virtual void time_update();
	virtual void set_boundary_force()=0;
	virtual void set_boundary_viscosity()=0;
	//
	double viscosity;
	double density;
	double dt;
	//
	double c,cs2;
	double w[15];
	int cx[15],cy[15],cz[15];
};

#endif // _LBM_3D_SIMULATOR_H_
