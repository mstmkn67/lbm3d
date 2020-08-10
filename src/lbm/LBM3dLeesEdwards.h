#ifndef _LBM_3D_LEES_EDWARDS_H_
#define _LBM_3D_LEES_EDWARDS_H_

#include "LBM3dSimulator.h"

class LBM3dLeesEdwards:public LBM3dSimulator{
public:
	LBM3dLeesEdwards(LBM3dLattice* lattice,
                  double viscosity,double density,double dt,
                  int* step,
                  double shear_rate=0.0);
	virtual ~LBM3dLeesEdwards();
	virtual void collide();
protected:
	virtual void collide_up_down_periodic(int i,int k);
	//virtual void collide_corner();
	//
	virtual void periodic_down(int n,int i,int k);//n=
	virtual void periodic_up(int n,int i,int k);//n=
	//virtual void set_boundary_distribution();
	virtual void set_boundary_force();
	virtual void set_boundary_viscosity();
private:
	int* step;
	double gdot;
	//
	double Du;
	double a0,a1;
	int Dn;
};


#endif // _LBM_3D_LEES_EDWARDS_H_
