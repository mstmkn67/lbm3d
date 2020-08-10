#ifndef _LBM_3D_PLATES_H_
#define _LBM_3D_PLATES_H_

#include "LBM3dSimulator.h"

class LBM3dPlates:public LBM3dSimulator{
public:
	LBM3dPlates(LBM3dLattice* lattice,
	            double viscosity,double density,double dt,
	            double delta_p=0.0,
	            double u_down=0.0,double v_down=0.0,double w_down=0.0,
	            double u_up=0.0,  double v_up=0.0,  double w_up=0.0);
	virtual ~LBM3dPlates();
	virtual void collide();
protected:
	virtual void collide_down_wall(int i,int k);
	virtual void collide_up_wall(int i,int k);
	virtual void nonslip_down_wall(int i,int k);
	virtual void nonslip_up_wall(int i,int k);
	//virtual void collide_edge_x(int i);
	//virtual void collide_edge_y(int j);
	//virtual void collide_edge_z(int k);
	//virtual void collide_corner();
	//virtual void set_boundary_distribution();
	virtual void set_boundary_force();
	virtual void set_boundary_viscosity();
private:
	double u0,v0,w0;
	double u2,v2,w2;
};


#endif // _LBM_3D_PLATES_H_
