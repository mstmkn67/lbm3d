#ifndef _LBM_3D_RECT_H_
#define _LBM_3D_RECT_H_

#include "LBM3dSimulator.h"

class LBM3dRect:public LBM3dSimulator{
public:
	LBM3dRect(LBM3dLattice* lattice,
	          double viscosity,double density,double dt,
	          double delta_p=0.0);
	virtual ~LBM3dRect();
	virtual void collide();
protected:
	virtual void collide_down_wall(int i,int k);
	virtual void collide_up_wall(int i,int k);
	virtual void collide_back_wall(int i,int j);
	virtual void collide_front_wall(int i,int j);
	virtual void nonslip_down_wall(int i,int k);
	virtual void nonslip_up_wall(int i,int k);
	virtual void nonslip_back_wall(int i,int j);
	virtual void nonslip_front_wall(int i,int j);
	//
	virtual void collide_down_back_edge(int i);
	virtual void collide_down_front_edge(int i);
	virtual void collide_up_back_edge(int i);
	virtual void collide_up_front_edge(int i);
	virtual void nonslip_down_back_edge(int i);
	virtual void nonslip_down_front_edge(int i);
	virtual void nonslip_up_back_edge(int i);
	virtual void nonslip_up_front_edge(int i);
	//
	virtual void set_boundary_force();
	virtual void set_boundary_viscosity();
private:
};


#endif // _LBM_3D_PLATES_H_
