#ifndef _LBM3D_H_
#define _LBM3D_H_

#include "udfmanager.h"
#include "lbm/LBM3dSimulator.h"
#include "Timer.h"
#include <cmath>
#include <sstream>
using namespace std;

class LBM3d{
public:
	LBM3d(UDFManager* in,UDFManager* out);
	virtual ~LBM3d();
	virtual void update();
protected:
	virtual void input_lbm();
	virtual void input_other();
	virtual void output();
	virtual void nan_check(double x,char c='0');
private:
	//virtual void field_test();
	UDFManager *in,*out;
	Timer* timer;
	int iteration;
	//flow
	LBM3dSimulator *lbm;
	LBM3dLattice *lattice;
	//output
	bool fou,fop,forho,fof;
};

#endif // _MIST_SIMULATOR3D_H_
