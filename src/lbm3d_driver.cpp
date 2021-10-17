#include "udf/gourmain.h"
#include "LBM3d.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void udfHeaderCheck()
{
	string version("1.0"),engine("lbm3d");
	cout << "**************************************************************" << endl;
	cout <<  "              " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: lbm3d -I input_udf [-O output_udf] [-n number of thread]" << endl;
}


int gourmain(UDFManager* in,UDFManager* out,int n){
	udfHeaderCheck();
#ifdef _OPENMP
	omp_set_num_threads(n);
#endif
	LBM3d* sim=new LBM3d(in,out);
	sim->update();
	delete sim;
	return 0;
}
