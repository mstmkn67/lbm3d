#include "udf/gourmain.h"
#include "LBM3d.h"

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
	cout << "usage: lbm3d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	LBM3d* sim=new LBM3d(in,out);
	sim->update();
	delete sim;
	return 0;
}
