CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o lbm3d  -L $(PF_FILES)/lib/linux_64 -lplatform

OBJECT= Timer.o gourmain.o cmdlineopts.o \
	LBM3dSimulator.o LBM3dPlates.o LBM3dRect.o LBM3dLeesEdwards.o \
	LBM3dLattice.o  \
	LBM3d.o lbm3d_driver.o

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Timer.o : Timer.cpp Timer.h
	$(CC) Timer.cpp $(CFLAG)

gourmain.o : udf/gourmain.cpp udf/gourmain.h
	$(CC) udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : udf/cmdlineopts.cpp udf/cmdlineopts.h
	$(CC) udf/cmdlineopts.cpp $(CFLAG)

LBM3dSimulator.o : lbm/LBM3dSimulator.cpp lbm/LBM3dSimulator.h
	$(CC) lbm/LBM3dSimulator.cpp $(CFLAG)

LBM3dPlates.o : lbm/LBM3dPlates.cpp lbm/LBM3dPlates.h
	$(CC) lbm/LBM3dPlates.cpp $(CFLAG)

LBM3dRect.o : lbm/LBM3dRect.cpp lbm/LBM3dRect.h
	$(CC) lbm/LBM3dRect.cpp $(CFLAG)

LBM3dLeesEdwards.o : lbm/LBM3dLeesEdwards.cpp lbm/LBM3dLeesEdwards.h
	$(CC) lbm/LBM3dLeesEdwards.cpp $(CFLAG)

LBM3dLattice.o : lbm/LBM3dLattice.cpp lbm/LBM3dLattice.h
	$(CC) lbm/LBM3dLattice.cpp $(CFLAG)

LBM3d.o : LBM3d.cpp LBM3d.h
	$(CC) LBM3d.cpp $(CFLAG)

lbm3d_driver.o : lbm3d_driver.cpp
	$(CC) lbm3d_driver.cpp $(CFLAG)

clean : 
	rm -f *.o lbm3d

install :
	strip lbm3d
	cp lbm3d ../bin
