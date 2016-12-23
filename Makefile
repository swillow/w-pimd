#use these variable to set if we will use mpi or not 

AR=ar
ARFLAGS=-qs
RANLIB=ranlib

EXECUTABLE = wpimd1.x wpimd2.x

CXX = c++
HOME = .
.SUFFIXES: .cc 


CFLAGS = -g     -ffast-math    -std=c++11  
LIBS =   \
	-larmadillo -lblas -llapack 


MD_SRC = wpimd1.cc wpimd2.cc 

MD_OBJ=$(MD_SRC:.cc=.o)

.cc.o :
	$(CXX) $(CFLAGS) -c $< -o $@

all	: wpimd1.x wpimd2.x


wpimd1.x : wpimd1.o
	$(CXX)   $(FLAGS) -o  wpimd1.x wpimd1.o $(LIBS) 

wpimd2.x : wpimd2.o
	$(CXX)   $(FLAGS) -o  wpimd2.x wpimd2.o $(LIBS) 

clean:
	rm *.o *~ *.x

# DO NOT DELETE
