CXX = /usr/local/Cellar/gcc\@9/9.3.0/bin/g++-9

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors #-Werror

dirLib = /Users/cgiocoli/lib/CosmoBolognaLib/
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/eigen-3.3.4/
dir_CCfits = /Users/cgiocoli/lib/CCfits/include/
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL -L/Users/cgiocoli/lib/gsl-2.3/lib/ -L/Users/cgiocoli/lib/fftw-3.3.4-mpi/lib/ -L/Users/cgiocoli/lib/gsl-2.3/lib/ -lgsl nfwLens.cpp
FLAGS_INC = -I./ -I/Users/cgiocoli/lib/gsl-2.3/include/ -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits) -I/Users/cgiocoli/lib/gsl-2.3/include/ -I/Users/cgiocoli/lib/fftw-3.3.4-mpi/include/ 

OBJ = fit_DeltaSigma-2.2.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

fit_DeltaSigma-2.2: $(OBJ) 
	$(CXX) $(OBJ) -o fit_DeltaSigma-2.2_1h_gamma $(FLAGS_LIB)  $(FLAGS_INC)

clean:
	rm -f *.o fit_DeltaSigma-2.2 *~ \#* temp* core*

fit_DeltaSigma-2.2.o: fit_DeltaSigma-2.2.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fit_DeltaSigma-2.2.cpp 

