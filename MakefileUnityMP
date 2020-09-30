DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic/src
DIR_H          = $(DIR_MAIN)rhic/include
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
OPTIMIZATION = -O3
FLOWTRACE =


#Different options for different compilers and optimizations

#OPTIONS = -fopenmp -march=native -fopt-info-vec #-funroll-loops #-static-libstdc++ # for g++ with vector info 
#OPTIONS = -qopenmp -std=c++11 -lhdf5 -lhdf5_cpp # for icpc
OPTIONS = -fopenmp -march=native -std=c++11 -lhdf5 -lhdf5_cpp # for g++


#link against libconfig and googletest

#LINK_OPTIONS = -L/home/everett.165/gperftools/lib -lprofiler  #link against gperftools for cpu profiling 
LINK_OPTIONS = -L/home/du.458/libconfig/lib/.libs -lconfig -L/home/du.458/googletest/googletest/mybuild/ -lgtest


CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)


#choose a compiler

#COMPILER = icpc
COMPILER = g++


#choose libs for different compilers

#LIBS = -lgsl -lgslcblas # for icpc  
LIBS = -lm -lgsl -lgslcblas -lgomp -lconfig -lgtest -lhdf5 -lhdf5_cpp # for g++ 


INCLUDES = -I rhic/include -I rhic/freezeout -I /home/du.458/libconfig/lib/ -I /home/du.458/googletest/googletest/include/


CPP := $(shell find $(DIR_SRC) -name '*.cpp')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ)


EXE =\
beshydro

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(OPTIONS) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:

