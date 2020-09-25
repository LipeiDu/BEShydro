DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic/src
DIR_H          = $(DIR_MAIN)include/
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
OPTIMIZATION = -O5  
FLOWTRACE =
OPTIONS =
LINK_OPTIONS = -link -L/home/du.458/libconfig/lib/.libs -lconfig -L/home/du.458/googletest/googletest/mybuild/ -lgtest
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)
COMPILER = nvcc
LIBS = -lm -lgsl -lgslcblas -lconfig -lgtest -lhdf5 -lhdf5_cpp
INCLUDES = -I rhic/include -I /home/du.458/libconfig/lib/ -I /home/du.458/googletest/googletest/include/ -I rhic/freezeout

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ) 

EXE =\
beshydro

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	       @[ -d $(DIR_OBJ) ] || find rhic/ -type d -exec mkdir -p ./build/{} \;
	       @echo "Compiling: $< ($(COMPILER))"
	       $(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:
