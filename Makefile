UNAME := $(shell uname)
NVCC := $(shell command -v nvcc 2> /dev/null)

DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic/src
DIR_H          = $(DIR_MAIN)rhic/include
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
FLOWTRACE =
LDFLAGS=
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS) -Wno-comment

ifdef NVCC
COMPILER = nvcc
OPTIMIZATION = -O5
OPTIONS := $(OPTIONS) --relocatable-device-code=true -Wno-deprecated-gpu-targets
LINK_OPTIONS := $(LINK_OPTIONS) --cudart static --relocatable-device-code=true -link -Wno-deprecated-gpu-targets
endif
ifndef NVCC
COMPILER = gcc
OPTIMIZATION = -O3
endif

ifeq ($(UNAME), Linux)
LIBS = -lm -lgsl -lgslcblas -lconfig
endif
ifeq ($(UNAME), Darwin)
LIBS = -L /usr/local/lib -lm -lgsl -lgslcblas -lconfig -largp -lc++
endif

INCLUDES = -I rhic/include -I freezeout

CPP := $(shell find $(DIR_SRC) -name '*.cpp' -and -not -name '*Test.cpp' )
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ)

EXE = cpu-vh

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

TEST_CPP := $(shell find $(DIR_SRC) -name '*.cpp' -and -not -name '*Run.cpp')
TEST_CPP_OBJ  = $(TEST_CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
TEST_OBJ = $(TEST_CPP_OBJ)

TEST_EXE = cpu-vh-test

$(TEST_EXE): $(TEST_OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

all: $(EXE)

test: $(TEST_EXE)
	echo "Testing:   $(TEST_EXE)"
	$(DIR_MAIN)$(TEST_EXE)

hydro: $(EXE)
	echo "Running Hydro:   $(EXE)"
	rm -rf output
	mkdir output
	$(DIR_MAIN)$(EXE) --config=rhic-conf --output=output --hydro

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(TEST_EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:
