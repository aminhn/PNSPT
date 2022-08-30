# --- SYSTEM ---

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
BASISDIR    = /opt/ibm/ILOG

# --- DIRECTORIES ---

CCC = g++ -std=gnu++11 -no-pie -Iincludes
BASISILOG  = $(shell find $(BASISDIR) -maxdepth 1 -type d -name "CPLEX_Studio*" | sort -V | tail -1)
CONCERTDIR = $(BASISILOG)/concert
CPLEXDIR   = $(BASISILOG)/cplex

# --- FLAGS ---

CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DIL_STD -Wno-deprecated-declarations -Wno-ignored-attributes
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

# --- OPTIMIZATION FLAGS ---

DEBUG_OPT = -DNDEBUG -O3
#DEBUG_OPT = -g3 -O0
#PROF = -pg
PROF =

CFLAGS += $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I./include $(DEBUG_OPT) -c $(PROF)

LDFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -ldl

# ---- COMPILE  ----
SRC_DIR1   := src
OBJ_DIR1   := obj

SRC_DIRS1  := $(shell find $(SRC_DIR1) -type d)
OBJ_DIRS1  := $(addprefix $(OBJ_DIR1)/,$(SRC_DIRS1))

SOURCES1   := $(shell find $(SRC_DIR1) -name '*.cpp')
OBJ_FILES1 := $(addprefix $(OBJ_DIR1)/, $(SOURCES1:.cpp=.o))

vpath %.cpp $(SRC_DIRS1)


# ---- TARGETS ----

EXECUTABLE1 = mkd

EXECUTABLES = $(EXECUTABLE1) 

all: $(EXECUTABLES)

$(EXECUTABLE1): makedir $(SOURCES1) $(OBJ_FILES1) 
	$(CCC) $(OBJ_FILES1) $(LDFLAGS) $(PROF) -o $@

$(OBJ_DIR1)/%.o: %.cpp
	$(CCC) $(CFLAGS) $< -o $@

makedir: $(OBJ_DIRS1)

$(OBJ_DIRS1):
	@mkdir -p $@

clean:
	@rm -rf $(OBJ_DIR1)
	@rm -rf $(EXECUTABLE1)




