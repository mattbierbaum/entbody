GCC = gcc
EXE = $(NAME) 
SRC = $(ENTBODY)/main.c $(ENTBODY)/util.c 
HDR = $(FULLPATH).h
FLAGS = -O3 -DHEADER=\"$(HDR)\" 
LIBFLAGS = -lm

# we want the compile line to be essentially
# nvcc main.cu -arch sm_12 -DPLOT -lGL -lGLU -lglut
ifeq ($(CUDA), 1)
    FLAGS += -x cu -DCUDA -arch sm_11
    GCC = nvcc
else
    SRC += $(HDR)
endif 

ifeq ($(DOPLOT), 1)
    SRC += $(ENTBODY)/addons/plot.c
    FLAGS += -DPLOT
    LIBFLAGS += -lGL -lGLU -lglut
endif

ifeq ($(OPENMP), 1)
    FLAGS += -fopenmp
    FLAGS += -DOPENMP
endif

ifeq ($(FPS),1)
    LIBFLAGS += -lrt
    FLAGS += -DFPS
endif

ifeq ($(POINTS), 1)
    FLAGS += -DPOINTS
endif

ifeq ($(TEMPS), 1)
    FLAGS += -save-temps
endif

# default super-target
all: $(EXE)

$(HDR): 

# the standard executable
$(EXE): $(SRC) 
	$(GCC) $(FLAGS) $^ -o $@ $(LIBFLAGS)

clean: 
	rm -rf $(EXE)

.PHONY: clean
