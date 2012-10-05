# change this to 1 or 0 depending if you want a pretty plot,
# not sure what you want otherwise
DOPLOT = 1
FPS	   = 1
POINTS = 0
IMAGES = 0

# standard compile options for the c++ executable
GCC = gcc
EXE = entbody
OBJS =  main.c
FLAGS = -O3 -Wall 
LIBFLAGS = -lm

ifeq ($(DOPLOT), 1)
	OBJS += plot.c
	FLAGS += -DPLOT
	LIBFLAGS += -lGL -lGLU -lglut 
endif

ifeq ($(FPS),1)
	LIBFLAGS += -lrt
	FLAGS += -DFPS
endif

ifeq ($(POINTS), 1)
	FLAGS += -DPOINTS
endif

ifeq ($(IMAGES), 1)
    FLAGS += -DIMAGES
    LIBFLAGS += -lIL -lILU -lILUT
endif

# default super-target
all: $(EXE)

# the standard executable
$(EXE): $(OBJS)
	$(GCC) $(FLAGS) $^ -o $@ $(LIBFLAGS)


.PHONY: clean

clean:
	rm -rf $(EXE)
