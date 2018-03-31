CC =g++;
CFLAGS = -lm -lpthread -lpng -lz -O3 -ftree-vectorize
DEBUG =  -pg
TARGETS = fract

# Mark the default target to run (otherwise make will select the first target in the file)
.DEFAULT: all
# Mark targets as not generating output files (ensure the targets will always run)
.PHONY: all debug clean

all: clean $(TARGETS)

clean:
	rm -f $(TARGETS).exe *.png

cleaner:
	rm -f *.png slurm*

# A debug target to update flags before cleaning and compiling all targets
debug: CFLAGS += $(DEBUG)
debug: clean $(TARGETS)

fract: fractal.cpp CImg.h
	$(CC) $(CFLAGS) fractal.cpp CImg.h -fopenmp -o fract.exe