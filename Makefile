# Installation directories for pylon
PYLON_ROOT ?= /opt/pylon

# Build tools and flags
CPPFLAGS   := $(shell $(PYLON_ROOT)/bin/pylon-config --cflags) 
GLIBS := gtk+-3.0
CPPFLAGS   += $(shell pkg-config --cflags $(GLIBS)) -I./include
CPPFLAGS   += -O3
LDFLAGS    := $(shell $(PYLON_ROOT)/bin/pylon-config --libs-rpath)
LDLIBS     := $(shell $(PYLON_ROOT)/bin/pylon-config --libs) 
LDLIBS     += $(shell pkg-config --libs $(GLIBS))
LDLIBS     += -lasdk -lpng -lrt -lmatio -lpthread -lGL -lGLEW -lgsl -lgslcblas -lfreetype 


OBJS := src/video-main.o src/gl-main.o src/deformable_mirror.o src/gettime.o src/writematlab.o src/text_helper.o

# Rules for building
all: shwfs

%.o : %.cpp
	g++ -c $(CPPFLAGS) $< -o $@

shwfs: $(OBJS)
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -rf $(OBJS)
	rm -rf shwfs

deps: 
	sudo apt-get install libpng-dev libgtk-3-dev libmatio-dev \
	libglew-dev libfreetype-dev libgsl-dev freeglut3-dev \
	libglm-dev 
	fallocate -l 1M shared_centroids.dat
	fallocate -l 1k shared_dmctrl.dat
