# Installation directories for pylon
PYLON_ROOT ?= /opt/pylon

# Build tools and flags
CPPFLAGS   := $(shell $(PYLON_ROOT)/bin/pylon-config --cflags) 
GLIBS := gtk+-3.0
CPPFLAGS   += $(shell pkg-config --cflags $(GLIBS)) -I./include
CPPFLAGS   += -g
LDFLAGS    := $(shell $(PYLON_ROOT)/bin/pylon-config --libs-rpath)
LDLIBS     := $(shell $(PYLON_ROOT)/bin/pylon-config --libs)
LDLIBS     += $(shell pkg-config --libs $(GLIBS))
LDLIBS     += -lpng -lrt -lmatio -lpthread -lGL -lGLEW -lgsl -lgslcblas -lfreetype 
LDLIBS	 += -L/AlpaoLib/x64/ -lasdk

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
