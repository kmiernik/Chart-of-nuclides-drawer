CPP = g++
CPPFLAGS = -Wall
#Source dir
SDIR = src
#Header dir
HDIR = include

#Rule to make .o from .cpp files
%.o: $(SDIR)/%.cpp
	$(CPP) $(CPPFLAGS) -I $(HDIR) -c $< -o $@

all: chartdraw 

chartdraw: nuclide.o
	$(CPP) $(CPPFLAGS) -o $@ nuclide.o

clean: 
	rm -f *.o *~ include/*~ src/*~ chartdraw
