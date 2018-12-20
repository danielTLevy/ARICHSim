# Declare
#
# C FLAGS
CC := g++
CC_FLAGS =-std=c++11 -O2 -ffloat-store #-ansi -Wall -g
#INC := -I/Users/Daniel/Documents/courses/phys449/boost/boost_1_67_0/
#LIB := -L. -lMinuit -L/Users/Daniel/Documents/courses/phys449/boost/boost_1_67_0/libs/iostreams
#ROOT FLAGS
R_LDFLAGS :=     `root-config --ldflags`
R_LIBS    :=     `root-config --glibs`
R_CFLAGS  :=     `root-config --cflags`
R_ALL     :=     $(R_LADFLAGS) $(R_LIBS) $(R_CFLAGS)

# File names
EXEC = run
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

$(EXEC): $(OJECTS)
	$(CC)  $(CC_FLAGS) -fPIC $(LDFLAGS) $(R_ALL) $(B_FLAGS)  -c $(SOURCES)
	$(CC)  $(P_CC_FLAG) $(CC_FLAGS) $(R_ALL) $(G_FLAGS) $(B_FLAGS) $(OBJECTS) -o final.out

# To remove generated files
clean:
	rm *.o

lib:
	make clean
	make
	"g++" -fPIC $(CC_FLAGS)$ $(R_ALL) $(LIB) $(INC)  $(OBJECTS) -shared -o libPID.so
