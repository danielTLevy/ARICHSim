INCDIR= $(shell pwd)/inc
SRCDIR= $(shell pwd)/src
OBJDIR= $(shell pwd)/obj
BINDIR= $(shell pwd)/bin

VPATH = $(SRCDIR)

# FLAGS
CC := g++
CFLAGS=-c -g -Wall `root-config --cflags` -I${INCDIR}
LDFLAGS=`root-config --glibs` -lHistPainter

# File names
EXEC = $(BINDIR)/final.out
FILES= $(wildcard $(SRCDIR)/*.cpp)
SOURCES=$(FILES)

OBJECTS = $(FILES:$(SRCDIR)/%.cpp=${OBJDIR}/%.o)
print-%  : ; @echo $* = $($*)

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

# To remove generated files
clean:
	rm $(BINDIR)/* $(OBJDIR)/*.o

