SRCS=src/util.cpp src/random.cpp src/block_file.cpp src/b_node.cpp \
src/b_tree.cpp src/qalsh.cpp src/ann.cpp src/main.cpp

OBJS=$(SRCS:.cpp=.o)

# CXX=g++
# CXXFLAGS=-g -Wall -MMD -std=c++11
# CPPFLAGS=-w -O0

COMPILER=g++
COMP_FLAGS=-w -O0 -g -Wall -MMD -std=c++11

.PHONY: clean

all: $(OBJS)
	$(COMPILER) $(COMP_FLAGS) -o qalsh $(OBJS)

util.o: src/util.h

random.o: src/random.h

block_file.o: src/block_file.h

b_node.o: src/b_node.h

b_tree.o: src/b_tree.h

qalsh.o: src/qalsh.h

ann.o: src/ann.h

main.o:

clean:
	-rm $(OBJS) qalsh
