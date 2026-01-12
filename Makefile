ROOT_CONFIG := /home/ubuntu/Vinh/root_6_28/bin/root-config

CXX := g++
CXXFLAGS := -std=c++20 -Wall -O2 -fPIC \
            -Iinclude \
            -MMD -MP \
            -DROOT_STRICT_OWNERSHIP \
            $(shell $(ROOT_CONFIG) --cflags)

LDFLAGS := $(shell $(ROOT_CONFIG) --ldflags)
LDLIBS  := $(shell $(ROOT_CONFIG) --libs)

TARGET := main
BUILDDIR := build

SRC := main.cpp $(wildcard src/*.cpp)

# Превращаем src/foo.cpp → build/src/foo.o
OBJ := $(SRC:%.cpp=$(BUILDDIR)/%.o)
DEP := $(OBJ:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# build/src/foo.o ← src/foo.cpp
$(BUILDDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(DEP)

clean:
	rm -rf $(BUILDDIR) $(TARGET)
