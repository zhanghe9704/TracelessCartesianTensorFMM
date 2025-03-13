# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -O3 -std=c++17

# Source files
SRC = $(wildcard *.cc)

# Object files
OBJ = $(SRC:.cc=.o)

# Executable name
TARGET = sdafmm

# Default rule
all: $(TARGET)

# Link the final executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJ) $(TARGET)

# Phony targets (not actual files)
.PHONY: all clean
