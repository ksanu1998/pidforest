# Makefile for compiling the project

# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -c -g

# Source files
SOURCES = TimeSeries.cpp Histogram.cpp Cube.cpp Node.cpp PointSet.cpp Forest.cpp HistogramMemory.cpp Split.cpp main.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Executable
EXECUTABLE = pidforest

# Default target
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	@echo " >> Linking files and creating executable"
	$(CXX) -std=c++17 -o $(EXECUTABLE) $(OBJECTS)
	@echo " >> Linking files and creating executable [DONE]"

%.o: %.cpp
	@echo " >> Compiling $<"
	$(CXX) $(CXXFLAGS) $< -o $@
	@echo " >> Compiling $< [DONE]"

clean:
	@echo "Cleaning up..."
	rm -rf $(EXECUTABLE) $(OBJECTS)
	@echo "Cleaned up [DONE]"
