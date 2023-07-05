echo "start compiling"
# g++ -std=c++17 -c -o Histogram.o ./Histogram.cpp
# g++ -std=c++17 -c -o Node.o ./Node.cpp
# g++ -std=c++17 -c -o Cube.o ./Cube.cpp
# g++ -std=c++17 -c -o PointSet.o ./PointSet.cpp
# g++ -std=c++17 -c -o Forest.o ./Forest.cpp
# g++ -std=c++17 -c -o HistogramMemory.o ./HistogramMemory.cpp
# g++ -std=c++17 -c -o Split.o ./Split.cpp
# g++ -std=c++17 -c -o main.o ./main.cpp
# g++ -std=c++17 -o pidforest Histogram.o Node.o Cube.o PointSet.o Forest.o HistogramMemory.o Split.o main.o
g++ -std=c++17 -o pidforest Histogram.o Node.o Cube.o PointSet.o Forest.o HistogramMemory.o Split.o main.o
