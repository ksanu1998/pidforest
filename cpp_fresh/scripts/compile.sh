rm -r *.o
rm -r pidforest

echo "start compiling"
g++ -std=c++17 -c ./TimeSeries.cpp -o TimeSeries.o
g++ -std=c++17 -c ./Histogram.cpp -o Histogram.o
g++ -std=c++17 -c ./Cube.cpp -o Cube.o
g++ -std=c++17 -c ./Node.cpp -o Node.o
g++ -std=c++17 -c ./PointSet.cpp -o PointSet.o
g++ -std=c++17 -c ./Forest.cpp -o Forest.o
g++ -std=c++17 -c ./HistogramMemory.cpp -o HistogramMemory.o
g++ -std=c++17 -c ./Split.cpp -o Split.o
g++ -std=c++17 -c ./main.cpp -o main.o
g++ -std=c++17 -o pidforest TimeSeries.o Histogram.o Cube.o Node.o PointSet.o Forest.o HistogramMemory.o Split.o main.o