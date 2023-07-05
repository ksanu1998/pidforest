echo "start compiling"
g++ -c -o Histogram.o ./src/Histogram.cpp
g++ -c -o HistogramMemory.o ./src/HistogramMemory.cpp
g++ -c -o Node.o ./src/Node.cpp
g++ -c -o Split.o ./src/Split.cpp
g++ -c -o Tree_Param.o ./src/Tree_Param.cpp
g++ -c -o main.o ./src/main.cpp
g++ -o pidforest Histogram.o HistogramMemory.o Node.o Split.o Tree_Param.o main.o -lm
