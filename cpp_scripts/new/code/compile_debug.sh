rm -r *.o
rm -r pidforest

echo " >> Compiling"
echo " >> Compiling TimeSeries"
g++ -g -std=c++17 -c ./TimeSeries.cpp -o TimeSeries.o
echo " >> Compiling TimeSeries [DONE]"
echo " >> Compiling Histogram"
g++ -g -std=c++17 -c ./Histogram.cpp -o Histogram.o
echo " >> Compiling Histogram [DONE]"
echo " >> Compiling Cube"
g++ -g -std=c++17 -c ./Cube.cpp -o Cube.o
echo " >> Compiling Cube [DONE]"
echo " >> Compiling Node"
g++ -g -std=c++17 -c ./Node.cpp -o Node.o
echo " >> Compiling Node [DONE]"
echo " >> Compiling PointSet"
g++ -g -std=c++17 -c ./PointSet.cpp -o PointSet.o
echo " >> Compiling PointSet [DONE]"
echo " >> Compiling Forest"
g++ -g -std=c++17 -c ./Forest.cpp -o Forest.o
echo " >> Compiling Forest [DONE]"
echo " >> Compiling HistogramMemory"
g++ -g -std=c++17 -c ./HistogramMemory.cpp -o HistogramMemory.o
echo " >> Compiling HistogramMemory [DONE]"
echo " >> Compiling Split"
g++ -g -std=c++17 -c ./Split.cpp -o Split.o
echo " >> Compiling Split [DONE]"
echo " >> Compiling main"
g++ -g -std=c++17 -c ./main_debug.cpp -o main.o
echo " >> Compiling main [DONE]"
echo " >> Linking files and creating executable"
g++ -g -std=c++17 -o pidforest TimeSeries.o Histogram.o Cube.o Node.o PointSet.o Forest.o HistogramMemory.o Split.o main.o
echo " >> Linking files and creating executable [DONE]"
echo " >> Compiling [DONE]"