rm -r pidforest

echo " >> Compiling"
echo " >> Compiling TimeSeries"
g++ -std=c++17 -c -g ./TimeSeries.cpp -o TimeSeries.o
echo " >> Compiling TimeSeries [DONE]"
echo " >> Compiling Histogram"
g++ -std=c++17 -c -g ./Histogram.cpp -o Histogram.o
echo " >> Compiling Histogram [DONE]"
echo " >> Compiling Cube"
g++ -std=c++17 -c -g ./Cube.cpp -o Cube.o
echo " >> Compiling Cube [DONE]"
echo " >> Compiling Node"
g++ -std=c++17 -c -g ./Node.cpp -o Node.o
echo " >> Compiling Node [DONE]"
echo " >> Compiling PointSet"
g++ -std=c++17 -c -g ./PointSet.cpp -o PointSet.o
echo " >> Compiling PointSet [DONE]"
echo " >> Compiling Forest"
g++ -std=c++17 -c -g ./Forest.cpp -o Forest.o
echo " >> Compiling Forest [DONE]"
echo " >> Compiling HistogramMemory"
g++ -std=c++17 -c -g ./HistogramMemory.cpp -o HistogramMemory.o
echo " >> Compiling HistogramMemory [DONE]"
echo " >> Compiling Split"
g++ -std=c++17 -c -g ./Split.cpp -o Split.o
echo " >> Compiling Split [DONE]"
echo " >> Compiling main"
g++ -std=c++17 -c -g ./main.cpp -o main.o
echo " >> Compiling main [DONE]"
echo " >> Linking files and creating executable"
g++ -std=c++17 -o pidforest TimeSeries.o Histogram.o Cube.o Node.o PointSet.o Forest.o HistogramMemory.o Split.o main.o
echo " >> Linking files and creating executable [DONE]"
echo " >> Compiling [DONE]"
rm -r *.o