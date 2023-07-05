#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../Histogram.h"
#include "../Split.h"
//#include "../HistogramMemory"

namespace py = pybind11;

class PForest
{
    public:
        PForest(double epsilon, int bucket) : epsilon(epsilon), bucket(bucket) {}
        ~PForest() {}
		int compute_histogram(py::array x) {
		
			py::buffer_info info = x.request();
			auto ptr = static_cast<double *>(info.ptr);
			if (info.ndim > 1) {
				std::cerr<<"For now can only compute one dimensional arrays";
				return -1;
			}
			int array_size = 1;
    			for (auto r: info.shape) {
      			array_size *= r;
    		}

		
			Histogram H = Histogram(bucket, epsilon);
			Split S;
			HistogramMemory M  = HistogramMemory(array_size, array_size, bucket);
			H.ComputeAppBuckets(ptr, array_size, &M, &S);
			std::cout<<"Printing from split "<<S.bucket_level_ + 1<<" , "<<S.bucket_error_<<std::endl;
    		for (int k = 0; k < S.bucket_level_ + 1; k++)
        		std::cout<<S.bucket_list_[k]<<", ";
    		std::cout<<"\nHurray\n";	
    		return 1;		
		}
		
		double get_epsilon() { return epsilon; }
        
    private:
    	double epsilon;
    	int bucket;
};


PYBIND11_MODULE(pforest, m) {
// optional module docstring
    m.doc() = "pybind11 example plugin";
    
	py::class_<PForest>(m, "PForest")
        .def(py::init<double , int >())
        .def("get_epsilon", &PForest::get_epsilon)
        .def("compute_histogram", &PForest::compute_histogram);
}

