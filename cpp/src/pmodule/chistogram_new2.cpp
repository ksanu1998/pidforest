#include <iostream>
#include <string>
#include <math.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include "../Histogram.h"
#include "../Split.h"
#include "../HistogramMemory.h"
#include "../util_functions.h"

namespace py = pybind11;
/**
* This class should compute a histogram for the python implementation
******************************/


class CHistogram_new2
{
    public:
        CHistogram_new2(double epsilon, int bucket) : epsilon_(epsilon), n_buckets_(bucket) {
			try{
        	if (n_buckets_ > MAX_BUCKETS) {
				py::print("number of buckets reduced to ", MAX_BUCKETS);
				n_buckets_ = MAX_BUCKETS;
			}
           /* arr_ = extract_array_ptr(arr);
            arr_size_ = extract_array_size(arr);
            counts_ = extract_array_ptr(counts);
           */ 

            hist_mem_ = NULL;
			/*double list_length;
    		if (utils::isEqual(this->epsilon_, 0)) list_length = arr_size_;
    		else list_length = log(arr_size_)/log(1 + epsilon_);
    		if (list_length > arr_size_) list_length = arr_size_; 
			hist_mem_  = new HistogramMemory(int(list_length), arr_size_, n_buckets_);
			
			total_counts_ = hist_mem_->total_count_;
			total_counts_[0] = 0;
			for (int k = 0; k < arr_size_; k++)
				total_counts_[k+1] = total_counts_[k] + counts_[k];
            */
			}//end try

			catch(const std::runtime_error& error) {
				py::print(error.what());
   			} 
		
        }

        ~CHistogram_new2() { 
			delete hist_mem_;
		}

        Split split_;
		
		HistogramMemory * hist_mem_;

		//double * total_counts_; // should be moved to hist_mem
        
    	/***************************************************
		* computes the histogram, fills the split object.
        * arr: array of unique data values
        * counts: array of counts per data value
        * b_list: array in which the output would be placed in place 
    	***************************************************/ 
		int compute_histogram(py::array_t<double> arr, py::array_t<double> counts, py::array_t<double> b_list) {

            py::buffer_info arr_info = arr.request();
			double * arr_ =  static_cast<double *>(arr_info.ptr);
            
            if (arr_info.ndim > 1) {
				throw std::runtime_error("Exception: array has to be one dimensional");
			}
			int arr_size_ = 1;
    			for (auto r: arr_info.shape) {
      			arr_size_ *= r;
    		}  
            py::buffer_info counts_info = counts.request();
            double * counts_ = static_cast<double *>(counts_info.ptr);
			verify_memory(arr_size_);
            
            double * total_counts_ = hist_mem_->total_count_;
			total_counts_[0] = 0;
			for (int k = 0; k < arr_size_; k++)
				total_counts_[k+1] = total_counts_[k] + counts_[k];
            
            double * sum = hist_mem_->sum_;
    		double * sq_sum = hist_mem_->sqsum_;
    		//double * sq_sum = new double[his_mem->max_length_];
    		int ** b_val = hist_mem_->b_val_; //buckets are [a,b]
    		double ** app_err_b = hist_mem_->app_err_b_; // the approximate error at b
    		int ** bucket_indices = hist_mem_->bucket_indices_;
    		int max_length = hist_mem_->max_length_;
    		int list_length[MAX_BUCKETS];
    		double curr_err_a[MAX_BUCKETS];

			/* Initializing variables */
    		for (int k = 0; k < n_buckets_; k++) {
        		list_length[k] = 0;
        		b_val[k][0] = 0;
        		app_err_b[k][0] = 0;
        		bucket_indices[k][0] = 0;
        		curr_err_a[k] = 0;
        		app_buckets_list_[k][0] = 0;
				//total_counts_[0] = 0;
    		}

    		sum[0] = arr_[0] * counts_[0];
    		sq_sum[0] = arr_[0] * arr_[0] * counts_[0];
			//sum[0] = arr_[0];
			//sq_sum[0] = arr_[0] * arr_[0];
    		double tmp = 0;

			/* Main loop on array items */
    		for (int j = 1; j < arr_size_; j++) {
        		//preparing the one bucket case. That's the recursion base.
				//total_counts_[j] = total_counts_[j-1] + counts_[j];
        		sum[j] = sum[j-1] + (arr_[j] * counts_[j]);
				//sum[j] = sum[j-1] + arr_[j] ;
        		sq_sum[j] = sq_sum[j - 1] + (arr_[j] * arr_[j] * counts_[j]);
				//sq_sum[j] = sq_sum[j - 1] + ( arr_[j] * arr_[j] );
        		//tmp = sq_sum[j] - ((sum[j] * sum[j]) / (double)(j + 1)); //error for 1 bucket
        		tmp = sq_sum[j] - ((sum[j] * sum[j]) / (double)total_counts_[j + 1]); //error for 1 bucket
				if (tmp < (1 + epsilon_) * curr_err_a[0]) { //segment ends in j
            		b_val[0][list_length[0]] = j;
        		} 
				else { // open a new segment
            		list_length[0]++;
            		b_val[0][list_length[0]] = j;
            		curr_err_a[0] = tmp;
            		bucket_indices[0][list_length[0]] = 0;
       	 		}
        		app_err_b[0][list_length[0]] = tmp;


        		// computing for buckets 2...
        		for (int k = 1; k < n_buckets_; k++) {
            		double tmp_err = sq_sum[j];
            		int curr_index = 0;
            		for (int i = 0; i <= list_length[k - 1]; i++) {
                		int b  = b_val[k-1][i];
                		if (tmp_err > app_err_b[k - 1][i] + IntervalError(sum, sq_sum, total_counts_, b + 1, j)) {
                    		tmp_err = app_err_b[k - 1][i] + IntervalError(sum, sq_sum, total_counts_, b + 1, j);
                    		curr_index = i;
                		}
            		}
            		if ((k < n_buckets_ - 1 ) and (tmp_err > (1 + epsilon_) * curr_err_a[k])) { //open a new segment
                		list_length[k]++;
                		int l = list_length[k];
                		b_val[k][l] = j;
                		app_err_b[k][l] = tmp_err;
                		bucket_indices[k][l] = curr_index;
                		curr_err_a[k] = tmp_err;
            		}
            		else {
                		int l = list_length[k];
                		b_val[k][l] = j;
                		app_err_b[k][l] = tmp_err;
                		bucket_indices[k][l] = curr_index;
            		}
        		}
    		}
    		for (int k = 0; k < n_buckets_; k++) {
        		app_buckets_err_[k] = app_err_b[k][list_length[k]];
   	 		}

   			// Now find the best split:
   			double var_red = 0;
   			int best_bucket = 0;
  			if (utils::isEqual(app_buckets_err_[0], 0)) {
      			split_.bucket_error_ = 0;
      			split_.bucket_level_ = 0;
  			}
  			else {
      			for (int k = 1; k < n_buckets_; k++) {
          			double err_red = (app_buckets_err_[0] - app_buckets_err_[k]) / app_buckets_err_[0];
          			if (err_red > var_red) {
              			var_red = err_red;
              			best_bucket = k;
          			}
      			}
          		split_.bucket_error_ = var_red;
          		split_.bucket_level_ = best_bucket;
  			}
  			ComputeAppList(bucket_indices, b_val, list_length[best_bucket], best_bucket, split_.bucket_list_);
            py::buffer_info blist_info = b_list.request();
			double * blist_ptr = static_cast<double *>(blist_info.ptr);
			for (int k = 0; k < split_.bucket_level_ + 1; k++)
        		blist_ptr[k] = split_.bucket_list_[k];
            blist_ptr[n_buckets_] = split_.bucket_error_;
  			return split_.bucket_level_ + 1;	
		}
		
		double get_epsilon() { 
			return epsilon_; 
		}

		void print_split() {
			py::print("Printing from split ", split_.bucket_level_ + 1, " , ", split_.bucket_error_);
    		for (int k = 0; k < split_.bucket_level_ + 1; k++)
        		py::print(split_.bucket_list_[k], ", ");
			py::print(" ");
		}

		int get_bucket() { return split_.bucket_level_ + 1; }

		double get_error() { return split_.bucket_error_; }

        void set_mem(int mem) { hist_mem_ = reinterpret_cast <HistogramMemory *> (mem); }

		/**
		 * Gets a numpy array as argument and returns the list of buckets by changing its values in place
		 */

		void get_bucket_list(py::array b_list) {
			py::buffer_info arr_info = b_list.request();
			double * list_ptr = static_cast<double *>(arr_info.ptr);
			for (int k = 0; k < split_.bucket_level_ + 1; k++)
        		list_ptr[k] = split_.bucket_list_[k];  
		}
    
        double epsilon_;
    	int n_buckets_;
    private:

		//double * arr_;
		//double * counts_;
		//int arr_size_;
		int  app_buckets_list_[MAX_BUCKETS][MAX_BUCKETS];
		double  app_buckets_err_[MAX_BUCKETS]; //an array with the app error accumulated per number of buckets used

		void ComputeAppList(int **indices, int **b_val, int list_length, int bucket, int *out_list) {
    		int curr_level = bucket;
    		int loc = list_length;
    		while (curr_level > 0) {
        		out_list[curr_level] = b_val[curr_level][loc];
        		loc = indices[curr_level][loc];
        		curr_level--;
    		}
    		out_list[curr_level] = b_val[curr_level][loc];
		}

		// computes the error of a bucket  in indexes j,...,i
		double IntervalError(double *sum, double *sqsum, double * total_counts_, int j, int i) {
    		double lb_s, lb_sq;
			//double total_pts  = i - j + 1;
			//double total_pts = 0;
    		if (i <= j) { return 0; }
    		if (j == 0) { 
				lb_s = 0; lb_sq = 0;
			}
        	else {
            	lb_s = sum[j-1]; 
            	lb_sq = sqsum[j-1]; 
        	}
			double total_pts = total_counts_[i+1] - total_counts_[j];   
    		double avg = (sum[i] - lb_s)/total_pts;
    		if (utils::isEqual(avg, 0))
        		avg = 0;
    		double ans = (sqsum[i] - lb_sq) - (total_pts * avg * avg);
    		if (utils::isEqual(ans,0))
        		ans = 0;
    		if (ans < 0) {
        		std::cout<<"ERROR in InteervalError "<<j<<i<<"..."<<sqsum[i]<<","<<lb_sq<<","<<sum[i]<<","<<lb_s<<","<<avg<<","<<ans<<std::endl;
				py::print("ERROR in IntervalError, got a negative error");
    		}
    		return ans;
		}

        void verify_memory(int arrsize) {
            double list_length;
    		if (utils::isEqual(this->epsilon_, 0)) list_length = arrsize;
    		else list_length = log(arrsize)/log(1 + epsilon_) + 1;
    		if (list_length > arrsize) list_length = arrsize; 
            if (hist_mem_ == NULL) {
			    hist_mem_  = new HistogramMemory(int(list_length), arrsize, n_buckets_);
                return;
            }
            
            bool sufficient_memory = true;
            if (arrsize < hist_mem_->arr_length_)
                sufficient_memory = false;
            if (n_buckets_ <  hist_mem_->n_buckets_)
                sufficient_memory  = false; 
            if (hist_mem_->max_length_ < int(list_length))
                sufficient_memory = false;
            if (!sufficient_memory) {
                delete hist_mem_;
                hist_mem_  = new HistogramMemory(int(list_length), arrsize, n_buckets_);
            }
        }
	
};


PYBIND11_MODULE(chistogram_new2, m) {
// optional module docstring
    m.doc() = "pybind11 plugin to compute best_split";
	py::class_<CHistogram_new2>(m, "chist")
        .def(py::init< double , int >())
        .def("get_epsilon", &CHistogram_new2::get_epsilon)
		.def("print_split", &CHistogram_new2::print_split)
		.def("get_bucket", &CHistogram_new2::get_bucket)
		.def("get_error", &CHistogram_new2::get_error)
		.def("get_bucket_list", &CHistogram_new2::get_bucket_list)
        .def("compute_histogram", &CHistogram_new2::compute_histogram)
        .def(py::pickle(
            [](const CHistogram_new2 &p) { // __getstate__
            /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(p.n_buckets_, p.epsilon_, reinterpret_cast<long long> (p.hist_mem_));
            },
            [](py::tuple t) { // __setstate__
            if (t.size() != 3)
                throw std::runtime_error("Invalid state!");

                /* Create a new C++ instance */
                CHistogram_new2 p(t[0].cast<double>(), t[1].cast<long long>());

                /* Assign any additional state */
                p.set_mem(t[2].cast<int>());

            return p;
        }
        ));
}

