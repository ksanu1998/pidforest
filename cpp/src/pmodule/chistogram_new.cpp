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

class CHistogram_new {
    public:
        CHistogram_new( double epsilon, int bucket) : epsilon_(epsilon), n_buckets_(bucket) {
            hist_mem_ = NULL; 
            if (n_buckets_ > MAX_BUCKETS) {
				py::print("number of buckets reduced to ", MAX_BUCKETS);
				n_buckets_ = MAX_BUCKETS;
			}
        }

        ~CHistogram_new() { 
            if (hist_mem_ != NULL)
                delete hist_mem_;
        }

        /**
         * Computes the histogram and fills the Split data member
         */
        void compute_histogram( py::array parr, py::array pcounts) {
            try {
                double * arr = extract_array_ptr(parr);
                double * counts = extract_array_ptr(pcounts);
                int arr_size = extract_array_size(parr);
                verify_memory(arr_size);
                double * total_counts;
                int  app_buckets_list[MAX_BUCKETS][MAX_BUCKETS];
		        double  app_buckets_err[MAX_BUCKETS]; //an array with the app error accumulated per number of buckets used
                total_counts = hist_mem_->total_count_;
			    total_counts[0] = 0;
			    for (int k = 0; k < arr_size; k++)
				    total_counts[k+1] = total_counts[k] + counts[k];
                

                double * sum = hist_mem_->sum_;
    		    double * sq_sum = hist_mem_->sqsum_;
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
        		    app_buckets_list[k][0] = 0;
    		    }
                sum[0] = arr[0] * counts[0];
    		    sq_sum[0] = arr[0] * arr[0] * counts[0];
                double tmp = 0;

                /***** MAIN LOOP *****/
                for (int j = 1; j < arr_size; j++) {
        		    //preparing the one bucket case. That's the recursion base.
        		    sum[j] = sum[j-1] + (arr[j] * counts[j]);
        		    sq_sum[j] = sq_sum[j - 1] + (arr[j] * arr[j] * counts[j]);
        		    tmp = sq_sum[j] - ((sum[j] * sum[j]) / (double)total_counts[j + 1]); //error for 1 bucket
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

                    for (int k = 1; k < n_buckets_; k++) {
            		    double tmp_err = sq_sum[j];
            		    int curr_index = 0;
            		    for (int i = 0; i <= list_length[k - 1]; i++) {
                		    int b  = b_val[k-1][i];
                            double can_err = app_err_b[k - 1][i] + IntervalError(sum, sq_sum, counts, b + 1, j);
                		    if (tmp_err > can_err) {
                    		    tmp_err = can_err;
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
                } // end main loop

                for (int k = 0; k < n_buckets_; k++) {
        		    app_buckets_err[k] = app_err_b[k][list_length[k]];
   	 		    }

                // Now find the best split:
   			    double var_red = 0;
   			    int best_bucket = 0;
  			    if (utils::isEqual(app_buckets_err[0], 0)) {
      			    split_.bucket_error_ = 0;
      			    split_.bucket_level_ = 0;
  			    }
  			    else {
      			    for (int k = 1; k < n_buckets_; k++) {
          			    double err_red = (app_buckets_err[0] - app_buckets_err[k]) / app_buckets_err[0];
          			    if (err_red > var_red) {
              			    var_red = err_red;
              			    best_bucket = k;
          			    }
      			    }
          		    split_.bucket_error_ = var_red;
          		    split_.bucket_level_ = best_bucket;
  			    }
                ComputeAppList(bucket_indices, b_val, list_length[best_bucket], best_bucket, split_.bucket_list_);	
            } catch (const std::runtime_error& error) {
				py::print(error.what());
   			} 
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

		/**
		 * Gets a numpy array as argument and returns the list of buckets by changing its values in place
		 */
		void get_bucket_list(py::array b_list) {
			py::buffer_info arr_info = b_list.request();
			double * list_ptr = static_cast<double *>(arr_info.ptr);
			for (int k = 0; k < split_.bucket_level_ + 1; k++)
        		list_ptr[k] = split_.bucket_list_[k];  
		}

    private:
        double epsilon_;
        int n_buckets_;
        HistogramMemory * hist_mem_;
        Split split_;

        double * extract_array_ptr( py::array arr ) {
            py::buffer_info arr_info = arr.request();
            double * arr_ptr;
			arr_ptr = static_cast<double *>(arr_info.ptr);
			if (arr_info.format.compare("d") != 0) {
				throw std::runtime_error("Exception: wrong format for data array");
			}
            return arr_ptr;
        }

        int extract_array_size( py::array arr ) {
            py::buffer_info arr_info = arr.request();
            if (arr_info.ndim > 1) {
				throw std::runtime_error("Exception: array has to be one dimensional");
			}
			int arr_size = 1;
    			for (auto r: arr_info.shape) {
      			arr_size *= r;
    		}
            return arr_size;
        }

        // computes the error of a bucket  in indexes j,...,i
		double IntervalError(double *sum, double *sqsum, double * total_counts, int j, int i) {
    		double lb_s, lb_sq;
    		if (i <= j) { return 0; }
    		if (j == 0) { 
				lb_s = 0; lb_sq = 0;
			}
        	else {
            	lb_s = sum[j-1]; 
            	lb_sq = sqsum[j-1]; 
        	}
			double total_pts = total_counts[i+1] - total_counts[j];   
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

        //computes the bucket boundaries
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

        // verifies hist_mem_ has sufficient memory. Allocates more if needed.
        void verify_memory(int arrsize) {
            double lb; //compute list length
            if (utils::isEqual(this->epsilon_, 0)) lb = arrsize;
            else lb = log(arrsize)/log(1 + this->epsilon_) + 1;
            if (lb > arrsize) lb = arrsize; 

            if (hist_mem_ == NULL) {
                hist_mem_ = new HistogramMemory((int) lb, arrsize, n_buckets_);
                return;
            }

            bool sufficient_memory = true;
            if (arrsize < hist_mem_->arr_length_)
                sufficient_memory = false;
            if (this->n_buckets_ <  hist_mem_->n_buckets_)
                sufficient_memory  = false;
            if (hist_mem_->max_length_ < int(lb))
                sufficient_memory = false;
            
            if  (!sufficient_memory) { //alocate new memory
                delete hist_mem_;
                hist_mem_ = new HistogramMemory((int) lb, arrsize, n_buckets_);
            }
        }   
};

PYBIND11_MODULE(chistogram_new, m) {
// optional module docstring
    m.doc() = "pybind11 plugin to compute best_split";
	py::class_<CHistogram_new>(m, "chist")
        .def(py::init<double , int>())
        .def("get_epsilon", &CHistogram_new::get_epsilon)
		.def("print_split", &CHistogram_new::print_split)
		.def("get_bucket", &CHistogram_new::get_bucket)
		.def("get_error", &CHistogram_new::get_error)
		.def("get_bucket_list", &CHistogram_new::get_bucket_list)
        .def("compute_histogram", &CHistogram_new::compute_histogram);
}