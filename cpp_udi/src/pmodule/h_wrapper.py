import numpy as np
import chistogram_new2 as chs


def testy(a):
    print(a)


class histWrap:

    def __init__(self, max_buckets, eps):
        try:
            self.max_buckets = max_buckets
            self.chist = chs.chist(eps, max_buckets)
        except():
            print("Exception")

    def compute_histogram(self, val, counts, A):
        counts = np.array(counts, dtype='double')
        return self.chist.compute_histogram(val, counts, A)

    def print_split(self):
        self.chist.print_split()

    def get_epsilon(self):
        return self.chist.get_epsilon()

    def get_bucket(self):
        return self.chist.get_bucket()

    def get_error(self):
        return self.chist.get_error()

    def get_bucket_list(self):
        A = np.ones(self.get_bucket())
        A = np.array(A, dtype='double')
        self.chist.get_bucket_list(A)
        A = A + 1
        A = np.array(A, dtype='int64')
        return A[:-1]

    #last two items of A are bucket and error.
    def best_split(self, val, counts):
        A = np.zeros(self.max_buckets + 1)
        bucket = self.compute_histogram(val, counts, A)
        the_error = A[-1]
        B = A[:bucket-1]
        B = np.array(B+1, dtype = 'int64')
        return bucket, the_error, B

