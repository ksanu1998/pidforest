import numpy as np
import cforest as cfr

arr_size = 3000
Arr = np.random.rand(arr_size)
Arr = np.reshape(Arr, (10,300))
kwargs = {'max_depth': 12, 'n_trees':20,  'max_samples': 300, 'max_buckets': 3, 'epsilon': 0.1, 'sample_axis': 1, 'threshold': 0}

cforest = cfr.Forest(**kwargs)
cforest.fit(Arr)