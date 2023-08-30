
from typing import List
import numpy as np
from numpy.core import ndarray
import matplotlib.pyplot as plt
from node import Node
import h_wrapper as hwp
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing import get_context
to_print = False


#def create_node(k_args):
#    return Node(**k_args)

class Forest:
    """This creates a forest of trees, given the following list of parameters:
    1) n_trees: the number of trees
    2) max_depth: the depth of the trees
    3) max_samples: number of samples per tree
    4) max_buckets: maximum number of buckets used by the histogram
    5) epsilon: accuracy of the histogram"""
    dim = ...  # type: int
    size = ...  # type: int
    points = ...  # type: ndarray
    start = ...  # type: ndarray
    end = ...  # type: ndarray

    def __init__(self, **kwargs):
        self.n_trees = kwargs['n_trees']
        self.max_depth = kwargs['max_depth']
        self.max_samples = kwargs['max_samples']
        self.max_buckets = kwargs['max_buckets']
        self.epsilon = kwargs['epsilon']
        self.sample_axis = kwargs['sample_axis']
        self.threshold = kwargs['threshold']
        self.tree = {}
        # np.random.seed(seed = 17)
        self.n_leaves = np.zeros(self.n_trees)
        self.hist_array = [hwp.histWrap(self.max_buckets, self.epsilon) for i in range(self.n_trees)]


    def create_node(self, k_args, i):
        sample = np.random.choice(self.size, self.max_sample_size, replace=False)
        root = Node(**k_args)
        root.compute_density(sample)
        self.tree[i] = root
        self.n_leaves[i] = root.compute_leaf_num()

    def test_pool(self):
        print("in test_pool")
        #return x * x

    def fit(self, pts):
        self.points = pts
        self.dim, self.size = np.shape(self.points)
        if int(self.sample_axis*self.dim) == 0:
            print("sample_axis is too low")
            return
        self.start = np.zeros(self.dim)
        self.end = np.zeros(self.dim)
        for axis in range(self.dim):
            val = np.unique(np.array(self.points[axis]))
            if len(val) <= 1:
                print("No entropy in dimension :", axis)
                return
            self.start[axis] = (3 * val[0] - val[1]) / 2
            self.end[axis] = (3 * val[-1] - val[-2]) / 2
        k_args = {'depth': 0, 'forest': self}
        self.max_sample_size = np.min((self.size, self.max_depth * 200))
        #sample = np.random.choice(self.size, max_sample_size, replace=False)

        with get_context("spawn").Pool() as pool:
            print("in pool")
            pool.apply(self.test_pool)
            pool.join()
            #print("mapping")
            #print (pool.map(self.test_pool, range(10)))

        for i in range(self.n_trees):
            k_args['indices'] = np.random.choice(self.size, self.max_samples, replace=False)
            k_args['histogram'] = self.hist_array[i]
            self.create_node(k_args, i)
            #if __name__ == '__main__':
            #    print("name is main")
            #    proc = Process(target=self.create_node, args=(k_args, i, ))
            #    proc.start()
            #    print("process started",)
            #    proc.join()
            #    print("process joined")
            #    procs.append(proc)

            #root_node = Node(**k_args)


            #root_node = create_node(k_args)
            #root_node.compute_density(sample)
            #self.tree.append(root_node)
            #self.tree[i] = root_node
            #self.n_leaves[i] = root_node.compute_leaf_num()
        #for proc in procs:
        #    print("process joining")
        #    proc.join()

    def plt_scores(self, pts):
        _, n_pts = np.shape(pts)
        n_show = int(2 * self.n_trees / 3)
        scores = np.zeros((self.n_trees, n_pts))  # need to do this streaming
        indices = [i for i in range(n_pts)]
        for i in range(self.n_trees):
            self.tree[i].compute_split(pts, indices, scores[i])
        for i in range(n_pts):
            plt.plot(np.sort(scores[:, i])[:n_show])
            plt.show()

    def predict(self, pts, err=0.1, pct=50):
        _, n_pts = np.shape(pts)
        scores = np.zeros((self.n_trees, n_pts))  # need to do this streaming
        indices = [i for i in range(n_pts)]
        for i in range(self.n_trees):
            self.tree[i].compute_split(pts, indices, scores[i])
        n_err = int(err * n_pts)
        min_score = np.percentile(scores, pct, axis=0)
        top_indices = np.argsort(min_score)[:n_err]
        anom_pts = {}
        anom_scores = {}
        anom_pct = {}
        for i in range(n_err):
            anom_pts[top_indices[i]] = pts[:, top_indices[i]]
            anom_scores[top_indices[i]] = scores[:, top_indices[i]]
            anom_pct[top_indices[i]] = min_score[top_indices[i]]
        return top_indices, anom_pts, anom_scores, anom_pct, min_score
    
    def predict_depth(self, pts, err=0.1, pct=50):
        _, n_pts = np.shape(pts)
        scores = np.zeros((self.n_trees, n_pts))  # need to do this streaming
        indices = [i for i in range(n_pts)]
        for i in range(self.n_trees):
            self.tree[i].compute_depth(pts, indices, scores[i])
        n_err = int(err * n_pts)
        min_score = np.percentile(scores, pct, axis=0)
        top_indices = np.argsort(min_score)[:n_err]
        anom_pts = {}
        anom_scores = {}
        anom_pct = {}
        for i in range(n_err):
            anom_pts[top_indices[i]] = pts[:, top_indices[i]]
            anom_scores[top_indices[i]] = scores[:, top_indices[i]]
            anom_pct[top_indices[i]] = min_score[top_indices[i]]
        return top_indices, anom_pts, anom_scores, anom_pct, min_score
    
    def find_min_cube(self, pts, n_err=1, pct1=50, pct2=50):
        _, n_pts = np.shape(pts)
        scores = np.zeros((self.n_trees, n_pts))  # need to do this streaming
        indices = [i for i in range(n_pts)]
        for i in range(self.n_trees):
            self.tree[i].compute_split(pts, indices, scores[i])
        min_score = np.percentile(scores, pct1, axis=0)
        top_indices = np.argsort(min_score)[:n_err]
        sort_trees = np.argsort(scores, axis=0)
        min_tree = sort_trees[int(self.n_trees*pct2/100.0)]
        boxes = []
        for i in range(0,n_err):
            print('----')
            best_tree = min_tree[top_indices[i]]
            boxes.append( self.tree[best_tree].compute_box(pts[:,top_indices[i]]) )
            # boxes = self.tree[best_tree].compute_box(pts[:,top_indices[i]])
        return top_indices, boxes








