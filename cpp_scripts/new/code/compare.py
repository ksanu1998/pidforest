import numpy as np
from sklearn import metrics
# Specify the file path
lab_file_path = "y_slice.txt"
sco_file_path = "our_scores.txt"
# Read the column of values into a NumPy array
y_slice = np.loadtxt(lab_file_path, delimiter='\n')
our_scores = np.loadtxt(sco_file_path, delimiter='\n')

precision_our, recall_our, thresholds_our = metrics.precision_recall_curve(y_slice, our_scores, pos_label=1)
print(max(2*precision_our*recall_our/(precision_our+recall_our)))