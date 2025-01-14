import logging
import numpy as np
import matplotlib.pyplot as plt
from scripts.forest import Forest
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
import scripts.timeseries as ts
import pandas as pd
import time
import scipy.io as sio
from sklearn import metrics
import os

formatter = logging.Formatter('%(message)s')

@profile
def setup_logger(name, log_file, level=logging.INFO):

    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

    
# datasets = ['nyc_taxi','ambient_temperature_system_failure','cpu_utilization_asg_misconfiguration','machine_temperature_system_failure','ec2_request_latency_system_failure']
datasets = ['nyc_taxi','ambient_temperature_system_failure','cpu_utilization_asg_misconfiguration','machine_temperature_system_failure'] # original
# datasets = ['nyc_taxi']
L = len(datasets)
trials = 5 # original
# trials = 1
run_lof_svm = 1

file_prefix = 'experiment_results/'
if not os.path.exists(file_prefix):
    os.mkdir(file_prefix)
file_prefix = os.path.join(file_prefix, 'time_logs')
if not os.path.exists(file_prefix):
    os.mkdir(file_prefix)

logger_e2e = setup_logger("logger_e2e", os.path.join(file_prefix, "e2e.csv"))
reference_time = time.time()
for i in range(0,L):
    
    logger_path = os.path.join(file_prefix, datasets[i])
    if not os.path.exists(logger_path):
        os.mkdir(logger_path)
    logger_isoforest = setup_logger("logger_isoforest_" + datasets[i], os.path.join(logger_path, "isoforest.csv"))
    logger_lof = setup_logger("logger_lof_" + datasets[i], os.path.join(logger_path, "lof.csv"))
    logger_onecsvm = setup_logger("logger_onecsvm_" + datasets[i], os.path.join(logger_path, "onecsvm.csv"))
    logger_pidforest = setup_logger("logger_pidforest_" + datasets[i], os.path.join(logger_path, "pidforest.csv"))
       

    data = pd.read_csv('../data/numenta/'+datasets[i]+'.csv')
    value = list(data["value"])
    arr = np.array(value)
    y = np.array(list(data["label"]))
    X = ts.shingle(arr, 10)
    X = np.transpose(X)
    t1, _ = np.shape(X)
    if not os.path.exists("./experiment_results"):
        os.mkdir("./experiment_results")
    file_name = './experiment_results/' + datasets[i] + '.txt'
    File_object = open(file_name,"w")   
    time_all = np.zeros((trials,4))
    precision_all = np.zeros((trials,4))
    auc_all = np.zeros((trials,4))

    for j in range(0,trials):
    
        print('\n\n******'+datasets[i]+' trial '+str(j+1)+'*******\n\n')
        
        isoforest_start = time.time()

        print('\n******Iso-Forest*******\n')
        start = time.time()
        # clf = IsolationForest(contamination = 0.1, behaviour = 'new') # behaviour arg DEPRECATED
        clf = IsolationForest(contamination = 0.1)
        clf.fit(X)
        end = time.time()
        time_all[j,0] = end - start
        iso_scores = clf.score_samples(X)
        
        isoforest_end = time.time()
        logger_isoforest.info(str(j) + "," + str(isoforest_end - isoforest_start))

        if run_lof_svm == 0:
            lof_scores = iso_scores
            osvm_scores = iso_scores
        elif j == 0:
            print('\n******LOF*******\n')

            lof_start = time.time()

            start = time.time()
            lof = LocalOutlierFactor()
            lof.fit(X)
            end = time.time()
            time_all[j,1] = end - start
            lof_scores = lof.negative_outlier_factor_

            lof_end = time.time()
            logger_lof.info(str(j) + "," + str(lof_end - lof_start))

            print('\n******1-class SVM*******\n')

            onecsvm_start = time.time()

            start = time.time()
            osvm = OneClassSVM(kernel='rbf')
            osvm.fit(X)
            end = time.time()
            time_all[j,2] = end - start
            osvm_scores = osvm.score_samples(X)
                
            onecsvm_end = time.time()
            logger_onecsvm.info(str(j) + "," + str(onecsvm_end - onecsvm_start))

        print('\n******Our Algo*******\n')

        pidforest_start = time.time()

        start = time.time()
        #n_samples = int(t1/50)
        n_samples = 100 # original
        kwargs = {'max_depth': 10, 'n_trees':50,  'max_samples': n_samples, 'max_buckets': 3, 'epsilon': 0.1, 'sample_axis': 1, 
          'threshold': 0} # original
        forest = Forest(**kwargs)
        forest.fit(np.transpose(X))
        indices, outliers, scores, pst, our_scores = forest.predict(np.transpose(X), err = 0.1, pct = 0)
        end = time.time()
        time_all[j,3] = end - start
        
        pidforest_end = time.time()
        logger_pidforest.info(str(j) + "," + str(pidforest_end - pidforest_start))

        precision_iso, recall_iso, thresholds_iso = metrics.precision_recall_curve(y[:t1], -iso_scores, pos_label=1)
        precision_lof, recall_lof, thresholds_lof = metrics.precision_recall_curve(y[:t1], -lof_scores, pos_label=1)
        precision_osvm, recall_osvm, thresholds_osvm = metrics.precision_recall_curve(y[:t1], -osvm_scores, pos_label=1)
        precision_our, recall_our, thresholds_our = metrics.precision_recall_curve(y[:t1], -our_scores, pos_label=1)
        precision_all[j,0] = max(2*precision_iso*recall_iso/(precision_iso+recall_iso))
        precision_all[j,1] = max(2*precision_lof*recall_lof/(precision_lof+recall_lof))
        precision_all[j,2] = max(2*precision_osvm*recall_osvm/(precision_osvm+recall_osvm))
        precision_all[j,3] = max(2*precision_our*recall_our/(precision_our+recall_our))
        
        auc_all[j,0] = metrics.roc_auc_score(y[:t1], -iso_scores)
        auc_all[j,1] = metrics.roc_auc_score(y[:t1], -lof_scores)
        auc_all[j,2] = metrics.roc_auc_score(y[:t1], -osvm_scores)
        auc_all[j,3] = metrics.roc_auc_score(y[:t1], -our_scores)
        
        for k in range(0,4):
            print('{:.4f}\t'.format( precision_all[j,k] ))
        print('\n')
        
        for k in range(0,4):
            print('{:.4f}\t'.format( auc_all[j,k] ))
        print('\n')
    

    File_object.write(str(kwargs))    

    File_object.write('\n\nIF\tLOF\tSVM\tOur-Algo\n\n')    
        
    for j in range(0,trials):
        for k in range(0,4):
            File_object.write('{:.4f}\t'.format( precision_all[j,k] ))
        File_object.write('\n')
        
    File_object.write('\n')    
        
    for k in range(0,4):
            File_object.write('{:.4f}\t'.format( np.mean(precision_all[:,k]) ))
    File_object.write('\n')    
    
    for k in range(0,4):
            File_object.write('{:.4f}\t'.format( np.std(precision_all[:,k]) ))
    File_object.write('\n')
    
    File_object.write('\nIF\tLOF\tSVM\tOur-Algo\n\n')    
        
    for j in range(0,trials):
        for k in range(0,4):
            File_object.write('{:.4f}\t'.format( auc_all[j,k] ))
        File_object.write('\n')
        
    File_object.write('\n')    
        
    for k in range(0,4):
            File_object.write('{:.4f}\t'.format( np.mean(auc_all[:,k]) ))
    File_object.write('\n')    
    
    for k in range(0,4):
            File_object.write('{:.4f}\t'.format( np.std(auc_all[:,k]) ))
    File_object.write('\n')
        
    File_object.close()

    file_name = 'experiment_results/' + datasets[i] + '_results.mat'
    sio.savemat(file_name, {'time_all':time_all, 'precision_all':precision_all, 'auc_all':auc_all})

logger_e2e.info(time.time() - reference_time) 