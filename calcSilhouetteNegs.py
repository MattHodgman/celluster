import pandas as pd
from sklearn.metrics import silhouette_samples
import numpy as np
import sys

# params: (1) clean, transformed feature table (2) cell cluster labels

# constants
CELLID='CellID'
CLUSTER='Cluster'

# get data
X = pd.read_csv(sys.argv[1], index_col=CELLID)
cluster_labels = pd.read_csv(sys.argv[2], index_col=CELLID)

# calculate silhouette values
sample_silhouette_values = silhouette_samples(X, cluster_labels[CLUSTER])

# calculate proportion of negative silhouette values
n_all = len(sample_silhouette_values)
n_neg = np.sum(sample_silhouette_values < 0, axis=0)
prop = n_neg / n_all

print(prop)