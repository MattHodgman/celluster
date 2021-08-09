# celluster
Consensus clustering from multiple methods used in [MCMICRO](https://mcmicro.org/) to cluster cell types. An implementation of the method described in [this paper](https://www.sciencedirect.com/science/article/pii/S131915781930597X). Celluster works by taking the clustering results from different algorithms or algorithms run with different parameters, identifying the algorithm `A` with the fewest clusters and then mapping each of the clusters from the other algorithms to one of the clusters in `A`, creating essentially a many-to-one mapping of cluster labels. Any outlier cells can be put into the cluster with label `-1`. Celluster outputs the consensus clustering results as well as can output the feature table (from the input `feature_table.csv`) that contains only the outlier cells, so as to be used in the next iteration of clustering by the individual algorithms and then celluster.

Example usage (Python3):
```
python3 celluster.py -i clusters1.csv clusters2.csv clusters3.csv -o . -c 'CellID' -d feature_table.csv -z
```

Example usage (NextFlow + Docker):
```
docker build -t celluster .
nextflow run main.nf --input unmicst-exemplar-001.csv
```

## input file format
Input files such as `clusters1.csv` in the example usage have three columns: `CellID,Cluster,Method` where `CellID` is a unique integer identifier for each cell, `Cluster` is a unique cluster label identifier (an integer), and `Method` is short title of the method used to obtain the clustering labels. `feature_table.csv` is the original feature table used as input for the original clustering methods. 

## Parameters
```
usage: celluster.py [-h] -i [INPUT [INPUT ...]] [-o OUTPUT] -c ID [-d DATA] [-t] [-v] [-z] [-k CLUSTERS] [-n NAME] [-r ITER]

Implementation of the consensus clustering method proposed in https://www.sciencedirect.com/science/article/pii/S131915781930597X.

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT [INPUT ...]], --input [INPUT [INPUT ...]]
                        Input cell cluster assignment files.
  -o OUTPUT, --output OUTPUT
                        The directory to which output files will be saved.
  -c ID, --id ID        The name of the column that contains the cell ID.
  -d DATA, --data DATA  The orignal file that was used for clustering.
  -t, --tab             Flag to indicate that input files are tab delimited.
  -v, --verbose         Flag to print information about the consensus clustering.
  -z, --outliers        Flag to keep outliers in output file under the cluster label: -1.
  -k CLUSTERS, --clusters CLUSTERS
                        The file that contains cluster labels from previous iterations.
  -n NAME, --name NAME  The name of this dataset for output file naming consistency.
  -r ITER, --iter ITER  The iteration number.
  ```