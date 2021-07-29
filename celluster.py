import argparse
import pandas as pd
from collections import Counter
import numpy as np


'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Implementation of the consensus clustering method proposed in https://www.sciencedirect.com/science/article/pii/S131915781930597X.')
    parser.add_argument('-i', '--input', help='Input cell cluster assignment files.', nargs='*', action='store', dest='input', required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved.', type=str, required=False)
    parser.add_argument('-c', '--id', help='The name of the column that contains the item ID.', type=str, required=True)
    parser.add_argument('-d', '--data', help='The orignal file that was used for clustering.', type=str, required=False)
    parser.add_argument('-t', '--tab', help='Flag to indicate that input files are tab delimited.', action='store_true', required=False)
    parser.add_argument('-v', '--verbose', help='Flag to print information about the consensus clustering.', action='store_true', required=False)
    parser.add_argument('-z', '--outliers', help='Flag to keep outliers in output file under the cluster label: 0.', action='store_true', required=False)
    parser.add_argument('-k', '--clusters', help='The file that contains cluster labels from previous iterations.', type=str, required=False)
    parser.add_argument('-n', '--name', help='The name of this dataset for output file naming consistency.', type=str, required=False)
    parser.add_argument('-r', '--iter', help='The iteration number.', type=str, required=False)
    args = parser.parse_args()
    return args


'''
Get input data file name
'''
def getDataName(path):
    fileName = path.split('/')[-1] # get filename from end of input path
    dataName = fileName[:fileName.rfind('.')] # get data name by removing extension from file name
    return dataName


'''
Get the header from an input csv file.
Returns a list where each element is a column header.
'''
def getHeader(file):
    # get first line of file aka the csv header
    with open(file) as f:
        header = f.readline().rstrip()
    f.close()

    # get the input file column headers as list
    header_columns = header.split(delimiter)

    return header_columns


'''
Look at first line of an in input file (the header of the csv) to assess if it is the correct format.
'''
def validInput(file):

    # definition of the columns needed to constitute a valid input file
    NEEDED_COLUMNS = [id, CLUSTER, METHOD]

    input_columns = getHeader(file) # get input file header columns
    valid = all(col in input_columns for col in NEEDED_COLUMNS) # check if all needed columns are present in input columns list
    
    return valid


'''
Read file into dataframe.
'''
def getData(file):

    data = pd.read_csv(file, delimiter=delimiter, index_col=id) # load data from input csv into dataframe
    method = data[METHOD].iloc[0] # get method name from first row (this is assuming at all rows are from the same method) NOTE: we may not want to assume this!!
    data = data.drop(METHOD, axis=1) # drop method column
    
    return data, method


'''
Compute the jaccard index of two clusters (sets of samples).
'''
def jaccard_index(cluster1, cluster2):
    i = len(cluster1.intersection(cluster2))
    u = len(cluster1.union(cluster2))
    jaccard_index = i / u
    return jaccard_index


'''
Get the method that has the fewest (min) clusters. If there is a tie, just use the method whose column comes first
'''
def getMinMethod(cluster_table):
    min = len(cluster_table) # start with the min at the number of cells
    min_method = '' # the method with the fewest clusters

    # iterate through methods 
    for method in cluster_table.columns:
        num_clusters = len(pd.unique(cluster_table[method])) # calucate the number of clusters
        if num_clusters < min: # if min, update values
            min_method = method
            min = num_clusters

    return min_method


'''
Get clustering method re-labeling key
'''
def getNewLabels(min_clusters_table, m_clusters_table):
    label_key = {} # a dict where the key is the cluster label in the 'm' method and the value is the cluster label in the 'min' method that has the highest jaccard index

    # lists of cluster labels in each method
    m_clusters = pd.unique(m_clusters_table)
    min_clusters = pd.unique(min_clusters_table)

    for c_label1 in m_clusters:
        indices = {}
        for c_label2 in min_clusters:
            c_label1_items = set(m_clusters_table[m_clusters_table == c_label1].index)
            c_label2_items = set(min_clusters_table[min_clusters_table == c_label2].index)
            i = jaccard_index(c_label1_items, c_label2_items)
            indices[i] = c_label2
        max_index = max(indices.keys()) # find the largest jaccard index
        label_match =  indices[max_index] # get the label in the min method with that max jaccard index
        label_key[c_label1] = label_match # add to key dict
    
    return label_key


'''
Get the cluster label that has the majority vote
'''
def vote(labels):
    c = Counter(labels) # count the most common label

    # check for majority
    label, count = c.most_common()[0] # get most common label and its count
    freq = count / len(labels) # the frequency of the most common label

    if freq > 0.5:
        return label # return mode label if it has majority
    else:
        return np.nan # if it is a tie return NaN


'''
Consensus cluster by matching clusters between methods and relabeling them appropriately. Items without majority vote for cluster label are outliers.
'''
def consensusCluster(cluster_table):

    # compare all methods against the one with the least number of clusters
    min_method = getMinMethod(cluster_table) # get the method with the fewest number of clusters to compare all others against
    methods = list(cluster_table.columns) # get list of methods
    methods.remove(min_method) # remove min method because we will not compare it against itself

    for m in methods:
        label_key = getNewLabels(cluster_table[min_method], cluster_table[m]) # a dict where the key is the cluster label in the 'm' method and the value is the cluster label in the 'min' method that has the highest jaccard index
        cluster_table[m] = cluster_table[m].map(label_key).fillna(cluster_table[m]) # relabel clusters FASTER

    # vote
    consensus_labels = [vote(row) for row in cluster_table.to_numpy()] # a list of clusters labels (or nan) ordered for each item ID in cluster_table
    cluster_table[CLUSTER] = consensus_labels # add results to table


'''
If a previous consensus cluster file is provided, update it to contain the latest iteration item cluster labels
'''
def updateConsensusClusters(file):
    prev_clusters = pd.read_csv(file, delimiter=delimiter, index_col=id) # get previous consensus cluster labels for items
    max = prev_clusters[CLUSTER].max() # get max cluster label
    min = cluster_table[CLUSTER].min() # get min cluster label from latest iteration

    # increment latest cluster labels appropriately (might start at 0 or 1?)
    offset = abs(min - max) # find the difference between the last max and the latest min
    increment = offset + 1 # set increment so that the min cluster label in the latest iteration becomes the next label from the previous ones
    cluster_table[CLUSTER] = cluster_table[CLUSTER] + increment # increment latest cluster labels
    if args.verbose:
        print(f'This iteration cluster labels were incremented by {increment}')

    clusters_to_add = cluster_table[CLUSTER].dropna().to_frame()

    # update previous cluster labels table
    # if latest items are in prev_clusters, update cluster label, otherwise append rows
    contains_all = all(elem in list(prev_clusters.index) for elem in list(cluster_table.index))
    if contains_all:
        prev_clusters.update(clusters_to_add)
    else:
        constains_any = any(item in list(prev_clusters.index) for item in list(cluster_table.index))
        if not constains_any:
            prev_clusters = prev_clusters.append(clusters_to_add).sort_index()
        else:
            print('ERROR. I have not coded this yet.')

    # write cluster labels as ints
    prev_clusters[CLUSTER] = prev_clusters[CLUSTER].astype(int)

    return prev_clusters
    

'''
Get and write the input data of outlier items that need to be re-consensus clustered in the next iteration.
'''
def writeOutlierData(recluster_ids):
    # get items that need to be reclustered
    if args.data != None and len(recluster_ids) > 0:
        data = pd.read_csv(args.data, delimiter=delimiter, index_col=id)
        recluster_data = data.loc[recluster_ids]

    # write data to recluster
    if args.data != None and len(recluster_ids) > 0:
        recluster_data.to_csv(f'{output}/{data_prefix}-iter-outliers.{extension}', sep=delimiter)
        if args.verbose:
            print(f'Items to recluster are in the file: {output}/{data_prefix}-iter-outliers.{extension}')


'''
Format consensus clustered items for output
'''
def writeConsensusClusters(cluster_table):
    # write item cluster labels
    if args.outliers:
        cluster_table[CLUSTER] = cluster_table[CLUSTER] + 1 # increment cluster labels
        cluster_table[CLUSTER] = cluster_table[CLUSTER].map({np.nan : 0}).fillna(cluster_table[CLUSTER]) # set outliers to label '0'
        if args.verbose:
            print('Outlier items have the cluster label: 0')
    else:
        cluster_table = cluster_table.dropna()

    cluster_table[CLUSTER] = cluster_table[CLUSTER].astype(int)
    cluster_table.to_csv(f'{output}/{cluster_prefix}-iter-consensus.{extension}', sep=delimiter, columns=[CLUSTER])
    if args.verbose:
        print(f'Consensus cluster labels are in the file: {output}/{cluster_prefix}-iter-consensus.{extension}')


'''
Main.
'''
if __name__ == '__main__':

    # constants
    CLUSTER = 'Cluster' # the header of the cluster assignmnet column
    METHOD = 'Method' # the header of the method column

    # parse arguments
    args = parseArgs()

    id = args.id # the header of the sample/cell/drug/item ID column

    # whether to read and write csv or tsv
    if args.tab:
        delimiter = '\t'
        extension = 'tsv'
    else:
        delimiter = ','
        extension = 'csv'

    # get iteration number if provided
    if args.iter != None:
        iter = '-' + args.iter

    # get prefix for output files
    if args.name == None:
        cluster_prefix = getDataName(args.input[0]) + iter
    else:
        cluster_prefix = args.name + iter
    if args.data != None:
        if args.name == None:
            data_prefix = getDataName(args.data) + iter
        else:
            data_prefix = args.name + iter

    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output

    if args.verbose:
        print('Reading files...')

    # if input file has correct format, make a matrix from it
    cluster_table = pd.DataFrame() # a df where cell ID is the index and each column is the cluster assignment from a different algorithm
    for file in args.input:
        if validInput(file): # check that input file has the needed columns
            data, method = getData(file)
            cluster_table[method] = data[CLUSTER]
        else:
            print(f'{file} is incorrectly formatted.')

    if args.verbose:
        print('Performing consensus clustering...')

    # perform consensus clustering
    consensusCluster(cluster_table)

    if args.verbose:
        print('Retrieving outliers...')

    # get item ID's that need to be re-clustered in the next iteration
    recluster_ids = cluster_table[cluster_table[CLUSTER].isnull()].index

    # update running consensus clusters file
    if args.clusters != None:
        complete_cluster_table = updateConsensusClusters(args.clusters)
        complete_cluster_table.to_csv(f'{output}/{cluster_prefix}-running-consensus.{extension}', sep=delimiter, columns=[CLUSTER])
        if args.verbose:
            print(f'Wrote running consensus clusters to: {output}/{cluster_prefix}-running-consensus.{extension}')
    
    if args.verbose:
        print('Writing output files...')

    writeOutlierData(recluster_ids) # write outlier data
    writeConsensusClusters(cluster_table) # write consensus cluster labels

    # final stats
    if args.verbose:
        print('Done.')
        print(f'{len(cluster_table.index)} ({(len(cluster_table.index) / (len(recluster_ids) + len(cluster_table.index))) * 100:.2f}%) items clustered into {len(pd.unique(cluster_table[CLUSTER]))} clusters.')
        print(f'{len(recluster_ids)} ({(len(recluster_ids) / (len(recluster_ids) + len(cluster_table.index))) * 100:.2f}%) outlier items to be re-clustered in next iteration.')