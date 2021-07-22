import argparse
import numpy as np
import pandas as pd
from itertools import combinations
from tqdm import tqdm
from scipy import sparse
import csv
from operator import itemgetter


'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Consensus clustering from multiple methods used in MCMICRO to cluster cell types.')
    parser.add_argument('-i', '--input', help='Input cell cluster assignment files.', nargs='*', action="store", dest="input", required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved.', type=str, required=False)
    args = parser.parse_args()
    return args

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
    header_columns = header.split(',')

    return header_columns


'''
Look at first line of an in input file (the header of the csv) to assess if it is the correct format.
'''
def validInput(file):

    # definition of the columns needed to constitute a valid input file
    NEEDED_COLUMNS = [CELLID, CLUSTER, METHOD]

    input_columns = getHeader(file) # get input file header columns
    valid = all(col in input_columns for col in NEEDED_COLUMNS) # check if all needed columns are present in input columns list
    
    return valid


'''
Read file into dataframe.
'''
def getData(file):

    data = pd.read_csv(file, delimiter=',', index_col=CELLID) # load data from input csv into dataframe
    method = data[METHOD].iloc[0] # get method name from first row (this is assuming at all rows are from the same method) NOTE: we may not want to assume this!!
    data = data.drop(METHOD, axis=1) # drop method column
    
    return data, method


'''
Reformat consensus clusters dict to a list of tuples where each tuple is a cell-cluster assignment.
'''
def formatConsensusClusters(consensus_clusters):

    ccc_assignments = [] # short for cell consensus cluster assignments
    cluster_assignment = 0 # arbitrary int to identify different clusters

    # iterate through dict
    for cluster,cells in consensus_clusters.items():
        for cell in cells:
            ccc_assignments.append((cell,cluster_assignment))
        cluster_assignment += 1

    # sort
    ccc_assignments = sorted(ccc_assignments, key=itemgetter(0))

    return ccc_assignments


'''
Get cluster assignments that have 100% agreement between methods.
Return a dict where the key is a tuple with the cluster assignments for each method and the value is a list of cell IDs
'''
def getConsensusClusters(cluster_table):
    cluster_table = cluster_table[cluster_table.duplicated(keep=False)] # drop rows that aren't duplicates (duplicate rows mean that there are at least 2 cells that have were assigned to the same clusters across all algorithms)
    consensus_clusters = cluster_table.groupby(list(cluster_table.columns), as_index=False).groups # dictionary where the key is the cluster assignments of all datasets and the value is a list of all the cells that received those assignments

    return consensus_clusters


'''
Make a dictionary where the key is a tuple with two cell IDs and the value is number of clustering methods that assign them to the same cluster
'''
def makeConsensusDict(cluster_table):

    # there may be a way to write off certain pairs as outliers while we are filling up the 'matrix'
    # if a pair hasn't had any agreement from other algorithms and we are already through most of the algorithms, we could just write it off as an outlier

    # here's another idea
    # literally just get the duplicate rows and their indices are the cells that have 100% agreement of clusters?

    cluster_table = cluster_table[cluster_table.duplicated(keep=False)] # drop rows that aren't duplicates
    # (duplicate rows mean that there are at least 2 cells that have were assigned to the same clusters across all algorithms)

    methods = list(cluster_table.columns) # get the method (column) names

    clusters_all = cluster_table.groupby(methods, as_index=False).groups # dictionary where the key is the cluster assignments of all datasets and the value is a list of all the cells that received those assignments

    return clusters_all

    # ok so there are like 250 'clusters' where 'all' datasets agree these cells belong together
    # but like each dataset has 30 clusters??? 30 * 30 * 25

    # i could group clusters together if they have a 2/3 match of cluster assignments (how do i adapt that threshold for variable number of datasets? >50%?)
    # but will this just merge all of them together into one? hahaha lets just try it ig hahah

    # clusters_partial = {} # key is ???? value is all cells in that partially supported cluster

    # what if we sorted the clusters by size, then going from biggest to smallest, compared it against all the others and combined those that were majority match?
    # clusters_all_sorted = sorted(clusters_all, key=lambda k: len(clusters_all[k]), reverse=True)

    # this may be time consuming but I dont think there should be THAT many clusters
    # for cluster1 in clusters_all_sorted:
    #     for cluster2 in clusters_all_sorted:
    #         if cluster1 != cluster2:
    #             if len(set(cluster1).intersection(set(cluster2))) > 0.5 * float(len(methods)): # if the number of matched cluster assignments is > 50%, merge them
                     

    # get all comparisons
    # comparisons = list(combinations(clusters_all.keys(),2)) # is pairwise comparisosn best here?

    # for c in comparisons:
    #     if len(set(c[0]).intersection(set(c[1]))) > 0.5 * len(c):



    


    # get methods
    # methods = list(cluster_table.columns)

    # consensus_dict = {}

    # for m in methods:

    #     # get dict of cells in each cluster
    #     clusters = cluster_table.groupby(m).groups

    #     # initial make dict of same-cluster comparisons    
    #     if len(consensus_dict.keys() == 0):
    #         for cluster,cells in clusters.items():
    #             consensus_dict.update(dict.fromkeys(list(combinations(cells, 2))),1)
    #     else:
    #         for cluster,cells in clusters.items():
    #             combos = list(combinations(cells, 2))
    #             for c in combos:
    #                 if c in consensus_dict.keys():
    #                     consensus_dict[c] += 1
    #                 else:
    #                     consensus_dict[c] = 1 # we could just not do this lolol that way we only keep track of the pairs that have agreement between all the clusters

    # return consensus_dict


'''
Make a conensus matrix from the cluster matrix. Not very efficient and takes up a lot of space.
'''
def makeConsensusMatrix(cluster_table):
    # Make the square matrix for N x N
    consensus_matrix = np.zeros(shape=(len(cluster_table), len(cluster_table)))
    consensus_matrix[:] = np.NaN  # Replace with NaN to keep the diagonal NaN 
    iteration_table = cluster_table.to_numpy()

    # Find all i,j combinations of patients that need a consensus index value
    comb = list(combinations(list(range(0, iteration_table.shape[0])), 2))

    for c in tqdm(comb):
        both_clustered = 0
        same_cluster = 0

        for i, j in zip(iteration_table[c[0]], iteration_table[c[1]]):
            if i >= 0 and j >= 0:
                both_clustered += 1

                if i == j:
                    same_cluster += 1

        res = same_cluster/both_clustered if both_clustered != 0 else 0

        consensus_matrix[c[0]][c[1]] = res
        consensus_matrix[c[1]][c[0]] = res

    return consensus_matrix


'''
Assigns each cell a cluster based on the consensus matrix/dict
'''
def assignClusters(consensus_dict):

    # for every pair
    # oof this is gonna be inefficient haha


    # for each row (cell):
    #   get all column headers (cells) that have a 1
    #   all those cells are in a cluster with this cell

    pass


'''
Create a matrix that shows if cells are in the same cluster.
Rows and columns are cell ID and the values are '1' if they are in the same cluster or '0' if they are not.

NOTE: should I use a pandas df or a scipy sparse matrix? numpy array? dictionaries?
'''
def makeMatrix(file):

    # get indices of CELLID and CLUSTER columns
    # cols = getHeader(file)
    # i_cell = cols.index(CELLID) # index of cell ID column
    # i_cluster = cols.index(CLUSTER) # indec of cluster assignment column

    # # read input
    # input = np.loadtxt(file, delimiter=',', skiprows=1, usecols=(i_cell,i_cluster))

    # read cell cluster assignments into dict
    # cells = {}
    # with open(file, mode='r') as f:
    #     reader = csv.reader(f)
    #     cells = {rows[i_cell]:rows[i_cluster] for rows in reader}
    
    # pairwise comparison...
    # for cell1,cluster1 in cells.items():
    #     for cell2,cluster2 in cells.items():
    #         if 

    '''
    hang on a faster way to do this might be to load it into a df, then group by cluster, then get the cell IDs for each group 
    and then we can already know the coordinates in our soon-to-be sparse matrix that will have a value of 1
    '''

    # load in data
    data = pd.read_csv(file, delimiter=',', index_col=CELLID)

    return data

    # # get dict of cells in each cluster
    # clusters = input.groupby(CLUSTER).groups

    # # make dict of same-cluster comparisons    
    # pairs = {}
    # for cluster,cells in clusters.items():
    #     pairs.update(dict.fromkeys((list(combinations(cells, 2))),1))


'''
Main.
'''
if __name__ == '__main__':

    # constants
    CELLID = 'CellID' # the header of the cell ID column
    CLUSTER = 'Cluster' # the header of the cluster assignmnet column
    METHOD = 'Method' # the header of the method column

    # parse arguments
    args = parseArgs()

    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output

    # if input file has correct format, make a matrix from it
    cluster_table = pd.DataFrame() # a df where cell ID is the index and each column is the cluster assignment from a different algorithm
    for file in args.input:
        if validInput(file): # check that input file has the needed columns
            data, method = getData(file)
            cluster_table[method] = data[CLUSTER]
        else:
            print(f'{file} is incorrectly formatted.')


    consensus_clusters = getConsensusClusters(cluster_table)
    ccc_assignments = formatConsensusClusters(consensus_clusters)

    # write to csv
    with open('cell_consensus_cluster_assignments.csv','w') as f:
        csv_out = csv.writer(f)
        csv_out.writerow([CELLID,CLUSTER])
        for row in ccc_assignments:
            csv_out.writerow(row)

    # table.to_csv('test_cluster.csv', sep=',')

    # make consensus matrix
    # consensus_matrix = sparse.csr_matrix(makeConsensusMatrix(cluster_table))
    # np.savetxt('test_consensus_matrix.csv', consensus_matrix, delimiter=',')

    # assign cells to clusters