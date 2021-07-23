import argparse
import pandas as pd
from itertools import combinations
import csv
from operator import itemgetter
import numpy as np


'''
Cluster class
'''
class Cluster:

    # constructor
    def __init__(self, assignments, cells):
        self.assignments = assignments # a tuple where each element is the cluster assignment of a different method 
        self.cells = cells # a list of cell ID's in this cluster
        self.consensus_scores = {} # a dictionary where keys are other clusters and values are the consensus score between them

    # compare this cluster to another to calculate a consensus score
    def addComparison(self, cluster):
        num_methods = len(self.assignments) # the number of clustering methods/algorithms used
        consensus_score = 0 # a score of what fraction of all methods agreed on cluster assignments

        # compare all cluster assignments and keep score
        for i in range(num_methods):
            if self.assignments[i] == cluster.assignments[i]:
                consensus_score += 1

        # only add consensus score if it isn't zero
        if consensus_score != 0:
            consensus_score = consensus_score / num_methods # calculate consensus score NOTE: does this math work correctly or not because I'm dividing ints?
            self.consensus_scores[cluster.assignments] = consensus_score # add comparison to consensus score sdict
            


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
Write cell consensus cluster assignments to csv
'''
def writeCSV(ccc_assignments):
    with open('cell_consensus_cluster_assignments.csv','w') as f:
        csv_out = csv.writer(f)
        csv_out.writerow([CELLID,CLUSTER])
        for row in ccc_assignments:
            csv_out.writerow(row)

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

    # make cluster objects
    clusters = []
    for cluster,cells in consensus_clusters.items():
        clusters.append(Cluster(cluster, cells))

    # compare all clusters against each other
    all_comparisons = list(combinations(clusters, 2))
    for comparison in all_comparisons:
        comparison[0].addComparison(comparison[1])
        comparison[1].addComparison(comparison[0])


    # for cluster in clusters:
    #     print(f'{cluster.assignments}\tNum Cells: {len(cluster.cells)}\t Num Comparisons: {len(cluster.consensus_scores.keys())}')

    # print(f'Num Clusters: {len(clusters)}')

    import matplotlib.pyplot as plt
    num_cells = [len(c.cells) for c in clusters]
    num_comps = [len(c.consensus_scores) for c in clusters]
    plt.scatter(num_cells,num_comps)
    plt.show()

    # print(f'{clusters[2].assignments}')
    # for c, score in clusters[2].consensus_scores.items():
    #     print(f'{c}:\t{score}')

    # ccc_assignments = formatConsensusClusters(consensus_clusters)

    # write to csv
    # writeCSV(ccc_assignments)