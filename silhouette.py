from sklearn.metrics import silhouette_score
import pandas as pd
import argparse


'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='calculate silhouette scores of various clustering methods')
    parser.add_argument('-f', '--features', help='feature table', type=str, required=True)
    parser.add_argument('-c', '--clusters', help='cluster labels from different methods', nargs='*', action='store', dest='clusters', required=False)
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
Main.
'''
if __name__ == '__main__':

    # constants
    CELLID = 'CellID' # the header of the cell ID column
    CLUSTER = 'Cluster' # the header of the cluster assignmnet column
    METHOD = 'Method' # the header of the method column

    # parse arguments
    args = parseArgs()

    # get feature table that is already appropriately transformed, if necessary
    feature_table = pd.read_csv(args.features, index_col=CELLID)

    # get cluster label table
    cluster_table = pd.DataFrame() # a df where cell ID is the index and each column is the cluster assignment from a different algorithm
    for file in args.clusters:
        if validInput(file): # check that input file has the needed columns
            data, method = getData(file)
            cluster_table[method] = data[CLUSTER]
        else:
            print(f'{file} is incorrectly formatted.')

    # calculate silhouette scores
    # what distance metric do I use???
    for method in cluster_table.columns:
        s = silhouette_score(feature_table, cluster_table[method])
        print(f'{method}:\t{s}')