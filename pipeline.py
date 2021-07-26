import argparse

'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Calculate silhouette scores for different cluster numbers and methods from exemplar-00X quantification file.')
    parser.add_argument('-i', '--input', help='Input exemplar-00X quantification file.', type=str, required=True)
    args = parser.parse_args()
    return args


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