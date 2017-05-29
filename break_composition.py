import argparse
import numpy as np
import pandas as pd
import scipy.optimize

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input file (rows = samples, columns = otus)', required=True)
parser.add_argument('-d', help='field delimiter', default='\t')
parser.add_argument('-t', help='input data type', default='norm', choices=['counts', 'norm', 'log'])
parser.add_argument('-a', help='min coefficient', default=.1, type=float)
parser.add_argument('-b', help='max coefficient', default=10, type=float)
parser.add_argument('-o', help='output file (col1 = sample, col2 = size)', default='')
args = parser.parse_args()

# read data
x = pd.read_table(args.i, sep=args.d, header=0, index_col=0)
if args.t == 'counts':
    x = 1.*x.divide(x.sum(axis=1), axis=0)
if args.t == 'log':
    x = np.exp(x)
nrowx, ncolx = np.shape(x)

# check data
if nrowx == 0 or ncolx == 0:
    quit('error: input file format (check delimiter)')
if x.min().min() < 0 or x.max().max() > 1:
    quit('error: input file format (check data type)')

# objective function
def f(a):
    # calculate size-adjusted otu abundances
    y = np.einsum('ij,i->ij', x, a)
    # standardize otu abundances
    y = ((y - y.mean(axis=0))/y.std(axis=0))
    y = y[:,~np.isnan(y).any(axis=0)]
    # calculate the total covariance
    y = y.sum(axis=1).var() - sum(y.var(axis=0))
    return y

# initial estimates
S = np.array([1. for i in range(nrowx)]) # sizes
b = np.array([[args.a, args.b] for i in range(nrowx)]) # bounds

# run optimization
soln = scipy.optimize.minimize(f, S, bounds=b, method='L-BFGS-B')
if soln.success == True:
    S = soln.x

# write results
outline = '\n'.join(['\t'.join(zi) for zi in zip(map(str, x.index), map(str, S))])
if args.o != '':
    out = open(args.o, 'w')
    out.write(outline + '\n')
    out.close()
else:
    print outline
