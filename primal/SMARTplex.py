from . import settings
from Bio.Seq import Seq
from primer3 import calcTm

def LCSubStr(X, Y):
    # LCSuff is the table with zero value initially in each cell
    LCSuff = [[0 for k in range(len(Y)+1)] for l in range(len(X)+1)]
    # To store the length of longest common substring
    result = 0
    # Following steps to build LCSuff[m+1][n+1] in bottom up fashion
    for i in range(len(X) + 1):
        for j in range(len(Y) + 1):
            if (i == 0 or j == 0):
                LCSuff[i][j] = 0
            elif (X[i-1] == Y[j-1]):
                LCSuff[i][j] = LCSuff[i-1][j-1] + 1
                result = max(result, LCSuff[i][j])
            else:
                LCSuff[i][j] = 0
    return result

def SMARTplex(right):
    seq = right.seq
    ref = right.alignments[0].aln_ref
    for i in range(5, len(seq)):
        RTprimer = settings.RLBseq + seq[-i:]
        lcs = LCSubStr(ref, RTprimer)
        thermo = calcTm(RTprimer[-lcs:], mv_conc=75, dv_conc=3, dntp_conc=0.5)
        if thermo > 40.0:
            break
    subseq = RTprimer[-lcs:]
    lensubseq = i
    lenmatch = lcs
    return RTprimer, thermo, subseq, lensubseq, lenmatch
