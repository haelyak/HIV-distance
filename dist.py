from hivSeqs import * 
from distHelper import *
import math

def jukes(propDiff):
    """Correction for multiple hits using Jukes-Cantor model.""" 
    if propDiff < 0:
        print("jukes was passed a negative number, which it doesn't know how to handle.")
        return
    elif propDiff == 0:
        return 0
    elif propDiff < 0.75:
        return -3.0*math.log(1-(4.0*propDiff/3))/4
    else:
        print("jukes was passed a number >= 0.75, which is too big.")
        return


def propDifferent(seqA,seqB):
    """propDifferent
    Inputs: two different sequences
    Outputs: proportion of differences 
    between a pair of aligned sequences considering only non-gap sites."""
    i = 0
    same = 0
    prop = 0
    gap = 0
    while i < len(seqA):
        if seqA[i] == '-' or seqB[i] == '-':
            gap += 1 
        elif seqA[i] == seqB[i]:
            same += 1
        i+=1
    prop = (len(seqA) - same - gap)/(len(seqA) - gap)
    return prop


def distances(strainNamesL,seqsL):
    """distances
    Inputs: strainNamesL, list of strain names to act as keys in the dictionary
            seqsL, list of sequences
    Outputs: a dictionary with keys of tuples of two strain names and the distances between pairs of sequences
    """
    D = {(strainNamesL[0], strainNamesL[0]): 0}
    k = 0
    
    while k < len(strainNamesL):
        l = 0
        while l < len(strainNamesL):
            D[(strainNamesL[k], strainNamesL[l])] = jukes(propDifferent (seqsL[k], seqsL[l]))
            l += 1
        k += 1
    return D
