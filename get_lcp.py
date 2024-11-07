import numpy as np
import os,sys
import matplotlib.pyplot as plt
import re
from matplotlib import cm
import matplotlib.colors as cl


"""
Determine the size and location of the largest positive and negative
charge patches by calculating l_cp for across a given sequence

also estimates f_t at physiological conditions (I=0.1M, T=300K)

usage: 
    python3 get_lcp.py *seq*.dat
"""


def ref_ft(patch):
    """
    calculate f_t30, the reference f_t if the patches were 30 residues apart
    """
    b1=0.49
    b2=8.6
    guess=1-(1/(1+np.exp(a1*(patch-a2))))
    return guess


def calc_ft(patch,sep):
    """
    use the reference f_t30 and sequence separation to calculate f_t
    """
    c=-0.019
    f0=ref_ft(patch)
    ft=f0*np.exp(y*(1-f0)*(sep-30))
    return ft


def calc_lcp(patch,pos):
    """
    calculate l_cp for a given patch from the original sequence
    """
    charge=[['R','K'],['D','E']]
    check=charge[pos]
    opp=charge[pos-1]
    N=len(patch)
    count=0
    for aa in patch:
        if aa in check:
            count+=1
        elif aa in opp:
            count-=1
    fc=count/N
    return fc*N**0.76


def window_scan(seq):
    """
    perform a scan of the input sequence calculating l_cp for patches
    of size ranging from 3 to 15 residues and save the largest positive
    and negative values as well as their positions

    input: protein sequence as alphabetic string

    output: l_cp, position, and N for strongest positive and negative patch
    """
    N=len(seq)
    results=[]
    charge=[['R','K'],['D','E']]
    wins=range(3,16)
    for c in range(2):
        value=0
        pos=0
        size=0
        for window in wins:
            for n in range(N-window):
                patch=seq[n:n+window]
                lcp=calc_lcp(patch,pos=c)
                if lcp>value:
                    value=lcp
                    pos=n+window//2
                    size=window
        results.append([value,size,pos])
    return results

if __name__ == "__main__":
    seq=sys.argv[1]
    with open(seq,'r') as fid:
        for line in fid:
            s=line
    p1,p2=window_scan(s)
    lcp1,N1,pos1=p1
    lcp2,N2,pos2=p2
    sep=abs(N1-N2)
    lcp=abs(lcp1*lcp2)
    ft=calc_ft(lcp,sep)
    print(u'$l_{cp,+}$= ',lcp1,', at residue:',pos1,', of size:',N1)
    print(u'$l_{cp,-}$= ',lcp2,', at residue:',pos2,', of size:',N2)
    print(u'$f_t$= ',ft)

