#!/usr/local/bin/python

# functions for creating a sequence matrix from a file according to markov order

import sys
import numpy as np
import pickle
import argparse
import math
from string import maketrans
import load_PWM


def get_pwm(seq_file, mk_ord):  # get pwm matrices

    smat_list = []  # list of sequence matrices
    pwm_width = 0
    pwm_len = 0
    fh = open(seq_file, 'r')  # open sequence file
    for line in fh:  # loops though file lines
        if line[0] != '>' and line[0] != '\t' and line[0] != ' ':
            pwm_width = (pwm_width+1)  # count sequences lines, for pwm width
            if pwm_len == 0:
                pwm_len = len(line)-2  # get pwm length from sequence length

            cur_mat = get_seqMat(line, mk_ord)  # translate sequence lines to sequence matrix
            smat_list.append(cur_mat)  # insert matrices to list
    return pwm_from_seq(pwm_width, pwm_len, smat_list)

    


          
def get_seqMat(seq, mk_ord):  # create position  matrix from sequence string
    col = 0
    pos = 0
    seq_len = len(seq)-1
    mklen=int(math.pow(4, mk_ord+1))  # current markov sequence permutations number
    seq_mat = np.zeros((seq_len, mklen), dtype=np.int)  # creats zero filled matrix according to markov order#
    DNA = "ACGT"  # dna proteins
    loc = "0123"  # locations
    trantab = maketrans(DNA, loc)
    seq = seq.translate(trantab)
    for cur in range(seq_len):  # go through sequence
        if mklen != 1:
             mklen = (mklen/4)  # decreas markov sequence length
             pos = (pos+int(seq[cur])*mklen)
             # get current digit representing current base for representing current sequence
             for i in range(mklen):
                 seq_mat[col][pos+i] = 1  # mark position on sequence matrix
        elif mklen == 1:  # erase last digit in number representing current sequence
             pos = (pos%(math.pow(4, mk_ord)))
             pos = (pos*4)
             pos = (pos+int(seq[cur]))
             seq_mat[col][int(pos)] = 1  # mark position on sequence matrix
        col = (col+1)  # go to the next coloumn in sequence matrix
    return seq_mat

