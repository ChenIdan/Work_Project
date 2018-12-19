import numpy as np
import math
from string import maketrans



def rc_seq(seq, alphabet, rc_alphabet):
	DNA = ''.join(alphabet)  # dna proteins
	RVERSE_DNA = ''.join(rc_alphabet)  # locations
	trantab = maketrans(DNA, RVERSE_DNA)
	rev_seq = seq.translate(trantab)  # convert the sequence into a string of numbers, according to the dictionary of tranbtab
	rev_seq = rev_seq[::-1]

	return rev_seq

def get_seqMat(seq, mk_ord,alphabet):  # get regular sequence matrix.
	col = 0
	pos = 0
	seq_len = len(seq)
	alphabet_len = len(alphabet)
	mklen = int(math.pow(alphabet_len, mk_ord + 1))  # current markov sequence permutations number
	seq_mat = np.zeros((seq_len, mklen), dtype=np.int)  # creates zero filled matrix according to markov order#
	DNA = ''.join(alphabet)  # dna proteins
	nums = range(0, alphabet_len)  # locations
	loc = ''.join(str(x) for x in nums)
	trantab = maketrans(DNA, loc)
	seq = seq.translate(trantab)
	for cur in range(seq_len):  # go through sequence
		if mklen != 1:
			mklen = (mklen / alphabet_len)  # decrease markov sequence length
			pos = (pos + int(seq[cur]) * mklen)
			# get current digit representing current base for representing current sequence
		elif mklen == 1:  # erase last digit in number representing current sequence
			pos = (pos % (math.pow(alphabet_len, mk_ord)))
			pos = (pos * alphabet_len)
			pos = (pos + int(seq[cur]))
			seq_mat[col][int(pos)] = 1  # mark position on sequence matrix
		col = (col + 1)  # go to the next column in sequence matrix

	#returns matrix without spare information
	return seq_mat[mk_ord+1:]


def get_mk_seqMat(seq, mk_ord, uniform, offset,alphabet):  # create markov position  matrix from sequence string.
	#  used for creating pwm
	col = 0
	pos = 0
	tmp = mk_ord + 1
	seq_len = len(seq)  - 2 * uniform
	alphabet_len = len(alphabet)
	mklen = int(math.pow(alphabet_len, tmp))  # current markov sequence permutations number
	seq_mat = np.zeros((seq_len, mklen))  # create zero filled matrix according to markov order#
	DNA = ''.join(alphabet) # dna proteins
	nums = range(0, alphabet_len)  # locations
	loc = ''.join(str(x) for x in nums)
	trantab = maketrans(DNA, loc)
	seq = seq.translate(trantab)
	for cur in range(uniform + offset, uniform + offset + seq_len):  # go through sequence
		if mklen != 1:
			mklen = (mklen / alphabet_len)  # decrease markov sequence length
			pos = (pos + int(seq[cur]) * mklen)
			# get current digit representing current base for representing current sequence
			for i in range(mklen):
				seq_mat[col][pos + i] = 1  # mark position on sequence matrix
		elif mklen == 1:  # erase last digit in number representing current sequence
			pos = (pos % (math.pow(alphabet_len, mk_ord)))
			pos = (pos * alphabet_len)
			pos = (pos + int(seq[cur]))
			seq_mat[col][int(pos)] = 1  # mark position on sequence matrix
		col = (col + 1)  # go to the next column in sequence matrix
	return seq_mat
