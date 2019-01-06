import numpy as np
import math
from string import maketrans
from scipy.sparse import coo_matrix


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
	seq = seq.split()[0]
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
			if mklen == 1:
				seq_mat[col][int(pos)] = 1  # mark position on sequence matrix
			# get current digit representing current base for representing current sequence
		elif mklen == 1:  # erase last digit in number representing current sequence
			pos = (pos % (math.pow(alphabet_len, mk_ord)))
			pos = (pos * alphabet_len)
			pos = (pos + int(seq[cur]))
			seq_mat[col][int(pos)] = 1  # mark position on sequence matrix
		col = (col + 1)  # go to the next column in sequence matrix

	# returns matrix without spare information
	return seq_mat[mk_ord:]


def get_mk_seqMat(seq, mk_ord, uniform, offset,alphabet):  # create markov position  matrix from sequence string.
	#  used for creating pwm
	col = 0
	pos = 0
	tmp = mk_ord + 1
	seq = seq.split()[0]
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


def get_sparse_SeqMat(seq, mk_ord,alphabet):
	seq = seq.split()[0]

	alphabet = ''.join(alphabet)
	nums = range(0, len(alphabet))  # locations
	loc = ''.join(str(x) for x in nums)
	trantab = maketrans(alphabet, loc)
	nuc_list = seq.translate(trantab)

	nuc_list = np.array(map(int, nuc_list))

	kmers = np.zeros(len(nuc_list)-mk_ord)

	for digit in range(mk_ord+1)[::-1]:
		power = np.power(len(alphabet), digit)
		window_start = mk_ord - digit
		window_end = window_start + len(seq) - mk_ord
		kmers = kmers + power*nuc_list[window_start:window_end]

	rows = np.arange(0, len(seq)-mk_ord, step=1)
	cols = kmers

	data = np.ones(len(seq) - mk_ord)

	sparse_seq_mat = coo_matrix((data, (rows, cols)), shape=(len(seq) - mk_ord, np.power(len(alphabet), mk_ord+1)))

	return sparse_seq_mat



