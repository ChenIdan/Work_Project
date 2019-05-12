import numpy as np
from string import maketrans
from scipy.sparse import coo_matrix
from sklearn.preprocessing import LabelEncoder
import scipy.sparse as sparse
import itertools


def rc_seq(seq, alphabet, rc_alphabet):
	DNA = ''.join(alphabet)  # dna proteins
	RVERSE_DNA = ''.join(rc_alphabet)  # locations
	trantab = maketrans(DNA, RVERSE_DNA)
	rev_seq = seq.translate(trantab)  # convert the sequence into a string of numbers, according to the dictionary of tranbtab
	rev_seq = rev_seq[::-1]

	return rev_seq


def get_seqMat(seq, mk_ord,uniform, offset, alphabet):  # get regular sequence matrix.

	seq = list(seq)

	seq = seq[(uniform + offset): (len(seq) - uniform + offset)]
	diags_mat = np.zeros((len(seq) , len(seq) - mk_ord))
	i=0
	alphabet_len = len(alphabet)

	while (i <= mk_ord):
		np.fill_diagonal(diags_mat[i:], np.power(alphabet_len, mk_ord - i))
		i= i+1

	label_encoder = LabelEncoder()
	integer_encoded = label_encoder.fit_transform(seq)

	kmers = np.dot(integer_encoded, diags_mat)

	rows = np.arange(0, len(seq) - mk_ord, step=1)
	cols = kmers

	data = np.ones(len(seq) - mk_ord)

	sparse_seq_mat = coo_matrix((data, (rows, cols)), shape=(len(seq) - mk_ord, np.power(len(alphabet), mk_ord + 1)))

	return sparse.lil_matrix(sparse_seq_mat).toarray()


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

def rc_mat( mk_ord, alphabet, rc_alphabet):
	c_kmers =list(itertools.product(rc_alphabet, repeat=mk_ord+1))

	rc_kmers =  [''.join(c)[::-1] for c in c_kmers]

	DNA = ''.join(alphabet)  # dna proteins
	NUMS_DNA = ''.join((np.arange(0,len(alphabet)).astype('S10')))  # locations
	trantab = maketrans(DNA, NUMS_DNA)

	num_kmers = [int(c.translate(trantab), len(alphabet)) for c in rc_kmers]

	rows = np.arange(0, len(num_kmers), step=1)
	cols = num_kmers

	data = np.ones(len(num_kmers))

	sparse_seq_mat = coo_matrix((data, (rows, cols)), shape=(len(rows), len(rows)))

	inv_mat =  sparse.lil_matrix(sparse_seq_mat).toarray()

	return inv_mat