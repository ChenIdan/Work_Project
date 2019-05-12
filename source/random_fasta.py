import numpy as np
from scipy import stats
from string import maketrans

# code for creating random fasta file and random sequence files

#function for creating random fasta file. seq_num is the number of sequences in the the file.

#  seq_len is the length of the sequences


def fasta_file(seq_num,seq_len, file_name,alphabet):
	seq_list = []
	xk = np.arange(0,len(alphabet))
	pk = np.array([1.0/len(alphabet)]*len(alphabet))
	custm = stats.rv_discrete(name='custm', values=(xk, pk))
	nums = range(0, len(alphabet))  # locations
	nums = ''.join(str(x) for x in nums)
	gens = alphabet

	trantab = maketrans(nums, gens)
	for num in xrange(seq_num):
		seq = "".join(map(str, custm.rvs(size=seq_len)))
		seq = seq.translate(trantab)
		seq_list.append(seq)



	np.array(seq_list).tofile(file=file_name, sep="\n>\n", format="%s")

	return seq_list



#  seq_len is the length of the sequences
def random_chrom(seq_len, file_name,alphabet, nc_num):
	seq_list = []
	xk = np.arange(0, len(alphabet))
	pk = np.array([1.0 / len(alphabet)] * len(alphabet))
	custm = stats.rv_discrete(name='custm', values=(xk, pk))
	nums = range(0, len(alphabet))  # locations
	nums = ''.join(str(x) for x in nums)
	gens = alphabet

	trantab = maketrans(nums, gens)
	seq = "".join(map(str, custm.rvs(size=seq_len)))
	seq = seq.translate(trantab)
	seq_list.append(seq)
	seq = "".join(map(str, custm.rvs(size=seq_len)))
	seq = seq.translate(trantab)
	np.array(seq).tofile(file=file_name, sep="\n>\n", format="%s")

	np.save("config_vec",np.random.randint(0,nc_num,size=seq_len-1))

	return seq




print random_chrom(10, "chromosome","ACGT",1)

fasta_file(3, 5, "/home/chenidan/nucleusome/train_files/random_fasta","ACGT")






