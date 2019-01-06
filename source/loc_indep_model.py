'''
function for calculating the position independent probabilities, from a vector that its i'th pace contains the number of nucleusome covering a this location on the chromosome
'''



import seq_funcs
import numpy as np




'''
function that calculates the number of the nucleusome covering every k- mere on the nuclesome
the k in the k-mere is the integer mk_order
confg_vec is the configuration vector
the summing is done with matrix manipulation
'''
def get_sums(confg_vec, mk_order):
	mat_sums = np.tile(confg_vec,(mk_order+1,1))

	mat_sums = np.resize(mat_sums , (mat_sums.shape[0], mat_sums.shape[1]+1))

	mat_sums = np.transpose(mat_sums)

	confg_sums = np.sum(mat_sums, axis = 1)[0: len(confg_vec) - mk_order -1 ]

	return confg_sums


'''
function that calculates the independent location probabilities on of every possible k-mere
the const variable is a constant added to the sum of configuration, to prevent division by zero

'''
def old_get_model(confg_vec, chromosome, mk_order, alphabet,const, protein_len):

	#get the number of nuck
	chromosome_len = len(chromosome)


	loc_sums= get_sums(confg_vec,mk_order)

	#get seq matrix (from the proterin sequence)

	seq_mat =seq_funcs.get_seqMat(chromosome, mk_order, alphabet)

	#get the number of nucleusomes on every k-mere by multiplying the sequence matrix and the configuration vector
	confg_sums = np.dot(loc_sums, seq_mat)

	#add the constant to every cell on the configuration sums vector, to prevent division by zero

	base_vec = np.ones(len(confg_sums)) * const

	confg_sums = confg_sums+ base_vec

	#calculate Joint probabilities by normalizing the vector

	probs = confg_sums / np.sum(confg_sums)

	#calculates the dependent probabilities by

	lower_probs = np.sum(np.reshape(probs, (len(probs)/len(alphabet) ,len(alphabet) )), axis=1)

	lower_probs = np.repeat(lower_probs, len(alphabet), axis= 0)

	conditional_probs = np.divide(probs,lower_probs)

	#creat the model matrix

	conditional_mat = np.tile(conditional_probs, (protein_len -1 - mk_order, 1) )

	model_mat = np.vstack((probs, conditional_mat))

	return model_mat


def get_model(confg_vecs_file,chromosome_file, mk_order, alphabet,const, protein_len):
	confg_vecs_sum = np.zeros(np.power(len(alphabet),mk_order+1))
	confg_vecs_file = open(confg_vecs_file)
	chromosome_file=open(chromosome_file)
	for vec_line,chromosome in zip(confg_vecs_file,chromosome_file):
		confg_vec = np.fromstring(vec_line, sep=";")

		chromosome_len = len(chromosome)

		loc_sums = get_sums(confg_vec, mk_order)
		# get seq matrix (from the proterin sequence)

		seq_mat = seq_funcs.get_sparse_SeqMat(chromosome, mk_order, alphabet)

		# get the number of nucleusomes on every k-mere by multiplying the sequence matrix and the configuration vector
		confg_sums = loc_sums*seq_mat

		# add the constant to every cell on the configuration sums vector, to prevent division by zero


		confg_vecs_sum =confg_vecs_sum + confg_sums

	base_vec = np.ones(len(confg_sums)) * const

	confg_vecs_sum = confg_vecs_sum + base_vec

	#calculate Joint probabilities by normalizing the vector

	probs = confg_vecs_sum/np.sum(confg_vecs_sum )

	#calculates the dependent probabilities by

	lower_probs = np.sum(np.reshape(probs, (len(probs)/len(alphabet) ,len(alphabet) )), axis=1)

	lower_probs = np.repeat(lower_probs, len(alphabet), axis= 0)

	conditional_probs = np.divide(probs,lower_probs)

	#creat the model matrix

	conditional_mat = np.tile(conditional_probs, (protein_len -1 - mk_order, 1) )

	model_mat = np.vstack((probs, conditional_mat))

	return model_mat