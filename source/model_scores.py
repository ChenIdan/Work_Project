# function for calculating the probability that the nucleosome will sit on a specific dna sequence

import numpy as np


def Pfunc(PWM, DNA, uniform):

	new_DNA= np.transpose(DNA)

	uniform_l = np.ones((uniform, DNA.shape[1]))
	uniform_r = np.ones((uniform, DNA.shape[1]))

	PWM = np.vstack([np.vstack([uniform_l , PWM]), uniform_r]) # add ones filled rows on the upper and lower part of the matrix
	product = np.dot(PWM, new_DNA)
	return product


def get_prob(array):  # calculates probability product
	product = 1
	arr_len = len(array)
	i = 0
	while i < arr_len:  # calculate product of all array members
		if array[i] > 1:
			i = i+1
			continue
		product = product*array[i]
		i = (i+1)

	return product


def get_col(mat, pos):
	return [row[pos] for row in mat]


def creatPlist(PPM, mk_ord):  # returns probability list. get the probability matrix from Pfunc as input

	plist = []  # list of probabilities
	ppm_width = PPM.shape[1]  # get probability matrix width
	ppm_len = PPM.shape[0]  # get probability matrix length
	plist_len = abs(ppm_width-ppm_len+1) +1 # get plist len, by counting number of diagonals in ppm
	new_width = ppm_width+1
	new_ppm = np.resize(PPM, (ppm_len, new_width))  # get new probability matrix, with with all diagonals as rows
	i = 0

	while i < plist_len:  # create probability list

		cur_col = np.array(get_col(new_ppm, i))  # get current column
		cur_col = cur_col[0:ppm_len]  # take a required number of steps on the vector
		cur_prob = get_prob(cur_col)  # get probability
		plist.append(cur_prob)  # add to list
		i = (i+1)

	return plist



def score_vec(PWM, DNA, mk_ord,uniform):
	matrix = Pfunc(PWM, DNA,uniform)
	return creatPlist(matrix, mk_ord)


# calculates log(pos_dep_plist/indep_pos_plist) for every sub sequence
def calc_log(pos_dep_plist, indep_pos_plist):

	log_plist = np.log(np.divide(pos_dep_plist, indep_pos_plist))

	return log_plist



