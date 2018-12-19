#all the algorithms are in the paper "feild et al" at pages 20 - 22

import numpy as np

#the probability of a nucleosome to start on the location of the input variable loc

#input: #loc - the starting location of the nucleosome
		#seq - the DNA sequence
		#score_vec: vector of probabities score. score_vec[i] is the probability of a nucleosome to sit on sqe[i,....,i+146]

#output: probability of a nucleusome to start on location in variable loc in the sequence seq


def prob_loc_start(cur_loc, seq_len, score_vec, forwards, backwards, log_nc_concentarion, temp_param):

	protein_length = seq_len - len(score_vec) + 1
	if cur_loc==0:
		return log_nc_concentarion+(score_vec[0]*temp_param)+backwards[protein_length] - backwards[0]
		#return (nc_concentarion*np.exp(score_vec[0]*temp_param)*backwards[147])/backwards[0]
	elif seq_len <=protein_length+cur_loc or cur_loc < 0:
		return np.NINF
	else:
		return forwards[cur_loc-1]+log_nc_concentarion+score_vec[cur_loc]*temp_param+backwards[cur_loc+protein_length] - backwards[0]
		#return (forwards[loc-1]*nc_concentarion*np.exp(score_vec[loc]*temp_param) * backwards[loc+147])/backwards[0]


def prob_loc_sit(cur_loc, seq_len, score_vec, forwards, backwards, log_nc_concentarion, temp_param):
	protein_length = seq_len - len(score_vec) +1
	prob = np.NINF
	for pos in range(cur_loc - protein_length + 1 , cur_loc +1):
		prob = np.logaddexp(prob, prob_loc_start(pos,seq_len,score_vec,forwards,backwards, log_nc_concentarion, temp_param))


	return prob


#method for summing the statitistical weight of all posibbile nucleosome configurations,
#so in the end, loc_prob[i] will contain the sum off the statistical weights of all the nuclesome configurations
#at seq[1,....,i]

#input: #loc_prob - initialize with the weights sums by the function
		#seq_len - sequence length
		#score_vec - the statistical score given to the nucleosome sitting at seq[i,...,i+146]
def sum_forwards(seq_len,score_vec, log_nc_concentarion, temp_param):
	probs_sum = np.zeros(seq_len)
	protein_len = seq_len -len(score_vec)
	for i in range(protein_len,seq_len):
		if i == protein_len:
			probs_sum[i] = log_nc_concentarion+score_vec[i-146]*temp_param
		else:
			probs_sum[i]=  np.logaddexp(probs_sum[i-1],log_nc_concentarion+probs_sum[i-147]+score_vec[i-146]*temp_param)
	return probs_sum

#method for summing the statitistical weight of all posibbile nucleosome configurations,
#so in the end, loc_prob[i] will contain the sum off the statistical weights of all the nuclesome configurations
#at seq[i,....,seq_len]

#input: #loc_prob -  initialize with the weights sums by the function
		#seq_len - sequence length
		#score_vec - the statistical score given to the nucleosome sitting at seq[i,...,i+146]
def sum_backwards(seq_len, score_vec, log_nc_concentarion, temp_param):
	probs_sum = np.zeros(seq_len)
	loop_range = range(0,len(score_vec))[::-1]
	for i in loop_range:
		if i == len(score_vec) - 1:
			probs_sum[i] = log_nc_concentarion +score_vec[i]*temp_param
		else:
			probs_sum[i] = np.logaddexp(probs_sum[i+1],log_nc_concentarion +probs_sum[i+147]+score_vec[i]*temp_param)
	return probs_sum




