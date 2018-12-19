'''
this module contains test function
'''


import pandas as pd
import system_funcs
import os
import matplotlib.pyplot as plt
import numpy as np
import pwm
import plotly.plotly as py
import plotly.graph_objs as go
import plotly
from string import maketrans


sep = os.sep
base_dir = system_funcs.prog_root_dir() +sep
segal_dir = base_dir +'segal_results'+sep
kaplan_dir = base_dir+'kaplan_results'+sep
segal_mats_dir =  base_dir +'segal_model_matrices'+sep
seq_dir = base_dir+'sequence_file'+sep
train_dir = base_dir + 'train_files' + sep


argv_input_gxw_file_check = 'python model.py -pos_indep_model_generate_method input_file  -pos_indep_model_file_method numpy -pos_indep_model_prob_output_method  numpy ' \
	   '-pos_indep_model_input_file ' +segal_mats_dir+ 'segal_pos_indep_pwm.npy ' \
'-pos_indep_model_markov_order 4 -pos_indep_model_prob_output_file ' + segal_dir +'pos_indep_probs ' \
'-pos_dep_model_file_method numpy' \
	   ' -pos_dep_model_input_file ' + segal_mats_dir+ 'segal_pos_dep_pwm.npy' \
	   ' -uniform 0 -offset 0 -seq_file ' + seq_dir+'seq_input ' \
'-pos_dep_model_generate_method input_file -pos_dep_model_prob_output_method numpy -pos_dep_model_markov_order 1' \
' -pos_dep_model_prob_output_file ' + segal_dir+ 'pos_dep_probs -temp_param 1 -nc_concentration  1 -row_binding_output_method numpy' \
' -row_binding_file_name ' \
 +segal_dir+'kaplan_scores_from_segal_pwm '



argv_for_train_file_check = 'python model.py -pos_indep_model_generate_method train_file -pos_indep_model_prob_output_method  numpy ' \
'-pos_indep_model_markov_order 4 -pos_indep_model_prob_output_file pos_indep_probs '  \
	   ' -uniform 0 -offset 0 -seq_file ' + seq_dir+ 'seq_input ' \
'-pos_dep_model_generate_method train_file -pos_dep_model_prob_output_method numpy -pos_dep_model_markov_order 1' \
' -pos_dep_model_prob_output_file pos_dep_probs -temp_param 1 -nc_concentration  1 -row_binding_output_method numpy' \
' -row_binding_file_name ' + kaplan_dir+ \
'kaplan_scores  -train_file ' +train_dir+ 'fasta_file1'


def get_segal_scores(scores_file):
	base_path = system_funcs.prog_root_dir()
	scores_file_path = base_path + os.sep + scores_file
	table = pd.read_csv(scores_file_path, sep='\t', lineterminator='\n')

	return table['Raw Binding (log ratio)']


def get_segal_probs(probs_file):
	base_path = system_funcs.prog_root_dir()
	probs_file_path = base_path + os.sep + probs_file
	table = pd.read_csv(probs_file_path, sep='\t', lineterminator='\n')
	return {'nucleosome_start_probs' : table['P start'],
	        'nucleosome_sit_probs': table['P occupied']}

def check_model_symmetry():
	os.system("python load_gxw.py")
	base_path = system_funcs.prog_root_dir()
	seq_file_path = base_path+ os.sep+"sequence_file" + os.sep + "seq_input"
	fh = open(seq_file_path,'r+')
	buf_size = 0
	for line in fh:
		if line[0] == 'A' or line[0] == 'C' or line[0] == 'G' or line[0] == 'T':
			seq = line.splitlines()[0]
			break
		else:
			buf_size = buf_size +len(line)
	os.system(argv_input_gxw_file_check)
	seq_graph = go.Scatter(x=range(len(np.load("/home/chen/Desktop/Work_Project/segal_results/pos_indep_probs.npy"))), y=np.load("/home/chen/Desktop/Work_Project/segal_results/pos_indep_probs.npy"),
	                           name='seq_probs')

	fh.seek(buf_size)
	fh.write(pwm.rc_seq(seq))
	fh.close()


	#seq_graph = go.Scatter(x=range(len(np.load("/home/chen/Desktop/Work_Project/segal_results/pos_dep_probs.npy"))), y=np.load("/home/chen/Desktop/Work_Project/segal_results/pos_dep_probs.npy"),
	 #                      name='seq_probs')



	os.system(argv_input_gxw_file_check)




	#rev_seq_graph = go.Scatter(x=range(len(np.load("/home/chen/Desktop/Work_Project/segal_results/pos_dep_probs.npy")))
	 #                          , y=np.load("/home/chen/Desktop/Work_Project/segal_results/pos_dep_probs.npy"),
	  #                     name='rev_seq_probs')

	rev_seq_graph = go.Scatter(x= range(len(np.load("/home/chen/Desktop/Work_Project/segal_results/pos_indep_probs.npy"))), y = np.load("/home/chen/Desktop/Work_Project/segal_results/pos_indep_probs.npy")[::-1],name='rev_seq_probs')
	# diff_vec = go.Scatter(x=range(start,end),
	#               y=diff_vec[start:end],
	#              name='differences between kaplan and segal model' )

	layout = go.Layout(title="symmetry_check_kaplan",
	                   xaxis=dict(title='sequence bases'),
	                   yaxis=dict(title='Nucleusome_probs'))

	data = [seq_graph, rev_seq_graph]

	fig = go.Figure(data=data, layout=layout)



	py.plot(fig, sharing='public', filename='kaplan_symmetry_check')



def check_segal_symmetry():
	cur_path = os.getcwd()
	os.chdir("/home/chen/nucleosome_model/")
	os.system("./change_parameters.sh 1 1")
	os.chdir(cur_path)
	fh = open("/home/chen/nucleosome_model/input.fa", 'r+')
	buf_size =0
	pos_dep_probs = np.load("/home/chen/Desktop/Work_Project/segal_results/pos_indep_probs.npy")
	for line in fh:
		if line[0] == 'A' or line[0] == 'C' or line[0] == 'G' or line[0] == 'T':
			seq = line.splitlines()[0]
			break
		else:
			buf_size = buf_size+len(line)
	segal_scores = get_segal_scores("scores.tab")
	fh.seek(0)
	fh.read(buf_size)

	fh.write(pwm.rc_seq(seq))
	fh.close()

	seq_graph = go.Scatter(x=range(len(segal_scores)), y=segal_scores ,
	                       name='segal_scores')

	cur_path = os.getcwd()
	os.chdir("/home/chen/nucleosome_model/")
	os.system("./change_parameters.sh 1 1")
	os.chdir(cur_path)
	rev_segal_scores = get_segal_scores("scores.tab")

	rev_seq_graph = go.Scatter(x=range(len(segal_scores)), y=rev_segal_scores[::-1],
	                       name='rev_segal_scores')




	#for line in fh:
		#if line[0] == 'A' or line[0] == 'C' or line[0] == 'G' or line[0] == 'T':
			#fh.write( pwm.rc_seq(seq))


	layout = go.Layout(title="symmetry_check_segal",
	                   xaxis=dict(title='sequence bases'),
	                   yaxis=dict(title='Nucleusome_probs'))

	data = [seq_graph, rev_seq_graph]

	fig = go.Figure(data=data, layout=layout)


	py.plot(fig, sharing='public', filename='segal_symmetry_check')


def kaplan_segal_score_cmp():
	os.system("python load_gxw.py")
	os.system(argv_input_gxw_file_check)
	cur_path = os.getcwd()
	os.chdir("/home/chen/nucleosome_model/")
	os.system("./change_parameters.sh 1 1")
	os.chdir(cur_path)
	segal_scores = get_segal_scores("scores.tab")

	kaplan_scores = np.load("/home/chen/Desktop/Work_Project/segal_results/kaplan_scores_from_segal_pwm.npy")

	print segal_scores - kaplan_scores
	kaplan_graph = go.Scatter(x=range(len(kaplan_scores)), y=kaplan_scores,
	                          name='kaplan_scores')

	segal_graph = go.Scatter(x=range(len(segal_scores)), y=segal_scores,
	                         name='segal_scores')

	layout = go.Layout(title="kaplan vs segal",
	                   xaxis=dict(title='sequence bases'),
	                   yaxis=dict(title='Nucleusome_probs'))

	data = [kaplan_graph, segal_graph]

	fig = go.Figure(data=data, layout=layout)

	py.plot(fig, sharing='public', filename='kaplan and segal scores')

def test_matmul_naive(pwm, seq):
	mk_order = int(np.log2(len(pwm[0]))/2 - 1)

	array_len = len(seq) - len(pwm) -mk_order + 1

	DNA = 'ACGT'  # dna proteins
	base = '0123'  # locations
	trantab = maketrans(DNA,base)
	seq_nums = seq.translate(
		trantab)  # convert the sequence into a string of numbers, according to the dictionary of tranbtab
	probs_vec = np.ones(array_len, dtype=np.float)

	for i in range(array_len):
		for nc in range(i,i+len(pwm)):
			cell = int(seq_nums[nc:nc+mk_order+1], base=4)
			probs_vec[i] = probs_vec[i]*pwm[nc - i, cell]



	return probs_vec

def test_scores_naive(pos_indep_pwm,pos_dep_pwm, seq):
	loc_dep_probs =test_matmul_naive(pos_dep_pwm,seq)

	loc_indep_probs = test_matmul_naive(pos_indep_pwm,seq)

	kaplan_naive_scores =  np.log((loc_dep_probs/loc_indep_probs))

	kaplan_scores = np.load("/home/chen/Desktop/Work_Project/kaplan_results/kaplan_scores.npy")

	naive_graph = go.Scatter(x=range(len(kaplan_naive_scores )), y=kaplan_naive_scores ,
                          name='model_naive_scores')

	kaplan_graph = go.Scatter(x=range(len(kaplan_scores)), y=kaplan_scores,
                         name='model_scores')

	layout = go.Layout(title="kaplan vs naive",
                   xaxis=dict(title='sequence bases'),
                   yaxis=dict(title='binding scores'))

	data = [naive_graph, kaplan_graph ]

	fig = go.Figure(data=data, layout=layout)

	py.plot(fig, sharing='public', filename='kaplan and segal binding scores')




def kaplan_segal_probs_cmp():
	os.system("python load_gxw.py")
	os.system("python model.py")
	cur_path = os.getcwd()
	os.chdir("/home/chen/nucleosome_model/")
	os.system("./change_parameters.sh 1 1")
	os.chdir(cur_path)
	segal_sit_probs= get_segal_probs("occ_probs.tab")['nucleosome_sit_probs']
	segal_start_probs = get_segal_probs("occ_probs.tab")['nucleosome_start_probs']
	kaplan_sit = np.load("/home/chen/Desktop/Work_Project/source/kaplan_nucleusome_sit.npy")
	kaplan_start = np.load("/home/chen/Desktop/Work_Project/source/kaplan_nusleusome_start.npy")

	kaplan_graph = go.Scatter(x=range(len(kaplan_sit )), y=kaplan_sit ,
	                          name='kaplan_sit_probs')

	segal_graph = go.Scatter(x=range(len(segal_sit_probs )), y=segal_sit_probs,
	                         name='segal_sit_probs')

	layout = go.Layout(title="kaplan vs segal",
	                   xaxis=dict(title='sequence bases'),
	                   yaxis=dict(title='Nucleusome_probs'))

	data = [kaplan_graph, segal_graph]

	fig = go.Figure(data=data, layout=layout)

	py.plot(fig, sharing='public', filename='kaplan and segal sit probs')



test_scores_naive(np.load("/home/chen/Desktop/Work_Project/trained_models/trained_pos_indep_model.npy"),np.load("/home/chen/Desktop/Work_Project/trained_models/trained_pos_dep_model.npy"),"TTTCCGGTGTTGGTCAATCCCTGAATAGCGAACGATTATAATGAGACTGGAGTGTAATCGCTAGCGACGACCGAGAGGTTCTTAGCGGCATCGTAGTTACCCCTCGGCTCCAGTAGCTCACTAATGGTCGCCGATCTGATAGATTGCTCGGTAATCTCCGTCTGGTCTCGTCTGCATGGCGGACTTTATAGATCTAGTATCGGTCTGGTTATAACGATGCTGAGCATAAGCGTGGCTGAAAACCGGCGCATAAAGGGCAAATCCGAAAGCAACGAGTGTCCTGTGGACCATATGGATAACTAGAGTACCCCGCACGCTAGCATGACACATGCCTCCGCGAGGCCATTCTTCGATACAGTAAGAATAGATTCATTCGCCTATTTCTCTGTTTTTTTGCGCTATACTTACGTTCCCCTTAATTGCCTGCCGGAGCTAAGGACGTTCCATTGCGTGGGGCCGGGTTCTCGGCAATTCGCTATATTAGCCCGAAGACCAGTCGGGCAGGACCCGAGAAGTTGGTTGCCAGTGCCCCCGAATAATCCTCCTAAAACTATGCCTGGGCGTTGGCAGCAGTCCGTCGAAACTTCGGGCAGAGAGGCCTAGCCGGCTCAAGCGAGTGGAAAAACTTAGGAGCGTCCGGTATTGATTGTCGGGCAAATACCTTCTCATACGTGACGGATTCCCGCCTCCCGGTTCTCAGGAGGGTTAAAGGTGTAATTTCTATAGACGTGACCGCAGGGGAGAAGAGCTCCTGCTCCCATTTTTTGGCGGTAACGTGATCCTGAATGGCGCCTTTGACGGAGCAACGTACATAACGTCGGAGTCCCTGCAGTCATCTAGTACCTCCATCTCGTACACAGGAGACGACAATTTTACACTCAAGACTTGTTTCATGGTTAACCAGACTTGAATGTGATTCCTCTAAACCATCCACGTGCCAGAGACAGCATATCTGGTGAAGTAGACTCTACAGAAGGCCATAGTCAGTTAGGGACTTATC")

