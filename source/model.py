#!/home/chen/miniconda2/bin/python

# main program of the model.
# input is pwm matrix in some file, and a position matrix, also given
# in a file.
# output is probability score vector

import argparse
import sys
import model_scores
import pickle
import numpy as np
import model_probs
from string import maketrans
import plotly
import pickle
import pandas as pd
import os
plotly.__version__
import system_funcs
import loc_indep_model
import seq_funcs
import loc_dep_model

import matplotlib.pyplot as plt



def main(argv):  # main program
	print argv


	sep = os.sep
	base_dir = system_funcs.prog_root_dir() + sep

	# code for getting program input
	parser = argparse.ArgumentParser(prog="model.py", description="calculate probabilities of biological molecules",
									 argument_default=None, add_help=True)  # creating a parser for terminal args


#arguments for generating the position independent model matrix

	parser.add_argument("-pos_indep_model_generate_method", type = str, dest="pos_indep_model_gen_method", choices=['input_file', 'train_file_pos_indep'],
						help= "this command will decide how to generate the position independent model matrix. you can give the input 'input file' "
						"and then the model matrix will be created from an input file you will insert by the command 'pos_indep_model_file or you can insert the input "
						"train_file' and then the model will be created from a fasta file you will insert")

	parser.add_argument("-pos_indep_model_file_method", type = str, dest = "pos_indep_model_file_method", choices = ['tabs', 'pickle', 'numpy'],
					   help = "if you decided creating the position independent matrix with an input file containing the model matrix, you need to choose the type of file - "
							  "tabs, pickle, or numpy")

	parser.add_argument("-pos_indep_model_input_file", type=str, dest="pos_indep_model_input_file",
					   help="file containing the position independent model matrix")

	parser.add_argument("-pos_indep_model_markov_order", type=int, dest="pos_indep_model_markov_order", help= "position independent model markov order", default = 1)

	parser.add_argument("-pos_indep_model_prob_output_file", type=str, dest="pos_indep_model_prob_output_file",
						help="name of the output file of that will contain the probabilities calculated from the position independent model")
	parser.add_argument("-pos_indep_model_prob_output_method", type=str, dest="pos_indep_model_prob_output_method",
						help="the method which the probabilities calculated by the position independent model will be saved" , choices=["tabs", "numpy", "pickle"])

	parser.add_argument("-pos_indep_model_prob_train_file", type=str, dest="pos_indep_model_prob_train_file",
						help="train file for generating the position independent model")

	# arguments for generating the position dependent model matrix

	parser.add_argument("-pos_dep_model_generate_method", type = str, dest="pos_dep_model_gen_method", choices=["input_file", "train_file"], help=
						"this command will decide how to generate the position independent model matrix. you can give the input 'input file' "
						"and then the model matrix will be created from an input file you will insert by the command 'pos_dep_model_file or you can insert the input "
						"train_file' and then the model will be created from a fasta file you will insert")

	parser.add_argument("-pos_dep_model_file_method", type = str, dest = "pos_dep_model_file_method", choices = ['tabs', 'pickle',                                                      'numpy'],
					   help = "if you decided creating the position dependent model matrix with an input file containing the model matrix, you need to choose the type of file - "
							  "tabs, pickle, or numpy")

	parser.add_argument("-pos_dep_model_input_file", type=str, dest="pos_dep_model_input_file",
					   help="file containing the position independent model matrix")

	parser.add_argument("-pos_dep_model_markov_order", type = int, dest="pos_dep_model_markov_order", help= "position dependent model markov order", default = 5)

	parser.add_argument("-pos_dep_model_prob_output_file", type=str, dest="pos_dep_model_prob_output_file",
						help="name of the output file of that will contain the probabilities calculated from the position dependent model")
	parser.add_argument("-pos_dep_model_prob_output_method", type=str, dest="pos_dep_model_prob_output_method",
						help="the method which the probabilities calculated by the position dependent model will be saved. you can choose numpy, tabs, or pickle" , choices=["tabs", "numpy", "pickle"])


	#arguments for calculating probabilities on sequence

	parser.add_argument("-seq_file", type=str, dest="seq_file", help = "file containing the dna sequence")
	parser.add_argument("-uniform", type = int, dest= "uniform", help= "uniform part of the input nucleusome  sequence")

	parser.add_argument("-offset", type = int, dest= "offset", help = "offset of the uniform part of the nucleusome dna sequence")

	parser.add_argument("-nc_concentration", type = float, help = "nucleusome concentration")

	parser.add_argument("-temp_param", type = float, help = "temperature parameter")

	parser.add_argument("-row_binding_file_name" , type = str, dest= "row_binding_file_name", help= "row binding score file name")
	parser.add_argument("-row_binding_output_method", type = str, dest= "row_binding_output_method",
						choices= ["tabs", "pickle", "numpy"], help= "row binding score output method. "
																"the possible method are pickle, numpy,"
																"or tabs delimited file ")

	parser.add_argument("-train_file", type=str, dest="train_file",
						help="train file for generating the position dependent model and position dependent model")

	parser.add_argument("-chromosome_seq", type=str, dest="chrome_seq",
	                    help="chromosome sequence, used for calculating position independet probabilities")
	parser.add_argument("-config_vec", type=str, dest= "config_vec", help= "each cell contains the number of configuration in which a nucleutid on the chromosome -choromosome_seq is covered by a nucleusome")

	parser.add_argument("-alphabet",type=str,dest="alphabet", help="nucleutide alphabet")

	parser.add_argument("-rc_alphabet", type=str, dest="rc_alphabet",help="the compiments of the alphabet")

	parser.add_argument("-loc_indep_const", type=float, help="constant preventing division by zero")



	dict_prob_method = {'loc_dep_prob': 0, 'loc_indep_prob': 1}  # dictionary for probability methods for model

	print sys.argv
	args = parser.parse_args(argv)  # get command line arguments
	uniform = args.uniform
	offset = args.offset
	print args
	alphabet = (args.alphabet).split(',')
	rc_alphabet = (args.rc_alphabet).split(',')

	if args.pos_indep_model_markov_order:
		loc_indep_mk_order = args.pos_indep_model_markov_order

	#code for getting position independent model probability vector

	#code for generating pwm from input file
	if args.pos_indep_model_gen_method == "input_file":
		pos_indep_model_file = args.pos_indep_model_input_file
		if args.pos_indep_model_file_method == "numpy":
			loc_indep_model_mat = np.load(pos_indep_model_file)
		elif args.pos_indep_model_file_method == "pickle":
			loc_indep_model_mat = pickle.load(pos_indep_model_file)
		elif args.pos_indep_model_file_method == "tabs":
			loc_indep_model_mat = np.loadtxt(pos_indep_model_file)
		else:
			sys.exit( "method for generating input file in unknown. pleas enter the file with "
					  "the methods: pickle, numpy, or tab dilimeted text file\n")


	#code for generating pwm from train file
	elif args.pos_indep_model_gen_method == "train_file_pos_indep":
		confg_vec = np.load(args.config_vec)

		with open(args.chrome_seq, 'r') as chrome_file:
			chromosme_seq= chrome_file.read()


		const = args.loc_indep_const

		loc_indep_model_mat =  loc_indep_model.get_model(confg_vec ,chromosme_seq,loc_indep_mk_order, alphabet, const, 147)
		np.save(base_dir+"trained_models"+ sep + "trained_pos_indep_model", loc_indep_model_mat )



	else:
		sys.exit("You haven't entered method for generating position independent model. the model is entered with the flags '-position_indep_gen_method ', and options are 'input_file' or 'train_file!\n ")

	pos_indep_model_row_len = loc_indep_model_mat.shape[1]
	pos_indep_mat_mk_order = loc_indep_mk_order

	if pos_indep_mat_mk_order != loc_indep_mk_order:
			exit("markov order of the position independent model is not equal to the markov order of the model input matrix!\n")


	#code for generating location independent probabilities from the pwm
	#loc_indep_model_mat = loc_indep_model_mat[uniform: 147 - uniform]
	if args.seq_file:
		with open(args.seq_file, 'r+') as f:
			seq_str = f.readline().splitlines()[0]
			loc_indep_seq_mat = seq_funcs.get_seqMat(seq_str, loc_indep_mk_order,alphabet)

		loc_indep_score_vec = model_scores.score_vec(loc_indep_model_mat, loc_indep_seq_mat, loc_indep_mk_order, uniform)

		f.close()
	else:
		print 1
		exit("pleas enter a dna sequence file!\n")

	if args.pos_indep_model_prob_output_file:
		pos_indep_prob_output_file_name = args.pos_indep_model_prob_output_file
		if args.pos_indep_model_prob_output_method == "pickle":

			with open((pos_indep_prob_output_file_name + ".pickle"), 'w') as fh:
				pickle.dump(loc_indep_score_vec , fh)
				fh.close()
		elif args.pos_indep_model_prob_output_method == "numpy":
			np.save(pos_indep_prob_output_file_name+".npy", loc_indep_score_vec)
		elif args.pos_indep_model_prob_output_method  == "tabs":
			np.savetxt(pos_indep_prob_output_file_name+".tabs" , loc_indep_score_vec)
		else:
			exit("pleas choose method for saving you position independent probability vector. posibile methods are tabs, pickle, or numpy")






	# code for getting position dependent model probability vector

		# code for getting position independent model probability vector

		# code for generating pwm from input file
	if args.pos_dep_model_gen_method == "input_file":
		pos_dep_model_file = args.pos_dep_model_input_file
		if args.pos_indep_model_file_method == "numpy":
			loc_dep_model_mat = np.load(pos_dep_model_file)
		elif args.pos_dep_model_file_method == "pickle":
			loc_dep_model_mat = pickle.load(pos_dep_model_file)
		elif args.pos_indep_model_file_method == "tabs":
			loc_dep_model_mat = np.loadtxt(pos_dep_model_file)
		else:
			exit("method for generating input file for position dependent model in unknown. pleas enter                           the file with "
						 "the methods: pickle, numpy, or tab dilimeted text file\n")






		# code for generating pwm from train file
	elif args.pos_dep_model_gen_method== "train_file":

			if args.train_file:
				loc_dep_model_mat = loc_dep_model.pos_dep_con_mat(args.train_file, args.pos_dep_model_markov_order, uniform, offset,alphabet,rc_alphabet)
				np.save(base_dir + "trained_models" + sep + "trained_pos_dep_model", loc_dep_model_mat)
			else:
				exit("you haven't entered any train file!\n")

	if args.pos_dep_model_markov_order:
			loc_dep_mk_order = args.pos_dep_model_markov_order

	else:
		sys.exit("You haven't entered method for generating position independent model. the model is entered "
		         " with the flags '-position_indep_gen_method ', and options are 'input_file' or "
		         " 'train_file!\n ")

	pos_dep_model_row_len = loc_dep_model_mat.shape[1]
	pos_dep_mat_mk_order = loc_dep_mk_order

	if pos_dep_mat_mk_order != loc_dep_mk_order :
		exit("markov order of the position dependent model is not equal to the markov order of the model input matrix!\n")

		# code for generating location dependent probabilities from the pwm
	loc_dep_model_mat = loc_dep_model_mat

	if args.seq_file:
		with open(args.seq_file, 'r+') as f:
			seq_str = f.readline().splitlines()[0]
			loc_dep_seq_mat = seq_funcs.get_seqMat(seq_str, loc_dep_mk_order,alphabet)

			loc_dep_score_vec = model_scores.score_vec(loc_dep_model_mat, loc_dep_seq_mat, loc_dep_mk_order, uniform)

			f.close()
	else:
		exit("pleas enter a dna sequence file!\n")

	if args.pos_dep_model_prob_output_file:
		pos_dep_prob_output_file_name = args.pos_dep_model_prob_output_file
		if args.pos_dep_model_prob_output_method == "pickle":
			with open((pos_dep_prob_output_file_name + ".pickle"), 'w') as fh:
				pickle.dump(loc_dep_score_vec, fh)
		elif args.pos_dep_model_prob_output_method == "numpy":
			np.save(pos_dep_prob_output_file_name,loc_dep_score_vec)
		elif args.pos_dep_model_prob_output_method == "tabs":
			np.savetxt(pos_dep_prob_output_file_name, loc_dep_score_vec)
		else:
			exit("pleas choose method for saving you position dependent probability vector. possible methods                        are tabs, pickle, or numpy")


	# now the score vector is calculated from the position independent model and the position dependent model

	score_vec = np.log(np.divide(loc_dep_score_vec, loc_indep_score_vec))

	if args.row_binding_file_name:
		raw_binding_file_name = args.row_binding_file_name
		if args.row_binding_output_method == "pickle":
			with open((raw_binding_file_name  + ".pickle"), 'w') as fh:
				pickle.dump(score_vec, fh)
		elif args.row_binding_output_method == "numpy":
			np.save(args.row_binding_file_name, score_vec)
		elif  args.row_binding_output_method == "tabs":
			np.savetxt(args.row_binding_file_name, score_vec, delimiter= '\t')
		else:
			exit("you entered no known output method for saving score vector ")

	#Now calculating probabilities with dynamic programming methods, calculated from the score vector, generted from the two probabilities model...
	seq_len = len(seq_str)

	start_probability = range(seq_len)

	sit_probability = range(seq_len)

	inv_temp_param = args.temp_param
	log_nc_concentration = np.log(args.nc_concentration)

	forwards = model_probs.sum_forwards(seq_len, score_vec, log_nc_concentration, inv_temp_param)  # initializing forward array

	backwards =model_probs.sum_backwards( seq_len, score_vec, log_nc_concentration,
						  inv_temp_param)  # initializing backwards array



	for i in range(seq_len):
		start_probability[i] = model_probs.prob_loc_start(i, seq_len   , score_vec, forwards, backwards,
		                                              log_nc_concentration, inv_temp_param)
		sit_probability[i] = model_probs.prob_loc_sit(i, seq_len, score_vec, forwards, backwards,
		                                          log_nc_concentration, inv_temp_param)

	np.save("kaplan_nucleusome_sit", np.exp(sit_probability))
	np.save("kaplan_nusleusome_start", np.exp(start_probability))

	plt.plot(range(seq_len), np.exp(sit_probability))
	np.save("pos_dep_model_mat",loc_dep_model_mat )





	#
	# fh = open("scores.tab",'r')
	#
	# lines = np.array(fh.read().split('\n')[1:])
	#
	# lines = lines[: len(lines)-1]
	# length = len(lines)
	# print length
	# score = np.zeros(length)
	#
	# for i in range(length):
	#
	# 	score[i] = lines[i].split('\t')[2]
	#
	# fh.close()
	#
	# np.save("segal_scores", score)
	#
	#
	# probs = range(len(score))
	#
	# seq_len =len(score)
	# forwards = np.arange(seq_len,
	# 						 dtype=np.float)  # array contains in the i place the sum off the scores of all configuration in seq[1,...,i]
	# backwards = np.arange(seq_len,
	# 						  dtype=np.float)  # array contains in the i place the sum off the scores of all configuration in seq[i,...,seq_len]
	# inv_temp_param = args.temp_param
	# log_nc_concentration = np.log(args.nc_concentration)
	# dynamic.sum_backwards(backwards, seq_len, score,log_nc_concentration, inv_temp_param)  # initializing backwards array
	# dynamic.sum_forwards(forwards, seq_len, score, log_nc_concentration, inv_temp_param)  # initializing forward array
	#
	# for i in range(seq_len):
	# 	probs[i] = dynamic.prob_loc_sit(i,seq_len, score, forwards, backwards, log_nc_concentration, inv_temp_param )
	#
	# #plt.plot(range(len(probs)), np.exp(probs))
	#
	# fh.close()
	#
	# np.save("kaplan_nucleusome_sit_probability_derived_from_segal_scores", np.exp(probs))
	#
	#
	# fh = open("occ_probs.tab",'r')
	#
	# lines = np.array(fh.read().split('\n')[1:])
	#
	# lines = lines[: len(lines)-1]
	# length = len(lines)
	# print length
	# probs = np.zeros(length)
	#
	# for i in range(length):
	#
	# 	probs[i] = lines[i].split('\t')[3]
	#
	# fh.close()
	#
	# np.save("segal_nucleusome_sit_probability_derived_from_segal_scores", probs)
	#


	plt.show()


#main(sys.argv[1:])



sep = os.sep
base_dir = system_funcs.prog_root_dir() +sep
segal_dir = base_dir +'segal_results'+sep
kaplan_dir = base_dir+'kaplan_results'+sep
segal_mats_dir =  base_dir +'segal_model_matrices'+sep
seq_dir = base_dir+'sequence_file'+sep
train_dir = base_dir + 'train_files' + sep

argv_input_gxw_file_check = 'model.py -pos_indep_model_generate_method input_file -pos_indep_model_file_method numpy -pos_indep_model_prob_output_method numpy ' \
	   '-pos_indep_model_input_file ' +segal_mats_dir+ 'segal_pos_indep_pwm.npy ' \
'-pos_indep_model_markov_order 4 -pos_indep_model_prob_output_file ' + segal_dir +'pos_indep_probs ' \
'-pos_dep_model_file_method numpy' \
	   ' -pos_dep_model_input_file ' + segal_mats_dir+ 'segal_pos_dep_pwm.npy' \
	   ' -uniform 0 -offset 0 -seq_file ' + seq_dir+'seq_input ' \
'-pos_dep_model_generate_method input_file -pos_dep_model_prob_output_method numpy -pos_dep_model_markov_order 1' \
' -pos_dep_model_prob_output_file ' + segal_dir+ 'pos_dep_probs -temp_param 1 -nc_concentration  1 -row_binding_output_method numpy' \
' -row_binding_file_name ' \
 +segal_dir+'kaplan_scores_from_segal_pwm '



argv_for_train_file_check = 'model.py -pos_indep_model_generate_method train_file_pos_indep -pos_indep_model_prob_output_method  numpy ' \
'-pos_indep_model_markov_order 4  -loc_indep_const 2 -config_vec config_vec.npy -chromosome_seq chromosome -alphabet A,C,G,T,B,D -rc_alphabet T,G,C,A,D,B -pos_indep_model_prob_output_file pos_indep_probs '  \
	   ' -uniform 0 -offset 0 -seq_file ' + seq_dir+ 'seq_input ' \
'-pos_dep_model_generate_method train_file -pos_dep_model_prob_output_method numpy -pos_dep_model_markov_order 1' \
' -pos_dep_model_prob_output_file pos_dep_probs -temp_param 1 -nc_concentration  1 -row_binding_output_method numpy' \
' -row_binding_file_name ' + kaplan_dir+ \
'kaplan_scores  -train_file ' +train_dir+ 'fasta_file1'

main(argv_for_train_file_check.split()[1:])