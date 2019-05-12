"""This module contains functions that translate the gxw format to a pwm matrix.
this matrix can be location dependent, and location independent"""

import xml.etree.ElementTree as ET



import numpy as np
from numpy.core.multiarray import ndarray
from string import maketrans
import system_funcs
import loc_dep_model




def get_position_independent_pwm(gxw_path, uniform):
	"""

	:param gxw_path:
	:return: position independent pwm matrix, created from the first weight matrix element in the gxw file
	(the name of that element is "nucleusome background")
	"""


	fh = open(gxw_path)
	if not fh:
		print "cant open gxw file in path" + gxw_path
		return -1


	#read gxw file into an xml data structure
	gxw_str = fh.read()
	fh.close()
	gxw_tree = ET.fromstring(gxw_str)

	#find the position independent probability matrix in the gxw file
	position_independent_mat_element = gxw_tree.find('WeightMatrix')

	#get the matrix markov order
	mk_order =int(position_independent_mat_element.attrib['Order'])
	#get the number of sub elements in the gxw element that represents the pwm matrix. the number is the sum of all powers of 4 that are smaller from the matrix markov order
	list_len = (np.power(4, mk_order + 1) - 1) / 3

	#get the list of each sub element weights
	weights_list = map(weights_to_arr, position_independent_mat_element[0:list_len])  # type: List[ndarray]

	#
	#   the next loop will translate the weights of the all the sub elements in the weight list
	#   (referenced by the pointer @weights_list) into a position independent pwm matrix
	#  each iteration will create a new row in the pwm matrix, starting from the first row
	#


	idx = 0  # type: int # use: idx will be used to reference the start of the sub_list in the weight list that it's members will be used



	pwm_mat = np.ones(np.power(4, mk_order+1))
	for cur_mk_ord in range(mk_order+1):

		# <editor-fold desc="get the current row vector">

		# get the weights from the weights list that are required to create the probability pwm row
		cur_row= np.array(weights_list[idx: idx + np.power(4,cur_mk_ord)])
		#
		# make the weights list to a one dimentional array/
		# for example, if the weight list is [[1,0,0,0],[1,0,0,0]] after this code line it will
		# be [1,0,0,0,1,0,0,0]
		#
		cur_row = np.reshape(cur_row , -1)
		#
		# repeats the row so it will fill the whole pwm row (because it's markov order is smaller from the pwm matrix
		# general markov order
		#
		cur_row = np.repeat(cur_row, np.power(4,mk_order - cur_mk_ord))
		# </editor-fold>

		#promots the idx
		idx = idx + np.power(4,cur_mk_ord)

		#becuase this pwm matrix is not position dependent, the last row will be repeated until the end of the matrix
		#the last row is, of course, the one which it's markov order in equall to the model markov order
		if cur_mk_ord == mk_order:
			cur_row= np.tile(cur_row, (146 - mk_order   - 2*uniform,1))

		#stack vertically the current row into the pwm matrix
		pwm_mat = np.vstack((pwm_mat, cur_row))



	model_mat = np.delete(pwm_mat, 0, 0)

	model_mat =loc_dep_model.kill_spare_lines(model_mat,mk_order)




	upper_uniform = np.ones((uniform, np.power(4,mk_order+1)))

	lower_uniform = np.ones((uniform, np.power(4,mk_order+1)))

	model_mat = np.vstack(( lower_uniform,model_mat))

	model_mat = np.vstack(( model_mat,upper_uniform ))

	return model_mat








def weights_to_arr(elements):
	weights = elements.attrib['Weights']
	array = np.array(weights.split(';'))
	return np.asfarray(array, np.int)

def get_position_dependent_pwm(gxw_path, mk_order):
	"""

	:param gxw_path:  path of the gxw file from which we will build the position dependent matrix
	:param mk_order: markov order of the position dependent matrix
	:return: position dependent matrix, translated from the gxw file
	"""
	fh = open(gxw_path)
	if not fh:
		print "cant open gxw file in path" + gxw_path
		return -1

	gxw_str = fh.read() #read xml file
	fh.close() #close pointer to xml file
	gxw_tree = ET.fromstring(gxw_str) #read the xml file (gxw file from now on) into a data structure

	low_pwm=np.zeros((74, np.power(4, mk_order+1)))  # type: ndarray # use: the lower half pwm matrix

	loc = mk_order # type: int # use: the location of the current row written on the pwm matrix

	#
	# loop for reading the lower half of the pwm matrix from the gxw file into the lower half of the low_pwm matrix.
	# this loop reades only the rows with markov order that is equals to the pwm markov order
	#


	# <editor-fold desc= loop for reading the lower half of the pwm matrix from the gxw file into the lower half of the low_pwm matrix.
	# 	 this loop reades only the rows with markov order that is equals to the pwm markov order">

	for sub_mat in gxw_tree.findall('WeightMatrix')[1+mk_order:75]:
		#
		# all the sub_mat elements in the gxw file representing the rows of the lower half of the pwm matrix
		# the sub elements weights will be normalized with the the weights of the first sub element
			for child in sub_mat:
				#get the first sub element weights
				if int(child.attrib['Markov']) == mk_order-1:
					norm_weights = np.array(child.attrib['Weights'].split(';'))
					norm_weights = np.asfarray(norm_weights, np.float)

				# get every four conditional probability values, for every dna base pare combination of markov order length
				else:

					#get the weights
					prob_weights = np.array(child.attrib['Weights'].split(';'))
					prob_weights = np.asfarray(prob_weights, np.float)

					#normalizing the weights
					prob_weights = norm_weights*prob_weights/(np.sum(norm_weights*prob_weights))

					#get the parents of the conditional probability
					parents = np.array(child.attrib['Parents'].split(';'))
					parents = ''.join(parents)
					parents = int(parents , base=4)

					#writes the conditional probability values into the pwm
					low_pwm[loc, 4*parents :4*parents +4] =+ prob_weights

			#promots location to the next row
			loc= loc+1
	# </editor-fold>



	# <editor-fold desc="this part is used for read the pwm rows that their probability markov order is smaller from the pwm markov order">

	loc = 0 # the location variable will be used to promot the current row on the pwm
	#
	#
	# loop on all the elements representing the pwm matrix rows in the gxw file
	# each rows will be normalized
	#
	for sub_mat in gxw_tree.findall('WeightMatrix')[1:1 + mk_order]:
		sub_mat_order = int(sub_mat.attrib['Order']) # get current row conditional probabilities markov order
		loc_prob_row = np.zeros(np.power(4, sub_mat_order)) # row of the conditional probabilities
		if sub_mat_order != 0:
			for child in sub_mat:
				if int(child.attrib['Markov']) == sub_mat_order - 1:

					#get the weights for normalizing the conditional probabilities
					norm_weights = np.array(child.attrib['Weights'].split(';'))
					norm_weights = np.asfarray(norm_weights, np.float)


				else:
					# <editor-fold desc="region get every four conditional probability values, for every dna base pare combination of markov order">

					#get the weights of the conditional probabilities
					prob_weights = np.array(child.attrib['Weights'].split(';'))
					prob_weights = np.asfarray(prob_weights, np.float)

					# normalize the weights with the norm weights
					prob_weights = (norm_weights*prob_weights) / (np.sum(norm_weights * prob_weights))

					# get the parents of the conditional probability
					parents = np.array(child.attrib['Parents'].split(';'))
					parents = ''.join(parents)
					parents = int(parents, base=4)

					#writes the conditional probabilities
					loc_prob_row[loc, 4 * parents: 4 * parents + 4] = + prob_weights
	# </editor-fold>
			# endregion

		# region if the markov order of the rows is zero, the weights will be normalize in the classical way
		else:
			for child in sub_mat:

				#get the weights
				prob_weights = np.array(child.attrib['Weights'].split(';'))
				prob_weights = np.asfarray(prob_weights, np.float)

				#normalize the weights
				prob_weights = prob_weights / (np.sum(prob_weights))
		# endregion


		#repeat the probabilities on the pwm matrix so they will fill all the row

		low_pwm[loc] = + np.repeat(prob_weights, np.power(4, mk_order - sub_mat_order))

		#promots row location
		loc = loc +1









	# get the mirror image of the low pwm matrix, and flip the direction of sequence reading, so our model matrix (the pwm) will be symmetric according to the

	np.save("low_pwm_example", low_pwm)

	return symmetrize_model(low_pwm)




def symmetrize_model(low_model_mat):
	"""

	:param low_model_mat:
	:return returns the model matrix as a symmetric model
	"""


	mk_order = int(np.log2(len(low_model_mat[0]))/2) -1 #get markov order of the model

	last_row = np.copy(low_model_mat[mk_order])

	low_model_mat = loc_dep_model.kill_spare_lines(low_model_mat)
	up_model_mat = np.flipud(low_model_mat).copy()  # get the matrix flipped



	up_model_mat = np.delete(up_model_mat, len(up_model_mat) - 1, axis=0) #delete last row, so the middle of the model wont be double and wont be reversed
	if (mk_order + 1) % 2 == 0:

		up_model_mat = np.vstack((up_model_mat, last_row))


	# top_row = low_model_mat[0]
	#
	# for i in range(mk_order):
	# 	top_row = np.multiply(low_model_mat[0], low_model_mat[1])
	# 	low_model_mat = np.delete(low_model_mat,0,0)
	# 	low_model_mat[0] = top_row


	for row in up_model_mat :
		for i in range(len(row)):
			# get current index on the current row to base 4 string representation with the right padding (this is the last parametr)
			num_str = np.base_repr(i,4, padding= mk_order+1 - len(str(np.base_repr(i,4,0))))
			#reverse the base 4 string representation
			num_str = num_str[::-1]
			DNA = '0123'  # dna proteins
			RVERSE_DNA = '3210'  # locations
			trantab = maketrans(DNA, RVERSE_DNA)
			num_str = num_str.translate(trantab)  # convert the sequence into a string of numbers, according to the dictionary of tranbtab

			if i>int(num_str,base=4): #this condition exists to prevent double swapping (AC to CA and CA to AC, for example)
				#swap sequence reading direction
				tmp_num = row[int(num_str, base=4)]
				row[int(num_str, base= 4)] = row[i]
				row[i] = tmp_num





	model_mat = np.vstack((up_model_mat, low_model_mat)) #now vertically stack the lower and the upper halves of the pwm model matrix and returns the result



	return model_mat


cur_dir = system_funcs.prog_root_dir()




model=get_position_independent_pwm("nucleosome_model_1208.gxw",0)
np.save("segal_pos_indep_pwm*", model)