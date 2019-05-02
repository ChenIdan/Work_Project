#!/usr/bin/python
'''
TODO: Write documentation to the method that creats the position independent model matrix
TODO: Make the length of the alphabet and the akphabet itself a parameter
TODO: Make the length of the molecule a parameter
TODO: Add possibility to symmetrize the model.
TODO: write summary of what we did so far
'''
import numpy as np
import seq_funcs
import re
from string import maketrans


def Jpwm_from_seq(pos_mats_sum, delta, mats_num, mk_ord,
                  alphabet,protein_len):  # function for creating Joint probability matrix from a list of position
    #  matrices
    Jpwm = pos_mats_sum
    alphabet_len = len(alphabet)
    i = 0




    for row in Jpwm:
            Jpwm[i] = (Jpwm[i] + (np.ones(np.size(row)) * delta*np.power(alphabet_len, protein_len - mk_ord - 1))) / \
                      (np.ones(np.size(row)) * (mats_num + delta*np.power(alphabet_len, protein_len)))
            i = i + 1


    return Jpwm

def get_Jpwm(seq_file, mk_ord, cur_mk_ord, uniform, offset, alphabet, rc_alphabet, seq_len):  # return joint probability matrix

    seq_mats_num = 0  # number of sequence matrices

    fh = open(seq_file, 'r')  # open sequence file

    alphabet_len = len(alphabet)

    nums = range(0, alphabet_len)  # locations
    loc = ''.join(str(x) for x in nums)
    trantab = maketrans(''.join(alphabet), loc)

    str_alphabet = ''.join(alphabet)

    regex = re.compile(r'[ACGT]+')

    summing_matrix = np.zeros((seq_len, seq_len - mk_ord))

    for i in range(0,mk_ord+1):
        np.fill_diagonal(summing_matrix[::][i:], np.power(len(alphabet),mk_ord - i))


    cur_mats_sum= np.zeros((seq_len - mk_ord, np.power(len(alphabet), mk_ord+1)))
    fh.seek(0)





    while fh:

        lines = fh.readlines(10000000)

        if lines == []:
            break

        protein_lines = filter(regex.search,lines)

        protein_lines = ''.join(protein_lines)

        protein_lines = ''.join(protein_lines.split())

        protein_lines = protein_lines.translate(trantab)

        protein_lines = np.array(list(protein_lines))

        protein_lines = protein_lines.astype(np.float64)

        protein_lines = np.reshape(protein_lines,(len(protein_lines)/ seq_len,  seq_len))

        protein_kmers = np.dot(protein_lines, summing_matrix)

        N = np.power(len(alphabet), mk_ord+1)

        protein_kmers = np.transpose(protein_kmers)
        id = protein_kmers + (N * np.arange(protein_kmers.shape[0], dtype = np.int64))[:,None]

        kmers_counts = np.bincount((id.ravel()).astype(np.int64),minlength=N*protein_kmers.shape[0]).reshape(-1,N)

        cur_mats_sum = cur_mats_sum + kmers_counts


        seq_mats_num = seq_mats_num + (4*(offset > 0)  + 2)*len(protein_lines )  # count the sequence matrices










    #for line in fh:  # loops though file lines
     #   if line[0] in alphabet:

      #      if (offset > 0):
       #         cur_mat = seq_funcs.get_seqMat(line.splitlines()[0], mk_ord, uniform, offset,
             #                                     alphabet)  # translate sequence lines to sequence matrix
        #        cur_mats_sum = cur_mats_sum + cur_mat  # add current sequence matrix to sum

         #       cur_mat = seq_funcs.get_seqMat(line.splitlines()[0], mk_ord, uniform, -offset,
            #                                      alphabet)  # translate sequence lines to sequence matrix
          #      cur_mats_sum = cur_mats_sum + cur_mat  # add current sequence matrix to sum


           # cur_mat = seq_funcs.get_seqMat(line.splitlines()[0], mk_ord, uniform, 0,
           #                                   alphabet)  # translate sequence lines to sequence matrix

            #cur_mats_sum = cur_mats_sum + cur_mat  # add current sequence matrix to sum

            #seq_mats_num = seq_mats_num + (offset > 0) * 4 + 2  # count the sequence matrices



    if cur_mats_sum == []:
        exit("there are no sequence matrices in you train file\n")

    cur_mat_sum = np.transpose(cur_mats_sum)

    cur_mats_sum = cur_mats_sum[uniform: len(cur_mats_sum)-uniform] + \
                   (offset>0)*(cur_mats_sum[uniform + offset: len(cur_mats_sum)-uniform+offset] + cur_mats_sum[uniform - offset: len(cur_mats_sum)-uniform-offset])

    delta = 1.0 / np.power(float(alphabet_len), seq_len -2*uniform - cur_mk_ord - 1)  # we convert to float to prevent division by zero
    fh.close()

    inv_mat = seq_funcs.rc_mat(mk_ord, alphabet, rc_alphabet)

    mats_sum = cur_mats_sum + np.dot(cur_mats_sum,inv_mat)[::-1]



    joint_mat = Jpwm_from_seq(mats_sum, delta, seq_mats_num, mk_ord, alphabet, seq_len)


    return joint_mat




def pos_dep_con_mat(seq_file, mk_ord, uniform, offset, alphabet,
                    rc_alphabet, seq_len):  # calculating model probabilities as a position dependent probabilities
    if mk_ord == 0:  # return joint probability matrix if markov order is zero
        return get_Jpwm(seq_file, mk_ord, mk_ord, uniform, offset, alphabet, rc_alphabet)

    # create conditional probability matrix from joint probability matrices from current markov order
    #  (mk_ord) and previous markov order (mk_ord-1)
    alphabet_len = len(alphabet)
    joint_mat1 = np.array(get_Jpwm(seq_file, mk_ord - 1, mk_ord, uniform, offset, alphabet,
                                   rc_alphabet, seq_len))  # get joint probability matrix of previous order

    joint_mat2 = np.array(get_Jpwm(seq_file, mk_ord, mk_ord, uniform, offset, alphabet,
                                   rc_alphabet, seq_len))  # get joint probability matrix from current markov
    #  order

    joint_vec = np.array(joint_mat2[0])

    np.save("Joint mat1", joint_mat1)
    np.save("Joint mat2", joint_mat2)

    joint_mat1 = np.delete(joint_mat1, joint_mat1.shape[0] - 1, axis=0)  # delete last row

    joint_mat2 = np.delete(joint_mat2, 0, axis=0)

    joint_mat1 = np.delete(joint_mat1, 0 ,axis=0)


    joint_mat1 = np.repeat(joint_mat1, alphabet_len, axis=1)  # repeat every coloumn 4 time



    con_mat = np.divide(joint_mat2, joint_mat1)

    model_mat = np.vstack((joint_vec, con_mat))


    return model_mat  # return conditional probability matrix


def kill_spare_lines(model_mat, mk_order):
    for i in range(0, mk_order+1):
        first_row = np.multiply(model_mat[0], model_mat[1])
        model_mat = np.delete(model_mat, 0, 0)
        model_mat[0] = first_row

    return model_mat