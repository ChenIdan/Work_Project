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


def Jpwm_from_seq(pos_mats_sum, const, mats_num, mk_ord,
                  alphabet,protein_len):  # function for creating Joint probability matrix from a list of position
    #  matrices
    Jpwm = pos_mats_sum
    alphabet_len = len(alphabet)
    i = 0

    delta = const

    for row in Jpwm:
            Jpwm[i] = (Jpwm[i] + (np.ones(np.size(row)) * delta*np.power(alphabet_len, protein_len - mk_ord - 1))) / \
                      (np.ones(np.size(row)) * (mats_num + delta*np.power(alphabet_len, protein_len)))
            i = i + 1

    return Jpwm

def get_Jpwm(seq_file, mk_ord, const, uniform, offset, alphabet, rc_alphabet):  # return joint probability matrix

    seq_mats_num = 0  # number of sequence matrices

    fh = open(seq_file, 'r')  # open sequence file
    cur_mats_sum = []
    cur_mat = []

    alphabet_len = len(alphabet)

    for line in fh:  # loops though file lines
        if line[0] in alphabet:
            cur_mats_sum = seq_funcs.get_mk_seqMat(line.splitlines()[0], mk_ord, uniform, offset,
                                                   alphabet)  # translate sequence line to sequence matrix
            break

    cur_mats_sum.fill(0)
    fh.seek(0)

    for line in fh:  # loops though file lines
        if line[0] in alphabet:

            if (offset > 0):
                cur_mat = seq_funcs.get_mk_seqMat(line.splitlines()[0], mk_ord, uniform, offset,
                                                  alphabet)  # translate sequence lines to sequence matrix
                cur_mats_sum = cur_mats_sum + cur_mat  # add current sequence matrix to sum

                cur_rc_mat = seq_funcs.get_mk_seqMat(seq_funcs.rc_seq(line.splitlines()[0], alphabet, rc_alphabet),
                                                     mk_ord, uniform, offset,
                                                     alphabet)  # get sequence matrix from the reverse compliment sequence
                cur_mats_sum = cur_mats_sum + cur_rc_mat  # add current sequence matrix to sum

                cur_mat = seq_funcs.get_mk_seqMat(line.splitlines()[0], mk_ord, uniform, -offset,
                                                  alphabet)  # translate sequence lines to sequence matrix
                cur_mats_sum = cur_mats_sum + cur_mat  # add current sequence matrix to sum

                cur_rc_mat = seq_funcs.get_mk_seqMat(seq_funcs.rc_seq(line.splitlines()[0], alphabet, rc_alphabet),
                                                     mk_ord,
                                                     uniform, -offset,
                                                     alphabet)  # get sequence matrix from the reverse compliment sequence
                cur_mats_sum = cur_mats_sum + cur_rc_mat  # add current sequence matrix to sum

            cur_mat = seq_funcs.get_mk_seqMat(line.splitlines()[0], mk_ord, uniform, 0,
                                              alphabet)  # translate sequence lines to sequence matrix

            cur_mats_sum = cur_mats_sum + cur_mat  # add current sequence matrix to sum

            cur_rc_mat = seq_funcs.get_mk_seqMat(seq_funcs.rc_seq(line.splitlines()[0], alphabet, rc_alphabet), mk_ord,
                                                 uniform, 0,
                                                 alphabet)  # get sequence matrix from the reverse compliment sequence

            cur_mats_sum = cur_mats_sum + cur_rc_mat  # add current sequence matrix to sum

            seq_mats_num = seq_mats_num + (offset > 0) * 4 + 2  # count the sequence matrices
        else:
            a= re.search('[0-9]+', line).group()
            if a != '1':
               break

    if cur_mats_sum == [] or cur_mat == []:
        exit("there are no sequence matrices in you train file\n")

    seq_len = np.shape(cur_mat)[0]
    const = 1.0 / np.power(float(alphabet_len), seq_len - const - 1)  # we convert to float to prevent division by zero
    fh.close()

    return Jpwm_from_seq(cur_mats_sum, const, seq_mats_num, mk_ord, alphabet, seq_len)


def pos_dep_con_mat(seq_file, mk_ord, uniform, offset, alphabet,
                    rc_alphabet):  # calculating model probabilities as a position dependent probabilities
    if mk_ord == 0:  # return joint probability matrix if markov order is zero
        return get_Jpwm(seq_file, mk_ord, mk_ord, uniform, offset, alphabet, rc_alphabet)

    # create conditional probability matrix from joint probability matrices from current markov order
    #  (mk_ord) and previous markov order (mk_ord-1)
    alphabet_len = len(alphabet)
    joint_mat1 = np.array(get_Jpwm(seq_file, mk_ord - 1, mk_ord, uniform, offset, alphabet,
                                   rc_alphabet))  # get joint probability matrix of previous order

    joint_mat2 = np.array(get_Jpwm(seq_file, mk_ord, mk_ord, uniform, offset, alphabet,
                                   rc_alphabet))  # get joint probability matrix from current markov
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
    for i in range(0, mk_order):
        first_row = np.multiply(model_mat[0], model_mat[1])
        model_mat = np.delete(model_mat, 0, 0)
        model_mat[0] = first_row

    return model_mat
