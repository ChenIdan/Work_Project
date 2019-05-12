# function for calculating the probability that the nucleosome will sit on a specific dna sequence

import numpy as np

def Pfunc(PWM, DNA, uniform, mk_order):

    uniform_l = 0.25*np.ones(((uniform - mk_order+1)*((uniform - mk_order) > 0), DNA.get_shape()[1]))
    uniform_r = 0.25*np.ones((uniform, DNA.get_shape()[1]))

    if uniform - mk_order > 0:

        uniform_l[0] = np.power(0.25, mk_order) * uniform_l[0]

    if uniform > 0:
        PWM = np.vstack([np.vstack([uniform_l, PWM]), uniform_r]) # add ones filled rows on the upper and lower part of the matrix
    product = DNA*np.transpose(PWM)
    return np.transpose(product)






def creatPlist(PPM, mk_ord):  # returns probability list. get the probability matrix from Pfunc as input

    plist = []  # list of probabilities
    tmp=PPM
    ppm_width = PPM.shape[1]  # get probability matrix width
    ppm_len = PPM.shape[0]  # get probability matrix length
    plist_len = abs(ppm_width-ppm_len+1) # get plist len, by counting number of diagonals in ppm
    new_width = ppm_width+1
    new_ppm = np.resize(PPM, (ppm_len, new_width))  # get new probability matrix, with with all diagonals as cols

    vecs_mat = new_ppm[0:ppm_len, 0:plist_len]

    return np.log(np.prod(vecs_mat, axis=0))



def score_vec(PWM, DNA, mk_ord,uniform):
    matrix = Pfunc(PWM, DNA,uniform,mk_ord)
    return creatPlist(matrix, mk_ord)


# calculates log(pos_dep_plist/indep_pos_plist) for every sub sequence



