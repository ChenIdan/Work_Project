

import numpy as np



import plotly.plotly as py
import plotly.graph_objs as go


def main():
	print 1






model_mat = np.load("/home/chen/Desktop/Work_Project/segal_model_matrices/segal_pos_dep_pwm.npy")



seq="AACCTCTTTAGTCTAAGTTCAGACTAGTTGGAAGTTTGTCTAGATCTCAGATTTTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCCTAGTGCAATGGGGCTTTTTTTCCATAGTCCTCGAGAGGAGGAGACGTCAGTCCAGATATCTTTGATGTCGTGATTGGAAGGACCCTTGGCCCTCCACCCTTAGGCAGTGTATACTCTTCCATAAACGGGCTATTAGTTATGAGGTCCGAAGATTGAAAAAGGTGAGGGAACTCGGCCGAACGGGAAAGACGGACATCTAGGCAACCTGACCACGGTTGCGCGTCCGTATCAAGGTCCTCTTAATAGGCCCCCGTTACTGTTGGTCGTAGAGCCCAGAACGGGTTGGCCAGATGTGCGACAATATCGCTTAGTCGCTCTTGGGCCGCGGTGCGCTACCTTGCAGGAATTGAGACCGTCCGTTAATTTCCCTTGCATATATATTGCGTTTCTTTGACCTTTTAACCGCTCTCTTAGAAGAGAGACAGATAGCTTCTTACCGGTGCGCCACCGTAGGCAGTACGATCGCACGCCCCATGTGAACGATTGGTAAACCCAGTGTCCTGTGAGCGACAAAAGCTTAAATGGGAAATACGCGCCCATAACTTGGTGCGAATACGGGTCGTAGCAATGTTCGTCTGACTATGATCTACATATTACAGGCGGTACGTCTGCTTTGGTCAGCCTCTAATGGCTCGTAAGATAGTGCAGCCGCTGGTGATCACTCGATGACCTCGGCTCCCCATTGCAACTACGGGGATTCTTGGAGAGCCAGCTGCGTTCGCTAATGTGAGGACAGTGTAGTATTAGCAAACGATAAGTCCCGAACTGGTTGTGACCTAACGAAAAGTGAACTTCATAATACGTGCTGTCCCACGCACATGGTAGATTTGGACAAAATTGAATGGAGTCTGATCA"

seq_mat = pwm.get_seqMat(seq,1)


mid_prob = np.multiply(model_mat[73], model_mat[74])

model_mat = np.delete(model_mat,73,0)

model_mat[73] = mid_prob

prob1= prob.prob_vec(model_mat, seq_mat, 1,0)



print pwm.rc_seq(pwm.rc_seq(seq)) == seq
print len(pwm.rc_seq(seq))



seq_mat = pwm.get_seqMat(pwm.rc_seq(seq),1)

prob2= prob.prob_vec(model_mat, seq_mat, 1, 0)[::-1]


seq_graph=go.Scatter(x= range(len(prob1)), y=prob1,name = "seq_pos_probs")
rev_seq_graph = go.Scatter(x=range(len(prob2))
	                           , y=prob2,
	                       name='rev_seq_probs')


layout = go.Layout(title="symmetry_check chen method",
	                   xaxis=dict(title='sequence bases'),
	                   yaxis=dict(title='Nucleusome_probs'))

data = [seq_graph, rev_seq_graph]

fig = go.Figure(data=data, layout=layout)

py.plot(fig, sharing='public', filename='my method')



print prob1
print prob2



model_mat = np.load("/home/chen/Desktop/Work_Project/trained_models/trained_pos_dep_model.npy")

mid_prob = np.multiply(model_mat[0], model_mat[1])

model_mat = np.delete(model_mat, 0, 0)
model_mat[0] = mid_prob

seq="AACCTCTTTAGTCTAAGTTCAGACTAGTTGGAAGTTTGTCTAGATCTCAGATTTTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCCTAGTGCAATGGGGCTTTTTTTCCATAGTCCTCGAGAGGAGGAGACGTCAGTCCAGATATCTTTGATGTCGTGATTGGAAGGACCCTTGGCCCTCCACCCTTAGGCAGTGTATACTCTTCCATAAACGGGCTATTAGTTATGAGGTCCGAAGATTGAAAAAGGTGAGGGAACTCGGCCGAACGGGAAAGACGGACATCTAGGCAACCTGACCACGGTTGCGCGTCCGTATCAAGGTCCTCTTAATAGGCCCCCGTTACTGTTGGTCGTAGAGCCCAGAACGGGTTGGCCAGATGTGCGACAATATCGCTTAGTCGCTCTTGGGCCGCGGTGCGCTACCTTGCAGGAATTGAGACCGTCCGTTAATTTCCCTTGCATATATATTGCGTTTCTTTGACCTTTTAACCGCTCTCTTAGAAGAGAGACAGATAGCTTCTTACCGGTGCGCCACCGTAGGCAGTACGATCGCACGCCCCATGTGAACGATTGGTAAACCCAGTGTCCTGTGAGCGACAAAAGCTTAAATGGGAAATACGCGCCCATAACTTGGTGCGAATACGGGTCGTAGCAATGTTCGTCTGACTATGATCTACATATTACAGGCGGTACGTCTGCTTTGGTCAGCCTCTAATGGCTCGTAAGATAGTGCAGCCGCTGGTGATCACTCGATGACCTCGGCTCCCCATTGCAACTACGGGGATTCTTGGAGAGCCAGCTGCGTTCGCTAATGTGAGGACAGTGTAGTATTAGCAAACGATAAGTCCCGAACTGGTTGTGACCTAACGAAAAGTGAACTTCATAATACGTGCTGTCCCACGCACATGGTAGATTTGGACAAAATTGAATGGAGTCTGATCA"

seq_mat = pwm.get_seqMat(seq,1)
prob1= prob.prob_vec(model_mat, seq_mat, 1,0)



print pwm.rc_seq(pwm.rc_seq(seq)) == seq
print len(pwm.rc_seq(seq))



seq_mat = pwm.get_seqMat(pwm.rc_seq(seq),1)

prob2= prob.prob_vec(model_mat, seq_mat, 1, 0)[::-1]



seq_graph=go.Scatter(x= range(len(prob1)), y=prob1,name = "seq_pos_probs")
rev_seq_graph = go.Scatter(x=range(len(prob2))
	                           , y=prob2,
	                       name='rev_seq_probs')


layout = go.Layout(title="symmetry_check noam method",
	                   xaxis=dict(title='sequence bases'),
	                   yaxis=dict(title='Nucleusome_probs'))

data = [seq_graph, rev_seq_graph]

fig = go.Figure(data=data, layout=layout)

py.plot(fig, sharing='public', filename='noam method')

print prob1
print prob2


