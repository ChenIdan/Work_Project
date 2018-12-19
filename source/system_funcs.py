
import os


def prog_root_dir():
	'''

	:return: program root directory in current system
	'''
	sep = os.sep

	cur_abs_path = os.getcwd().split(sep)[:-1]
	return sep.join(cur_abs_path)
