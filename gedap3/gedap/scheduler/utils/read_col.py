"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
import numpy as np

def read_col (file_name, sep=' ', i_head=0):
	f = open(file_name, 'r')
	i=0
	while i<i_head:
		f.readline()
		i+=1

	head = f.readline()
	tags = head.split(sep)
	results = {}
	nb_col = len(tags)
	for i_col in range(nb_col):
		tags[i_col] = tags[i_col].replace('\n', '')
		results[tags[i_col]] = []

	for line in f.readlines():
		line = line.split(sep)
		for i_col in range(nb_col):
			try:
				results[tags[i_col]].append(int(line[i_col]))
			except:
				if line[i_col] == 'N/A': results[tags[i_col]].append(np.nan)
				else:			    	 results[tags[i_col]].append(line[i_col])
	return results
