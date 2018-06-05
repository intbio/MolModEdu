#!/usr/bin/env python2.7
"""
This is a module which provides the analysis of atomic displacements
in molecular system, and calculates RMSF.
"""


import numpy as np
from collections import OrderedDict

__author__="Alexey Shaytan"




def atomic_displ(A_coordinates,A_indices,B_coordinates,B_indices,d_threshold=3.9,exclude_bonded=False):
	""" Function calculates atomic displacements between strucure and reference and outputs it as B-factor.

	Parameters
	----------
	A_coordinates - list of tuples with (x,y,z) of group A.
	A_indices - list of indices of atoms corresponding to the list of coordinates. This facilitates using and avoids unneeded conversion by user.
	B_coordinates - same as for A.
	B_indices - same as for A.
	d_threshold - distance threshold for contacts
	returns - list of tuples (A_ind, B_ind, distance)
	"""
	bond_threshold=1.7

	A_tree=cKDTree(A_coordinates)
	B_tree=cKDTree(B_coordinates)

	#We can calculate neighbors but it is not needed since it is automatically done by distance matrix calculation
	#neighbors=A_tree.query_ball_tree(B_tree,d_threshold)
	#print neighbors #it is a list of lists

	dok_dist=A_tree.sparse_distance_matrix(B_tree,d_threshold)
	# Now we have a sparse distance matrix in a dictionary of keys format (see scipy.sparse.dok_matrix)

	