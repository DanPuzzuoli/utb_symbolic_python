#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 16:40:36 2019

@author: dpuzzuol


some functions to test out construction of 
"""

from Sym_UTB_Matrix import Sym_UTB_Matrix

def basic_dyson(G, Alist):
    '''
    Given a generator symbol G, and a list of A symbols, return a Sym_UTB_Matrix
    object representing a generator to compute Dyson term for the A list
    
    Need to construct this as a list of list, not a matrix, so that when its
    fed into Sym_UTB_Matrix, all of the type checking happens
    '''
    
    #if it is not a list of lists, make it so
    if type(Alist) != list:
        Alist = [[Alist]]
    elif type(Alist[0]) != list:
        Alist = [Alist]
    
    # construct the mat_list
    mat_list = [UTB_from_diag([[G]*(len(A)+1), A]) for A in Alist ]
    
    
    return Sym_UTB_Matrix(mat_list)
    

def UTB_from_diag(diag):
    """
    constructs a list of lists in full form, representing a symbolic upper
    triangular block matrix from a specification of its diagonals
    """
    
    #make sure it is sufficiently populated
    diags = diag.copy()
    
    # if there are less diagonals specified than there are dimensions,
    # fill it out
    while len(diags[-1]) > 1:
        diags.append( [0 for k in range(len(diags[-1]) -1)] )
    
    return UTB_from_rows(triangle_list_transpose(diags))

def UTB_from_rows(rows):
    """
    Constructs a list of lists in full form, representing a symbolic upper 
    triangular block matrix from a specification of itsrows.
    
    E.g. for X,Y,Z symbols or strings, maps [[X,Y], [Z]] to
    [[X,Y], [0, Z]]
    
    Parameters
    ----------
    rows :
        a list of lists containing strings, symbols, or numbers

    """
    
    # get the length of the first row, which is the longest
    d = len(rows[0])
    
    mat = [[0 for k in range(d)] for i in range(d) ]
    
    # loop through the rows
    for i in range(len(rows)):
        # loop through the columns
        for j in range(len(rows[i])):
            mat[i][i+j] = rows[i][j]
    
    return mat


def triangle_list_transpose(triangle):
    """ 
    Performs transposition on a triangular list. The purpose here is to
    transform the specification of an upper triangular block matrix in terms of
    the rows into one in terms of the diagonals, and vice versa
    
    Input: 
        a list of lists with `triangular shape', i.e. each successive list 
        has length one less than the previous list, and the last list has
        length 1
    Output:
        the transpose of the list
        e.g. [['a', 'b'], ['c']] -> [['a', 'c'], ['b']]
    """
    rows = []
    for i in range(0, len(triangle)):
        rows.append([x[i] for x in triangle[0: len(triangle)-i]])
    
    return rows