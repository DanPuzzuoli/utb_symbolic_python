#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 12:52:53 2019

@author: dpuzzuol

Notes:
    This class stores a list of symbolic utb matrices, along with some
    functionality for working with them. This file also contains some functions
    outside of the class that could potentially be useful elsewhere.
"""

from sympy.matrices import BlockMatrix,Matrix, zeros
from sympy.core.symbol import Symbol

class Sym_UTB_Matrix:
    
    ###################
    # Constructors
    ###################
    
    def __init__(self, mat_list,spec = 'full'):
        self.mat_list = to_noncomm_matrix_list(mat_list, spec = spec)
        self.shape = [A.shape[0] for A in self.mat_list]
    
    @classmethod
    def from_concise_description(cls, shape, nc_sym, nc_idx, c_sym, c_idx):
        '''
        constructs the object from the description format output by
        the function concise_description()
        '''
        
        # need to construct a list of matrices from the concise description
        # instantiate a list of zero Matrix objects
        mat_list = []
        for d in shape:
            mat_list.append(zeros(d,d))
        
        # populate the non-commutative symbols
        for k in range(len(nc_sym)):
            sym = nc_sym[k]
            
            for mat_idx, i, j in nc_idx[k]:
                
                mat_list[mat_idx][i,j] = sym
        
        # populate hte list of commutative symbols
        for k in range(len(c_sym)):
            sym = c_sym[k]
            
            for mat_idx, i, j in c_idx[k]:
                
                mat_list[mat_idx][i,j] = sym
        
        return Sym_UTB_Matrix(mat_list)
            
    
    
    ##################
    # object methods
    ##################
    
    def exp_deriv_generators(self, params):
        '''
        constructs the generators for computing the first derivative of the
        time ordered exponential of self, with respect to params
        
        Note: I'll need to either add in stuff that keeps track of where the
        blocks go, or just have a convention that is maintained everywhere.
        
        Convention is easier, but it might be nice to also return some data
        structure with this info
        '''
        
        D,C = self.parameter_decomposition(params)
        
        
        new_mats= []
        
        for sys_idx in range(len(self.mat_list)):
            G = self.mat_list[sys_idx]
            m0 = zeros(self.shape[sys_idx])
            
            for d_idx in range(len(params)):
                
                A = C[d_idx].mat_list[sys_idx]
                new_mats.append(Matrix(BlockMatrix(
                        [[G, A], 
                         [m0, G]]
                        )))
        
        return Sym_UTB_Matrix(new_mats)
        
    
    def parameter_decomposition(self, params):
        '''
        Given a list of parameters, which should be a list of commutative
        symbols, returns a list of Sym_UTB_Matrix objects, representing
        the parameter decomp.
        
        That is, for parameters [a[0], ..., a[k]], returns a list 
        [D, M0, ..., Mk] s.t. self = M0 + a[0]*M0 + ... + a[k]*Mk
        
        Note: this function assumes that each matrix element of self
        is an affine combination with coefficients a[0], ..., a[k], which
        will be true for these control problems, but of course is not true
        in general
        
        May want to change to affine control decomp (or something)
        '''
        
        D = self.subs({a : 0 for a in params})
        
        M = [self.diff(a) for a in params]
        
        return D,M
    
    def subs(self, sub_dict):
        '''
        substitute symbols in mat_list
        '''
        return Sym_UTB_Matrix([A.subs(sub_dict) for A in self.mat_list])
    
    
    def diff(self, sym):
        '''
        differentiation with respect to a symbol
        '''
        return Sym_UTB_Matrix([A.diff(sym) for A in self.mat_list])
        
    
    def direct_sum(self,other):
        '''
        direct sum two Sym_UTB_Matrix objects
        '''
        new_list = self.mat_list.copy()
        new_list.extend(other.mat_list)
        
        return Sym_UTB_Matrix(new_list)
    
    def add(self,other):
        '''
        if the shapes are equal, add with another like object
        '''
        if self.shape == other.shape:
            return Sym_UTB_Matrix([A+B for A,B in zip(self.mat_list, other.mat_list)])
        
        
    def mult(self, other):
        '''
        if the shapes are equal, multiply with another like object
        '''
        if self.shape == other.shape:
            return Sym_UTB_Matrix([A*B for A,B in zip(self.mat_list, other.mat_list)])
    
    def power(self, k):
        '''
        multiply with itself
        '''
        return Sym_UTB_Matrix([A**k for A in self.mat_list])
    
    # check to see if these are the same matrices
    # NOTE!!!: this comparison doesn't ignore symbols
    def compare(self,other):
        return self.mat_list == other.mat_list
    
    def concise_description(self):
        '''
        Compiles a reduced description of the matrix. In particular, it
        compiles a list of every unique non-commutative and commutative symbol, 
        along with corresponding sets of indices
        '''
        
        # list of non-commutative symbols, and list of sets of indices
        # for these symbols
        non_comm_symbols = []
        non_comm_idx = []
        
        # list of commutative symbols, along with corresponding sets of
        # indices
        comm_symbols = []
        comm_idx = []
        
        # iterate through successive off-diagonals
        for off_diag in range(max(self.shape)):
            
            # iterate through each matrix
            for mat_idx in range(len(self.shape)):
                A = self.mat_list[mat_idx]
                d = A.shape[0]
                if off_diag < d:
                    for j in range(d-off_diag):
                        entry = A[j,j+off_diag]
                        
                        # check if its commutative
                        if entry.is_commutative:
                            if entry not in comm_symbols:
                                comm_symbols.append(entry)
                                comm_idx.append({(mat_idx, j,j+off_diag)})
                            else:
                                comm_idx[comm_symbols.index(entry)].add((mat_idx, j,j+off_diag))
                        # if it's not a number, assume it is a symbolic sympy
                        # expression
                        elif entry not in non_comm_symbols:
                            non_comm_symbols.append(entry)
                            non_comm_idx.append({(mat_idx, j,j+off_diag)})
                        else:
                            non_comm_idx[non_comm_symbols.index(entry)].add((mat_idx, j,j+off_diag))
        
        return self.shape, non_comm_symbols,non_comm_idx, comm_symbols,comm_idx
    
    def concise_display(self):
        '''
        Print the output of self.concise_description()
        '''
        shape, nc_symbols, nc_idx, c_symbols, c_idx = self.concise_description()
        
        print('Shape:')
        print('------')
        print(self.shape)
        
        print('Non-commutative symbols:')
        print('------------------------')
        for k in range(len(nc_symbols)):
            print('{} {}'.format(nc_symbols[k], nc_idx[k]))
        
        print('Commutative symbols:')
        print('--------------------')
        for k in range(len(c_symbols)):
            print('{} {}'.format(c_symbols[k], c_idx[k]))
        
    
    
'''
Non-class functions
'''

  
def to_noncomm_matrix_list(mat_list, spec = 'full'):
    '''
    Converts input into a list of Matrix objects with non-commutative
    symbols. Handles several cases for ease of input:
        -   if the input type is a Matrix, assume it is of the correct form
            and just needs to be put into a list
        -   if the input type is list, change behaviour based on whether
            its first entry is a Matrix, if its a single matrix specified
            using lists, or if it's a list of matrices each specified with
            lists
    '''
    if type(mat_list) == Matrix:
        return [mat_list]
    elif type(mat_list) == list:
        # if the first entry is a Matrix, assume it is a list
        # of matrices
        if type(mat_list[0]) == Matrix:
            return mat_list
        # if the first entry isn't a list, assume it is specified
        # as a single matrix in list form
        elif type(mat_list[0][0]) != list:
            return [to_noncomm_matrix(mat_list)]
        # otherwise, assume it is a list of matrices each specified
        # as a list
        else:
            return [to_noncomm_matrix(A, spec = spec) for A in mat_list]
                

        
def to_noncomm_matrix(matrix, spec = 'full'):
    '''
    Converts input into a sympy Matrix object, converting all string
    entries into non-commutative matrix symbols. Also takes input type
    Matrix, in which case it simply returns the matrix.
    '''
    if type(matrix) == Matrix:
            return matrix
    elif type(matrix) == list:
        
        if spec == 'full':
            
            d = len(matrix)
            new_mat = zeros(d,d)
            for i in range(d):
                for j in range(i,d):
                    new_mat[i,j] = to_noncomm_sym(matrix[i][j])
            
            return new_mat

def to_noncomm_sym(a):
    '''
    If the input is str, convert to a non-commutative sympy Symbol,
    otherwise return it
    '''
    
    if type(a) == Symbol:
        return a
    elif type(a) == str:
        return Symbol(a, commutative = False)
    else:
        return a