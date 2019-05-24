#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 12:52:26 2019

@author: dpuzzuol

Notes:
    This class represents particular direct sums of subspaces of UTB matrices.
    Specifically, the data stored by this class is:
        - matrix_blocks -   a set of sets, where each lower level set represents
                            a repeated block that could take the form of any 
                            matrix, with the set containing all of the indices 
                            for where that block is repeated.
        - scalar_blocks -   the same as matrix_blocks, except meant to represent
                            blocks that are scalars, i.e. some multiple of the
                            identity.
        - zero_blocks -     the set of all zero blocks
        - shape -           a list containing the block dimension of each
                            direct sum
                            
        
    This class contains functionality for working with and manipulating these
    sets, and works with 
    
    
    NOTE: index lists are sets of frozensets
    NOTE2: I've changed it so that nc_blocks and c_blocks are stored as lists,
    which are then later turned into a set when set-like operations are needed
"""

from Sym_UTB_Matrix import Sym_UTB_Matrix, to_noncomm_sym
from sympy import Symbol

class UTB_Subspace:
    
    ###################
    # Constructors
    ###################
    
    def __init__(self, nc_blocks,c_blocks,z_blocks, shape):
        self.nc_blocks = nc_blocks
        self.c_blocks = c_blocks
        self.z_blocks = z_blocks
        self.shape = shape
    
    @classmethod
    def from_Sym_UTB_Matrix(cls, matrix):
        '''
        instantiate the object from a Sym_UTB_Matrix object. Basically,
        needs to convert it into one that 'forgets' the symbols
        '''
        shape, nc_sym,nc_idx,c_sym,c_idx = matrix.concise_description()
        
        
        # extract the nc_blocks, converting each inner set into a frozenset
        nc_blocks = [frozenset(s) for s in nc_idx]
        
        # check if 0 is in the list of commutative symbols. If it is,
        # pop the indices of the zeros (which is already a set of tuples),
        # otherwise create an empty set
        if 0 in c_sym:
            z_blocks = frozenset(c_idx.pop(c_sym.index(0)))
        else:
            z_blocks = frozenset()
        
        # get the c_blocks, converting inner set into a frozen set
        c_blocks = [frozenset(s) for s in c_idx]
        
        return cls(nc_blocks, c_blocks, z_blocks, shape)
    
    @classmethod
    def from_matrix(cls,matrix):
        '''
        instantiate it from an input that is accepted by the default 
        constructor of Sym_UTB_Matrix
        '''
        
        return cls.from_Sym_UTB_Matrix(Sym_UTB_Matrix(matrix))
    
    
    ##################
    # object methods
    ##################
    
    def issuperset(self,other):
        '''
        checks if the subspace is a superset of another subspace. We do the
        following checks:
            - check shape, if not equal, return false
            - check zeros, if the zeros of self are not a subset of zeros
            of other, return false
            - every self nc_block index set must be a subset of an nc_block
            index set in other
            - every self c_symbol index set must be contained in either an
            nc_block index set of other, or c_block index set of other 
        '''
        
        # first, if they are not the same shape, return false immediately.
        if self.shape != other.shape:
            return False
        
        #check zeros
        if not self.z_blocks.issubset(other.z_blocks):
            return False
        
        # first, the nc blocks 
        for snc_block in set(self.nc_blocks):
            found= False
            
            # for each nc_block in self, check if its contained in some blcok in
            # other
            for o_block in set.union(set(other.nc_blocks),set(other.c_blocks), {other.z_blocks}): 
                
                # if it is found to be a subset, set found to True and stop
                # looping
                if snc_block.issubset(o_block):
                    found = True
                    break
            
            # if no containing set was found, return False
            if found != True:
                return False
        
        # next, the c blocks
        for sc_block in set(self.c_blocks):
            found = False
            
            # this time, don't check nc_blocks of other
            for o_block in set.union(set(other.c_blocks), {other.z_blocks}):
                
                if sc_block.issubset(o_block):
                    found = True
                    break
            
            if found != True:
                return False
        
        # if it made it to the end without returning False, return True
        return True
        
    def issubset(self,other):
        return other.issuperset(self)
    
    def isequal(self,other):
        return self.issubset(other) and self.issuperset(other)
    
    
    #############################
    # Need to test this!!!
    #############################
    def subspace_add(self,other):
        srep = self.get_subspace_member(nc_symbol = 'A', c_symbol = 'a')
        orep = other.get_subspace_member(nc_symbol = 'B', c_symbol = 'b')
        
        return UTB_Subspace.from_Sym_UTB_Matrix(srep.add(orep))
    
    def subspace_product(self,other):
        srep = self.get_subspace_member(nc_symbol = 'A', c_symbol = 'a')
        orep = other.get_subspace_member(nc_symbol = 'B', c_symbol = 'b')
        
        return UTB_Subspace.from_Sym_UTB_Matrix(srep.mult(orep))
    
    def containing_algebra(self):
        '''
        returns a UTB_subspace that is an algebra containing self.
        '''
        
        
        # first candidate
        new_sub = self
        is_algebra = False
        
        while not is_algebra:
            curr_sub = new_sub
            # define a new subspace by multiplying the current with itself
            # and then adding it to itself
            new_sub = curr_sub.subspace_add(curr_sub.subspace_product(curr_sub))
            
            is_algebra = new_sub.isequal(curr_sub)
        
        return new_sub
    
    def is_algebra(self):
        alg = self.containing_algebra()
        
        return self.isequal(alg)
        
    
    def get_subspace_member(self, nc_symbol = 'A', c_symbol = 'a'):
        '''
        Generate a symbolic representative from the list of symbols.
        This uses the from_concise_description() constructor of Sym_UTB_Matrix.
        
        This also works if nc_symbol and c_symbol are given as lists, giving
        the option of specifying every symbol independently
        
        Note: there is actually some weirdness here in how everything has been
        defined. The issue stems from storing the sets of indices in a set.
        As a result, the ordering is lost, so if you specify a UTB_Subspace
        from a Sym_UTB_Matrix, then convert back (using the same list of
        symbols), you may actually get a different Sym_UTB_Matrix object (with
        the symbols permuted).
        '''
        
        # get the list of noncommutative symbols
        if type(nc_symbol) == list:
            nc_sym = [to_noncomm_sym(sym) for sym in nc_symbol]
            
        elif type(nc_symbol) == str:
            nc_sym = [Symbol('{}{}'.format(nc_symbol, str(i)), commutative = False) for i in range(len(self.nc_blocks))] 
        
        # get the list of commutative symbols
        if type(c_symbol) == list:
            c_sym = c_symbol
        elif type(c_symbol) == str:
            c_sym = [Symbol('{}{}'.format(c_symbol, str(i)), commutative = True) for i in range(len(self.c_blocks))] 
        
        return Sym_UTB_Matrix.from_concise_description(self.shape, nc_sym, list(self.nc_blocks), c_sym, list(self.c_blocks))
        
    
    
    def display(self):
        '''
        print a description of the object
        '''
        
        print('Unique non-commutative blocks:')
        print('------------------------------')
        k=1
        for s in self.nc_blocks:
            print('{}. {}'.format(str(k),s))
            k+=1
        
        print('Unique commutative blocks:')
        print('--------------------------')
        k=1
        for s in self.c_blocks:
            print('{}. {}'.format(str(k),s))
            k+=1
        
        print('Zero blocks:')
        print('------------')
        print(self.z_blocks)
        
    ############
    #Old stuff I don't want to get rid of yet
    ############
    
    """
    This is old, realized I could do it way more simply using symbol stuff
    
    def subspace_sum(self,other):
        '''
        Get the subspace corresponding to the addition of two subspaces. For
        now assumes they are of the same shape
        '''
        
        # the new zero blocks are just the intersection of the old ones
        new_z_blocks = self.z_blocks.intersection(other.z_blocks)
        
        # construct the new nc_blocks. These arise from intersections of
        # nc_blocks of one subspace with intersections of any other type of
        # block from the other
        int1 = index_set_intersection(self.nc_blocks, set.union(other.nc_blocks, other.c_blocks, {other.z_blocks}))
        int2 = index_set_intersection(other.nc_blocks, set.union(self.nc_blocks, self.c_blocks, {self.z_blocks}))
        new_nc_blocks = set.union(int1,int2)
        
        #construct new z_blocks. these arise from intersections of 
        # c_blocks of one subspace with intersections of c_blocks or z_blocks
        # of the other subspace
        int1 = index_set_intersection(self.c_blocks, set.union(other.c_blocks, {other.z_blocks}))
        int2 = index_set_intersection(other.c_blocks, set.union(self.c_blocks, {self.z_blocks}))
        new_c_blocks = set.union(int1,int2)
        
        return UTB_Subspace(new_nc_blocks, new_c_blocks, new_z_blocks, self.shape)
    """

def index_set_intersection(set1, set2):
    '''
    Returns the set whose elements are a.intersection(b) for all a in set1 and
    b in set2
    '''
    
    newset = set()
    
    for a in set1:
        for b in set2:
            c = a.intersection(b)
            
            if len(c) > 0:
                newset.add(c)
    
    return newset