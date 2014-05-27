# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 15:03:32 2013

@author: Chih-Hao Chen
"""

import numpy as np

def Cal(dim, disp, Coor, arrM, matrix, Coor_matrix):
    zeroM = np.zeros((dim+1, dim))
    non = np.zeros((dim+1, dim+1))
    P1st = np.ones((arrM.BoudBef.shape[0], 1))
    Pmatrix = np.concatenate((P1st, arrM.BoudBef), axis=1)
    PmatrixT = np.transpose(Pmatrix)
    MmatrixUp = np.concatenate((matrix, Pmatrix), axis=1)
    MmatrixDn = np.concatenate((PmatrixT, non), axis=1)
    Mmatrix = np.concatenate((MmatrixUp, MmatrixDn), axis=0)
    dispCom = np.concatenate((disp, zeroM), axis=0)
    #The interpolant matrix obtained vis numpy's function of solve    
    Coeff = np.linalg.solve(Mmatrix, dispCom)
    
    Pmat1col = np.ones((Coor.shape[0], 1))
    Tmatrix = np.concatenate((Coor_matrix, Pmat1col, Coor[:, 1:dim+1]), axis=1)
    return Tmatrix.dot(Coeff)
    
    