# -*- coding: utf-8 -*-
"""
Created on Sat Jul 06 15:41:50 2013

@author: Chih-Hao Chen
"""
import numpy as np

def displacement(arrayM,arrayT,tranX,tranY,pivot_x, pivot_y, rotation_angle):

    ''' Old version of displacement by finding the leftmost point
    if np.argmin(arrayM.InBoud[:, 0]) == np.argmin(arrayM.InBoud[:, 1]):
        pivot = np.argmin(arrayM.InBoud[0:arrayM.trackHighInt, 0])
    else:
        pivot = np.argmin(arrayM.InBoud[0:arrayM.trackHighInt, 1])
    '''
    #New version of displacement is done by given x, y coordinate directly
    RotRad = rotation_angle * np.pi / 180
    RotVec = [[np.cos(RotRad), -1*np.sin(RotRad)], [np.sin(RotRad), np.cos(RotRad)]]

    '''
    for Roti in range(arrayM.trackHighInt):
        if Roti == pivot:
            Vec[Roti, :] = 0
        else:
            Vec[Roti, :] = np.subtract(arrayM.InBoud[Roti, :], arrayM.InBoud[pivot, :])
    '''
    ##The calculation of new displacement
    pivot = np.hstack((pivot_x * np.ones((arrayM.trackHighInt, 1)), pivot_y * np.ones((arrayM.trackHighInt, 1))))
    vector = arrayM.InBoud - pivot
    RotCal = np.subtract(np.dot(vector, RotVec), vector)
    InBoudMoveX = np.multiply(np.ones((arrayM.trackHighInt, 1)), tranX)
    InBoudMoveY = np.multiply(np.ones((arrayM.trackHighInt, 1)), tranY)
    OutBoudMov = np.zeros((arrayM.trackHighOut, 2))

    dispMv = np.add(RotCal, np.hstack((InBoudMoveX, InBoudMoveY)))
    dispAl = np.subtract(arrayT.InBoud[0:arrayT.trackHighInt, :], arrayM.InBoud[0:arrayM.trackHighInt, :])
    disp = np.add(dispAl, dispMv)
    dispTot = np.concatenate((disp, OutBoudMov), axis=0)

    return dispTot
######################End of Function of displacement##########################