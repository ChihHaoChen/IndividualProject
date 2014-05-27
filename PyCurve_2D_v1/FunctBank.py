# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:09:15 2013

@author: Chih-Hao Chen
"""

import numpy as np
import matplotlib.pyplot as plt

class fun:
    def __init__(self, distSq, Radius):
        self.distSq = distSq
        self.Radius = Radius
        self.Eps = np.sqrt(self.distSq)/Radius
        self.EpsSq = np.square(self.Eps)

    def CPC0(self):
        return np.power(np.subtract(1, self.Eps), 2)
     
    def CPCsq(self):
        return np.multiply(np.power(np.subtract(1, self.Eps), 4), np.add(1, np.multiply(4, self.Eps)))
        
    def CPC4(self):
        return np.multiply(np.power(np.subtract(1, self.Eps), 6),\
            np.add((35/3)*self.EpsSq, np.add(1, np.multiply(6, self.Eps))))
            
    def CPC6(self):
        return np.multiply(np.power(np.subtract(1, self.Eps), 8),\
            np.add(32*np.multiply(self.Eps, self.EpsSq),\
            np.add(25*self.EpsSq, np.add(1, np.multiply(8, self.Eps)))))
            
    def CTPSC0(self):
        return np.power(np.subtract(1, self.Eps), 5)
        
    def CTPSC1(self):
        return 1+80*self.EpsSq/3-40*np.power(self.Eps, 3)+15*np.power(self.EpsSq, 2)-\
                (8/3)*np.power(self.Eps, 5)+20*np.multiply(self.EpsSq, np.log10(self.Eps))
    
    def CTPSC2A(self):
        return 1-30*self.EpsSq-10*np.power(self.Eps, 3)+45*np.power(self.EpsSq, 2)-\
                6*np.power(self.Eps, 5)-60*np.multiply(np.power(self.Eps, 3), np.log10(self.Eps))
    
    def CTPSC2B(self):
        return 1-20*self.EpsSq+80*np.power(self.Eps, 3)-45*np.power(self.EpsSq, 2)-\
                16*np.power(self.Eps, 5)+60*np.multiply(np.power(self.EpsSq, 2), np.log10(self.Eps))
                
    def TPS(self):
        return np.multiply(self.distSq, np.log10(np.sqrt(self.distSq)))
        
    def MQB(self):
        return np.sqrt(np.add(np.square(0.001), self.distSq))
        
    def IMQB(self):
        return np.divide(1, np.sqrt(np.add(np.square(0.001), self.distSq)))
        
    def QB(self):
        return np.add(1, self.distSq)
        
    def IQB(self):
        return np.divide(1, np.add(1, self.distSq))
        
    def Gaussian(self):
        return np.exp(np.multiply(-1, self.distSq))


class qplot:
    def __init__(self, parent=None):
        super(self).__init__(parent)

    def quality_plot(self, skewMaxAve, NonOrthAngleAve, skewMax, NonOrthAngle):
        meanStrSkew = 'Skewness with average of '
        meanStrSkew += skewMaxAve.astype(str)
        meanStrNonOr = 'Non-Orthognality with average of '
        meanStrNonOr += NonOrthAngleAve.astype(str)

        self.skew_dist = skewMax
        self.non_orth_angle_dist = NonOrthAngle
        # Two subplots, the axes array is 1-d
        f, axarr = plt.subplots(2) 
        
        axarr[0].hist(skewMax, bins=100, range=(0, 1))
        axarr[0].set_title(meanStrSkew)
        axarr[0].grid()
        axarr[1].hist(NonOrthAngle, bins=100, range=(0, 90))
        axarr[1].set_title(meanStrNonOr)
        axarr[1].grid()

        plt.show()


          


