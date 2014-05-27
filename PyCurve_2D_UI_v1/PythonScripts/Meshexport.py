# -*- coding: utf-8 -*-
"""
Created on Sat Jul 06 15:47:35 2013

@author: Chih-Hao Chen
"""

##Function of export
def meshexport(meshM,meshFinal,savefile):
    
    genid = open(savefile, 'w')
    genid.writelines(meshM.arrayString[0:meshM.CoorFlagStart])
    genid.writelines(meshFinal)
    genid.writelines('\n')
    genid.writelines(meshM.arrayString[meshM.EleFlagStart-3:meshM.EleFlagStart-1])
    genid.writelines(meshM.ElementNumb)
    genid.writelines('\n')
    genid.writelines(meshM.NewNodePrintStr)
    genid.writelines('\n')
    genid.writelines(meshM.arrayString[meshM.ElePrintStart:len(meshM.arrayString)])
    genid.close()
