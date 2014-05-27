# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:07:53 2013

@author: Chih-Hao Chen
"""
import numpy as np



class Quality:
    def __init__(self, meshArr, meshManip):
        self.meshArr = meshArr
        self.meshManip = meshManip
        self.skewCal = 0
        self.NonOrthAngle = 0
        self.mesh_quality_re = np.zeros((1, 2))
        self.skewMax = 0
        self.skewMaxAve = 0
        self.NonOrthAngleAve = 0


    def Skewness(self):
        ##This section is for calculations of mesh qualities
        #Re-map the elements grid to the new coordinates for update and for manipulation of mesh qualities
        self.GridPosX=self.meshManip[self.meshArr.NodeName-1,0]
        self.GridPosY=self.meshManip[self.meshArr.NodeName-1,1]
    
        GridBD=np.add(np.square(np.subtract(self.GridPosX[:,1],self.GridPosX[:,3])),\
                    np.square(np.subtract(self.GridPosY[:,1],self.GridPosY[:,3])))
        GridAC=np.add(np.square(np.subtract(self.GridPosX[:,0],self.GridPosX[:,2])),\
                    np.square(np.subtract(self.GridPosY[:,0],self.GridPosY[:,2])))
        self.DiffPosX=np.subtract(self.GridPosX,np.roll(self.GridPosX,-1,axis=1))
        self.DiffPosY=np.subtract(self.GridPosY,np.roll(self.GridPosY,-1,axis=1))
        self.DiffSqPosX=np.square(self.DiffPosX)
        self.DiffSqPosY=np.square(self.DiffPosY)
    
        DiffPosXposY=np.add(self.DiffSqPosX,self.DiffSqPosY)
        Angle1=np.divide(np.subtract(np.add(DiffPosXposY[:,3],DiffPosXposY[:,0]),GridBD)\
                , 2*(np.multiply(np.sqrt(DiffPosXposY[:,3]),np.sqrt(DiffPosXposY[:,0]))))
        Angle2=np.divide(np.subtract(np.add(DiffPosXposY[:,0],DiffPosXposY[:,1]),GridAC)\
                , 2*(np.multiply(np.sqrt(DiffPosXposY[:,0]),np.sqrt(DiffPosXposY[:,1]))))    	
        Angle3=np.divide(np.subtract(np.add(DiffPosXposY[:,1],DiffPosXposY[:,2]),GridBD)\
                , 2*(np.multiply(np.sqrt(DiffPosXposY[:,1]),np.sqrt(DiffPosXposY[:,2]))))
        Angle4=np.divide(np.subtract(np.add(DiffPosXposY[:,2],DiffPosXposY[:,3]),GridAC)\
                , 2*(np.multiply(np.sqrt(DiffPosXposY[:,2]),np.sqrt(DiffPosXposY[:,3])))) 
        Angle=np.transpose(np.vstack((np.arccos(Angle1),np.arccos(Angle2),np.arccos(Angle3), \
                    np.arccos(Angle4))))
        halfpi=np.pi/2
        
        self.skewCal=np.zeros((self.meshArr.nbHighQuads, 4))
        #Find out the skewness
        self.skewCal=np.abs((Angle-halfpi)/halfpi)
        self.skewMax=self.skewCal.max(axis=1)
        self.skewMaxAve=np.average(self.skewMax)
        return self.skewMax
        
    def NonOrth(self):
        NodeNameStr=self.meshArr.NodeName.astype(str)
        NodeNameBWStr=np.roll(NodeNameStr,-1,axis=1) 
        LineName11=np.core.defchararray.add(np.core.defchararray.add(NodeNameStr, '00'), NodeNameBWStr).astype(float)
        LineName12=np.core.defchararray.add(np.core.defchararray.add(NodeNameBWStr, '00'), NodeNameStr).astype(float)
        searchlist11=np.ravel(LineName11)
        searchlist12=np.ravel(LineName12)
        BingoIndx1=[[],[]]
        BingoIndx2=[[],[]] 
        BingoReg=[[],[]]
        
        #X and Y Coordinate of Centre of Elements
        CentreOfEleX=np.average(self.GridPosX, axis=1)
        CentreOfEleY=np.average(self.GridPosY, axis=1)
        CentreVec=[[],[]]
        FaceVec=[[],[]]
          
        LineCount=0
        for seIndx in range(len(searchlist11)):
            if searchlist11[seIndx]!=0:
                BingoReg=np.where(LineName11==searchlist12[seIndx],True,False)
                if BingoReg.any():
                    BingoIndx1=np.column_stack((BingoIndx1,np.where(LineName11==searchlist11[seIndx])))
                    BingoIndx2=np.column_stack((BingoIndx2,np.where(LineName11==searchlist12[seIndx])))
                    searchlist11[4*BingoIndx2[0,LineCount]+BingoIndx2[1,LineCount]]=0
                    LineCount+=1
                                
        CentreVec=np.column_stack((np.subtract(CentreOfEleX[BingoIndx1[0, :].astype(int)], \
                                CentreOfEleX[BingoIndx2[0, :].astype(int)]), \
                                np.subtract(CentreOfEleY[BingoIndx1[0, :].astype(int)], \
                                CentreOfEleY[BingoIndx2[0, :].astype(int)])))
        FaceVec=np.column_stack((self.DiffPosX[BingoIndx1[0, :].astype(int), BingoIndx1[1, :].astype(int)],\
                                self.DiffPosY[BingoIndx1[0, :].astype(int), BingoIndx1[1, :].astype(int)]))
                                
        CentreVecLen=np.sqrt(np.sum(np.square(CentreVec), axis=1))
        FaceVecLen=np.sqrt(np.sum(np.square(FaceVec), axis=1))
        self.NonOrthAngle=np.abs(np.arccos((CentreVec[:, 0]*FaceVec[:, 0]+CentreVec[:, 1]*FaceVec[:, 1])\
                /(CentreVecLen*FaceVecLen))-(np.pi/2))*180/np.pi
        self.NonOrthAngleAve=np.average(self.NonOrthAngle, axis=0)
        return self.NonOrthAngle
        
    def MeshQuality(self):
        self.Skewness()
        self.NonOrth()
        #self.mesh_quality_re = [[self.skewMaxAve], [self.NonOrthAngleAve]]
        self.mesh_quality_re = np.vstack((self.skewMaxAve, self.NonOrthAngleAve))
        return self.mesh_quality_re


        





    