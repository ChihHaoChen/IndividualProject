# -*- coding: utf-8 -*-
"""
Created on Sat Jul 06 15:16:54 2013

@author: Chih-Hao Chen
"""
#############Import Libraries
import numpy as np
from FunctBank import fun
from MeshQuality import Quality
import Matrix_Cal

##Function Displacement##
##Function of Interpolation##
class interpolation_compute(Quality):
    def __init__(self, arrM, disp, dim, meshIndex, radius, radiusCTPS):
        self.arrM = arrM
        self.disp = disp
        self.dim = dim
        self.meshIndex = meshIndex
        self.radius = radius
        self.radiusCTPS = radiusCTPS

        self.skew_max_ave = 0
        self.non_orth_angle_ave = 0

        self.skew_max = np.zeros([])
        self.nonorth_angle = np.zeros([])
        self.optimise_method = ""
        self.Coor = self.arrM.CoorImportArray.astype(float)
        self.CoorNum = self.arrM.CoorImportArray[:, 0:1].astype(int)
        self.interpolation()

    def interpolation_TPS(self):
        CoorTPS = np.zeros((self.Coor.shape[0], self.arrM.BoudBef.shape[0]))
        meshTPS = np.zeros((self.arrM.BoudBef.shape[0], self.arrM.BoudBef.shape[0]))
        xcoor = self.arrM.BoudBef[:, 0]
        ycoor = self.arrM.BoudBef[:, 1]
        xcoorT = np.transpose(xcoor)
        ycoorT = np.transpose(ycoor)
        distSqX = np.square(xcoor[:, np.newaxis]-xcoorT)
        distSqY = np.square(ycoor[:, np.newaxis]-ycoorT)
        distSq = distSqX+distSqY

        ##TPS index
        nonzerobou = distSq != 0

        Funct_Int = fun(distSq[nonzerobou], self.radius)
        meshTPS[nonzerobou] = Funct_Int.TPS()

        CoorAllX = self.Coor[:, 1]
        CoorAllY = self.Coor[:, 2]
        distNdX = np.square(CoorAllX[:, np.newaxis]-xcoorT)
        distNdY = np.square(CoorAllY[:, np.newaxis]-ycoorT)
        distNd = distNdX+distNdY

        nonzeroindx = distNd != 0
        FunctOut = fun(distNd[nonzeroindx], self.radius)
        CoorTPS[nonzeroindx] = FunctOut.TPS()
        meshCal_TPS = Matrix_Cal.Cal(self.dim, self.disp, self.Coor, self.arrM, meshTPS, CoorTPS)
        self.meshManipTPS = np.add(self.Coor[:, 1: self.dim+1], meshCal_TPS)
        self.MeshQuTPS = Quality(self.arrM, self.meshManipTPS)

    def interpolation_CPC2(self):
        meshCPCsq = np.zeros((self.arrM.BoudBef.shape[0], self.arrM.BoudBef.shape[0]))
        xcoor = self.arrM.BoudBef[:, 0]
        ycoor = self.arrM.BoudBef[:, 1]
        xcoorT = np.transpose(xcoor)
        ycoorT = np.transpose(ycoor)
        distSqX = np.square(xcoor[:, np.newaxis]-xcoorT)
        distSqY = np.square(ycoor[:, np.newaxis]-ycoorT)
        distSq = distSqX+distSqY
        distSqCompact = np.sqrt(distSq)/self.radius

        ##Local index
        localindex = distSqCompact <= 1
        Funct_IntCompact = fun(distSq[localindex], self.radius)
        meshCPCsq[localindex] = Funct_IntCompact.CPCsq()
        CoorCPCsq = np.zeros((self.Coor.shape[0], self.arrM.BoudBef.shape[0]))

        CoorAllX = self.Coor[:, 1]
        CoorAllY = self.Coor[:, 2]
        distNdX = np.square(CoorAllX[:, np.newaxis]-xcoorT)
        distNdY = np.square(CoorAllY[:, np.newaxis]-ycoorT)
        distNd = distNdX+distNdY
        distNdCompact = np.sqrt(distNd)/self.radius
        localout = distNdCompact <= 1

        FunctOutCompact = fun(distNd[localout], self.radius)
        CoorCPCsq[localout] = FunctOutCompact.CPCsq()
        meshCal_CPCsq = Matrix_Cal.Cal(self.dim, self.disp, self.Coor, self.arrM, meshCPCsq, CoorCPCsq)
        self.meshManipCPCsq = np.add(self.Coor[:, 1: self.dim+1], meshCal_CPCsq)
        self.MeshQuCPC2 = Quality(self.arrM, self.meshManipCPCsq)

    def interpolation_CTPSC2B(self):
        meshCTPSC2B = np.ones((self.arrM.BoudBef.shape[0], self.arrM.BoudBef.shape[0]))
        xcoor = self.arrM.BoudBef[:, 0]
        ycoor = self.arrM.BoudBef[:, 1]
        xcoorT = np.transpose(xcoor)
        ycoorT = np.transpose(ycoor)
        distSqX = np.square(xcoor[:, np.newaxis]-xcoorT)
        distSqY = np.square(ycoor[:, np.newaxis]-ycoorT)
        distSq = distSqX+distSqY
        distSqCompactCTPS = np.sqrt(distSq)/self.radiusCTPS
        localindexCTPS = distSqCompactCTPS <= 1

        ##Local index for CTPS Series
        nonzerobou = distSq != 0
        localFinal = np.logical_and(nonzerobou, localindexCTPS)
        ##For RBF with local support only CTPS
        Funct_IntCompactCTPS = fun(distSq[localFinal], self.radiusCTPS)
        meshCTPSC2B[localFinal] = Funct_IntCompactCTPS.CTPSC2B()

        CoorCTPSC2B = np.ones((self.Coor.shape[0], self.arrM.BoudBef.shape[0]))
        CoorAllX = self.Coor[:, 1]
        CoorAllY = self.Coor[:, 2]
        distNdX = np.square(CoorAllX[:, np.newaxis]-xcoorT)
        distNdY = np.square(CoorAllY[:, np.newaxis]-ycoorT)
        distNd = distNdX+distNdY
        nonzeroindx = distNd != 0
        distNdCompactCTPS = np.sqrt(distNd)/self.radiusCTPS
        localoutCTPS = distNdCompactCTPS <= 1
        localoutFinal = np.logical_and(nonzeroindx, localoutCTPS)
        ##for RBF with local support only for CTPS
        FunctOutCompactCTPS = fun(distNd[localoutFinal], self.radius)
        CoorCTPSC2B[localoutFinal] = FunctOutCompactCTPS.CTPSC2B()
        meshCal_CTPSC2B = Matrix_Cal.Cal(self.dim, self.disp, self.Coor, self.arrM, meshCTPSC2B, CoorCTPSC2B)
        self.meshManipCTPSC2B = np.add(self.Coor[:, 1: self.dim+1], meshCal_CTPSC2B)
        self.MeshQuCTPSC2B = Quality(self.arrM, self.meshManipCTPSC2B)

    def interpolation_CPC4(self):
        meshCPC4 = np.zeros((self.arrM.BoudBef.shape[0], self.arrM.BoudBef.shape[0]))
        #To calculate the interpolants for mesh generations, and apply the obtained
        #interpolatns to all the nodes in the meshes, and calculate the new coordinates
        xcoor = self.arrM.BoudBef[:, 0]
        ycoor = self.arrM.BoudBef[:, 1]
        xcoorT = np.transpose(xcoor)
        ycoorT = np.transpose(ycoor)
        distSqX = np.square(xcoor[:, np.newaxis]-xcoorT)
        distSqY = np.square(ycoor[:, np.newaxis]-ycoorT)
        distSq = distSqX+distSqY
        distSqCompact = np.sqrt(distSq)/self.radius

        ##Local index
        localindex = distSqCompact <= 1
        Funct_IntCompact = fun(distSq[localindex], self.radius)

        meshCPC4[localindex] = Funct_IntCompact.CPC4()
        CoorCPC4 = np.zeros((self.Coor.shape[0], self.arrM.BoudBef.shape[0]))

        CoorAllX = self.Coor[:, 1]
        CoorAllY = self.Coor[:, 2]
        distNdX = np.square(CoorAllX[:, np.newaxis]-xcoorT)
        distNdY = np.square(CoorAllY[:, np.newaxis]-ycoorT)
        distNd = distNdX+distNdY
        distNdCompact = np.sqrt(distNd)/self.radius
        localout = distNdCompact <= 1
        FunctOutCompact = fun(distNd[localout], self.radius)
        CoorCPC4[localout] = FunctOutCompact.CPC4()
        meshCal_CPC4 = Matrix_Cal.Cal(self.dim, self.disp, self.Coor, self.arrM, meshCPC4, CoorCPC4)
        self.meshManipCPC4 = np.add(self.Coor[:, 1: self.dim+1], meshCal_CPC4)
        self.MeshQuCPC4 = Quality(self.arrM, self.meshManipCPC4)

    def interpolation_CPC6(self):

        meshCPC6 = np.zeros((self.arrM.BoudBef.shape[0], self.arrM.BoudBef.shape[0]))
        #To calculate the interpolants for mesh generations, and apply the obtained
        #interpolatns to all the nodes in the meshes, and calculate the new coordinates
        xcoor = self.arrM.BoudBef[:, 0]
        ycoor = self.arrM.BoudBef[:, 1]
        xcoorT = np.transpose(xcoor)
        ycoorT = np.transpose(ycoor)
        distSqX = np.square(xcoor[:, np.newaxis]-xcoorT)
        distSqY = np.square(ycoor[:, np.newaxis]-ycoorT)
        distSq = distSqX+distSqY
        distSqCompact = np.sqrt(distSq)/self.radius
        nonzerobou = distSq != 0
        ##Local index
        localindex = distSqCompact <= 1
        Funct_IntCompact = fun(distSq[localindex], self.radius)
        meshCPC6[localindex] = Funct_IntCompact.CPC6()

        CoorCPC6 = np.zeros((self.Coor.shape[0], self.arrM.BoudBef.shape[0]))

        CoorAllX = self.Coor[:, 1]
        CoorAllY = self.Coor[:, 2]
        distNdX = np.square(CoorAllX[:, np.newaxis]-xcoorT)
        distNdY = np.square(CoorAllY[:, np.newaxis]-ycoorT)
        distNd = distNdX+distNdY
        distNdCompact = np.sqrt(distNd)/self.radius
        #To find nonzero terms, then using Broadcast to modify all the elements
        #By broadcast function, the for loops are replaced
        nonzeroindx = distNd != 0
        localout = distNdCompact <= 1
        ##for RBF with local suppor except CTPS
        FunctOutCompact = fun(distNd[localout], self.radius)
        CoorCPC6[localout] = FunctOutCompact.CPC6()
        meshCal_CPC6 = Matrix_Cal.Cal(self.dim, self.disp, self.Coor, self.arrM, meshCPC6, CoorCPC6)
        self.meshManipCPC6 = np.add(self.Coor[:, 1: self.dim+1], meshCal_CPC6)
        self.MeshQuCPC6 = Quality(self.arrM, self.meshManipCPC6)


    def interpolation(self):
        if self.meshIndex == 1 or self.meshIndex == 2:
            self.interpolation_TPS()
            self.interpolation_CPC2()
            self.interpolation_CTPSC2B()
            self.interpolation_CPC4()
            self.interpolation_CPC6()
            MeshResultsArray = np.hstack((self.MeshQuTPS.MeshQuality(), self.MeshQuCPC2.MeshQuality(), \
                                      self.MeshQuCTPSC2B.MeshQuality(), self.MeshQuCPC4.MeshQuality(), \
                                      self.MeshQuCPC6.MeshQuality()))

            Index = np.argmin(MeshResultsArray, axis=1)
        if self.meshIndex == 1:
            IndexPickup = Index[0]
            if IndexPickup == 0:
                mesh_final_out = np.concatenate((self.meshManipTPS, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuTPS.skewMax
                self.nonorth_angle = self.MeshQuTPS.NonOrthAngle
                self.optimise_method = " by TPS"

            elif IndexPickup == 1:
                mesh_final_out = np.concatenate((self.meshManipCPCsq, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCPC2.skewMax
                self.nonorth_angle = self.MeshQuCPC2.NonOrthAngle
                self.optimise_method = " by CPC2"

            elif IndexPickup == 2:
                mesh_final_out = np.concatenate((self.meshManipCTPSC2B, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCTPSC2B.skewMax
                self.nonorth_angle = self.MeshQuCTPSC2B.NonOrthAngle
                self.optimise_method = " by CTPSC2b"

            elif IndexPickup == 3:
                mesh_final_out = np.concatenate((self.meshManipCPC4, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCPC4.skewMax
                self.nonorth_angle = self.MeshQuCPC4.NonOrthAngle
                self.optimise_method = " by CPC4"

            elif IndexPickup == 4:
                mesh_final_out = np.concatenate((self.meshManipCPC6, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCPC6.skewMax
                self.nonorth_angle = self.MeshQuCPC6.NonOrthAngle
                self.optimise_method = " by CPC6"

            else:
                pass

            self.skew_max_ave = MeshResultsArray[0, IndexPickup]
            self.non_orth_angle_ave = MeshResultsArray[1, IndexPickup]

        elif self.meshIndex == 2:
            IndexPickup = Index[1]
            if IndexPickup == 0:
                mesh_final_out = np.concatenate((self.meshManipTPS, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuTPS.skewMax
                self.nonorth_angle = self.MeshQuTPS.NonOrthAngle
                self.optimise_method = " by TPS"

            elif IndexPickup == 1:
                mesh_final_out = np.concatenate((self.meshManipCPCsq, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCPC2.skewMax
                self.nonorth_angle = self.MeshQuCPC2.NonOrthAngle
                self.optimise_method = " by CPC2"

            elif IndexPickup == 2:
                mesh_final_out = np.concatenate((self.meshManipCTPSC2B, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCTPSC2B.skewMax
                self.nonorth_angle = self.MeshQuCTPSC2B.NonOrthAngle
                self.optimise_method = " by CTPSC2b"

            elif IndexPickup == 3:
                mesh_final_out = np.concatenate((self.meshManipCPC4, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCPC4.skewMax
                self.nonorth_angle = self.MeshQuCPC4.NonOrthAngle
                self.optimise_method = " by CPC4"

            elif IndexPickup == 4:
                mesh_final_out = np.concatenate((self.meshManipCPC6, self.Coor[:, 3:4]), axis=1)
                self.skew_max = self.MeshQuCPC6.skewMax
                self.nonorth_angle = self.MeshQuCPC6.NonOrthAngle
                self.optimise_method = " by CPC6"

            else:
                pass

            self.skew_max_ave = MeshResultsArray[0, IndexPickup]
            self.non_orth_angle_ave = MeshResultsArray[1, IndexPickup]

        elif self.meshIndex == 3:
            self.interpolation_TPS()
            mesh_final_out = np.concatenate((self.meshManipTPS, self.Coor[:, 3:4]), axis=1)
            self.skew_max_ave = self.MeshQuTPS.MeshQuality()[0, 0]
            self.non_orth_angle_ave = self.MeshQuTPS.MeshQuality()[1, 0]
            self.skew_max = self.MeshQuTPS.skewMax
            self.nonorth_angle = self.MeshQuTPS.NonOrthAngle
            self.optimise_method = " by TPS"

        elif self.meshIndex == 4:
            self.interpolation_CPC2()
            mesh_final_out = np.concatenate((self.meshManipCPCsq, self.Coor[:, 3:4]), axis=1)
            self.skew_max_ave = self.MeshQuCPC2.MeshQuality()[0, 0]
            self.non_orth_angle_ave = self.MeshQuCPC2.MeshQuality()[1, 0]
            self.skew_max = self.MeshQuCPC2.skewMax
            self.nonorth_angle = self.MeshQuCPC2.NonOrthAngle
            self.optimise_method = " by CPC2"

        elif self.meshIndex == 5:
            self.interpolation_CTPSC2B()
            mesh_final_out = np.concatenate((self.meshManipCTPSC2B, self.Coor[:, 3:4]), axis=1)
            self.skew_max_ave = self.MeshQuCTPSC2B.MeshQuality()[0, 0]
            self.non_orth_angle_ave = self.MeshQuCTPSC2B.MeshQuality()[1, 0]
            self.skew_max = self.MeshQuCTPSC2B.skewMax
            self.nonorth_angle = self.MeshQuCTPSC2B.NonOrthAngle
            self.optimise_method = " by CTPSC2b"

        elif self.meshIndex == 6:
            self.interpolation_CPC4()
            mesh_final_out = np.concatenate((self.meshManipCPC4, self.Coor[:, 3:4]), axis=1)
            self.skew_max_ave = self.MeshQuCPC4.MeshQuality()[0, 0]
            self.non_orth_angle_ave = self.MeshQuCPC4.MeshQuality()[1, 0]
            self.skew_max = self.MeshQuCPC4.skewMax
            self.nonorth_angle = self.MeshQuCPC4.NonOrthAngle
            self.optimise_method = " by CPC4"

        elif self.meshIndex == 7:
            self.interpolation_CPC6()
            mesh_final_out = np.concatenate((self.meshManipCPC6, self.Coor[:, 3:4]), axis=1)
            self.skew_max_ave = self.MeshQuCPC6.MeshQuality()[0, 0]
            self.non_orth_angle_ave = self.MeshQuCPC6.MeshQuality()[1, 0]
            self.skew_max = self.MeshQuCPC6.skewMax
            self.nonorth_angle = self.MeshQuCPC6.NonOrthAngle
            self.optimise_method = " by CPC6"

        else:
            pass

        meshFinal2 = mesh_final_out.tolist()
        CoorNumList = self.CoorNum.tolist()

        Container=[]
        for Namei in range(len(CoorNumList)):
            Container.append(CoorNumList[Namei])
            Container[Namei].append(meshFinal2[Namei])

        meshFinal3 = "\n".join(repr(x) for x in Container)

        meshFinal = meshFinal3.replace(',', ' ').replace('[', '').replace(']', '')
        return meshFinal




