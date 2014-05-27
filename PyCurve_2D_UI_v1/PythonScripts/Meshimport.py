# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 15:19:15 2013
This for mesh generation
@author: Chih-Hao Chen
"""

#############Import Libraries
import numpy as np
###########################################################
###############################################################################
#Lists of boundaries
class mesh:
    def __init__(self,dim,order,phyLine,phyLineOut):
        self.dim=dim
        self.order=order
        self.phyLine=phyLine
        self.phyLineOut=phyLineOut
        
    def openfile(self,x):
        self.fid=open(x, "r")
        self.arrayString=[]  
        for line in self.fid:
            self.arrayString.append(line)
                        
        self.fid.close()
        
        #Use str.index to find the index of certain strings
        self.CoorFlagStart=self.arrayString.index('$Nodes\n''')+2
        self.CoorFlagEnd=self.arrayString.index('$EndNodes\n''')-1
        self.EleFlagStart=self.arrayString.index('$Elements\n''')+2
        self.EleFlagEnd=self.arrayString.index('$EndElements\n''')-1
        ###########################################################################
        #Delcaration of an array to store all the information of coordinates
        self.CoorImport=[[] for i in range(self.CoorFlagEnd-self.CoorFlagStart+1)]
        for CoorI in range(self.CoorFlagEnd-self.CoorFlagStart+1):
            self.CoorImport[CoorI]=self.arrayString[self.CoorFlagStart+CoorI].strip().split()
         
        self.CoorImportArray=np.array(self.CoorImport)
        ###########################################################################
        #Delcaration of an array to store all the information of coordinates
        EleImport=[[] for j in range(self.EleFlagEnd-self.EleFlagStart+1)]
        for EleI in range(self.EleFlagEnd-self.EleFlagStart+1):
            EleImport[EleI]=self.arrayString[self.EleFlagStart+EleI].strip().split()
           
        self.EleImportArray=np.array(EleImport,ndmin=2)
        nbHighLine=0
        self.nbHighQuads=0
        self.nbNode=0
        nbHighInLine=0
        nbHighOutLine=0
        NodeIn=[]
        NodeOut=[]
        Node=[]
        NodeEle=[]
        NodeCom=[]
        GridName=np.zeros((self.EleFlagEnd-self.EleFlagStart,4))
           
        for iii in range(self.EleImportArray.size):
            EleContainer=np.ravel(self.EleImportArray[0,iii])
            if EleContainer[1] == '8'  or  EleContainer[1] == '26' or\
                EleContainer[1] == '27' or  EleContainer[1] == '28' or\
                EleContainer[1] == '62' or  EleContainer[1] == '63' or\
                EleContainer[1] == '64' or  EleContainer[1] == '65' or\
                EleContainer[1] == '1' :
                nbHighLine+=1
                NodeCom+=np.column_stack(EleContainer[5:6])
                for boui in range(len(self.phyLine)):
                   if EleContainer[4] == self.phyLine[boui]:
                        nbHighInLine+=1
                        NodeIn+=np.column_stack(EleContainer[5:])

                for bouj in range(len(self.phyLineOut)):
                    if EleContainer[4] == self.phyLineOut[bouj]:
                        nbHighOutLine+=1
                        NodeOut+=np.column_stack(EleContainer[5:])                     
                
            #To find out the node of grids for quadrilateral elements
            elif EleContainer[1] == '10' or  EleContainer[1] == '36' or\
                EleContainer[1] == '37' or  EleContainer[1] == '38' or\
                EleContainer[1] == '47' or  EleContainer[1] == '48' or\
                EleContainer[1] == '49' or  EleContainer[1] == '50' or\
                EleContainer[1] == '3' :  
                GridName[self.nbHighQuads,0:4]=EleContainer[5:9]
                self.nbHighQuads+=1
            #This added part is used to detect nodes
            elif EleContainer[1]=='15':
                Node.append(EleContainer[4])
                NodeEle.append(EleContainer)

        self.NodeName=GridName[0:self.nbHighQuads,:].astype(int)
        ##Deletion of redundant nodes in outer boundaries 
        ##By assinging all impossible values, then compare and replace if
        ##new element         
        NodeOutArray=np.ravel(np.array(NodeOut).astype(int))
        NodeInArray=np.ravel(np.array(NodeIn).astype(int))
        
        self.NodeOutBoud=-1000*np.ones(NodeOutArray.shape[0])
        self.trackHighOut=0
        for OutIndxi in range(NodeOutArray.shape[0]):
            NodeOK=0
            for OutIndxj in range(self.trackHighOut):
                if self.NodeOutBoud[OutIndxj]==NodeOutArray[OutIndxi]:
                    NodeOK+=1
            if NodeOK==0:
                self.NodeOutBoud[self.trackHighOut]=NodeOutArray[OutIndxi]
                self.trackHighOut+=1
         
        self.NodeInBoud=-1000*np.ones(NodeInArray.shape[0])
        self.trackHighInt=0
        
        for InIndxi in range(NodeInArray.shape[0]):
            NodeOK=0
            for InIndxj in range(self.trackHighInt):
                if self.NodeInBoud[InIndxj]==NodeInArray[InIndxi]:
                    NodeOK+=1
            if NodeOK==0:
                self.NodeInBoud[self.trackHighInt]=NodeInArray[InIndxi]
                self.trackHighInt+=1
        
        self.OutBoud=np.column_stack((\
        np.array(self.CoorImportArray[self.NodeOutBoud[0:self.trackHighOut].\
        astype(int)-1,1]).astype(float),\
        np.array(self.CoorImportArray[self.NodeOutBoud[0:self.trackHighOut].\
        astype(int)-1,2]).astype(float)))
        self.InBoud=np.column_stack((\
        np.array(self.CoorImportArray[self.NodeInBoud[0:self.trackHighInt].\
        astype(int)-1,1]).astype(float),\
        np.array(self.CoorImportArray[self.NodeInBoud[0:self.trackHighInt].\
        astype(int)-1,2]).astype(float)))
        self.BoudBef=np.concatenate((self.InBoud,self.OutBoud))
        
        NodeCom1=np.ravel(np.array(NodeCom).astype(int))               
        NodeArray=np.array(Node).astype(int)
        Coincid=NodeArray[:,np.newaxis]-NodeCom1
        NewNodeEle=[]
        CoincidIndx=np.argwhere(Coincid==0)
            
        #For processing of redundant nodes
        for Coini in range(CoincidIndx.shape[0]):
            NewNodeEle.append(NodeEle[CoincidIndx[Coini][0]])
        
        NewNodeEleArr=np.array(NewNodeEle).astype(int)
        self.NewNodePrint="\n".join(repr(x) for x in NewNodeEleArr)
        self.ElementNumb=str(int(self.arrayString[self.EleFlagStart-1])\
                        -len(Node)+CoincidIndx.shape[0])
        self.ElePrintStart=self.EleFlagStart+len(Node)
        self.NewNodePrintStr=self.NewNodePrint.replace('array([','').replace(', ', ' ').replace('])','')

        
        