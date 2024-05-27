import numpy as np
import os,sys
from .parameters import Parameters

class Shifter:

    '''
    Shift the original POSCAR to desired stackings.
    '''

    inDir = ""  
    outDirSTR = ""      #Directory to store shifted structures.
    latVec = []
    nSecs = 0           #number of sectors to shift
    splitter = None   #split function split(data, i), return the bool value deciding if the atom is on top or on bottom
    delta = []
    Header = ''         #Header of POSCAR
    minHeader = ''
    data = []
    splitDat = []
    elements = []       #element contained in POSCAR
    nElement = {        #Number of atoms of each element
    }
    mode = {
        "Cartesian":0,
        "Direct":1,
    }
    def update(self, nSecs=0, inDir:str=None, splitter=None):

        '''
        Update parameters of shifter.
        '''

        if splitter is not None:
            self.splitter = splitter
        else:
            self.splitter = self.splz
        if inDir is not None:
            self.inDir = inDir
        if nSecs > 0 :
            self.nSecs = nSecs

    def __init__(self, param:Parameters, splitter=None):
        self.update(param.nSecs, param.originDir, splitter)
        self.readData(self.inDir)

    def readData(self, inDir):

        '''
        Read POSCAR data from Infile.
        '''

        infile = open(inDir + "/POSCAR_top", 'r')
        for i,line in enumerate(infile):
            words = line.split()
            if i<8 :
                if i in range(2,5):
                    #read lattice vectors.
                    self.latVec.append([float(words[0]),float(words[1]),float(words[2])])
                    self.delta.append([float(words[0])/self.nSecs,float(words[1])/self.nSecs])
                if i == 5:
                    #read elements
                    self.elements = [element for element in words]
                if i == 6:
                    for i,word in enumerate(words):
                        self.nElement[self.elements[i]] = int(word)
                if i == 7:
                    #read mode, Cartesian or Direct
                    readMode = self.mode[words[0]]
                    self.Header += "Cartesian\n"
                    self.minHeader += "Selective dynamics\nCartesian\n"
                    continue
                self.Header += line
                self.minHeader += line
            else:
                #read the main part of POSCAR
                if readMode == 0:
                    x = float(words[0])
                    y = float(words[1])
                    z = float(words[2])
                else:
                    a = float(words[0])
                    b = float(words[1])
                    c = float(words[2])
                    x = a*self.latVec[0][0] + b*self.latVec[1][0] + c*self.latVec[2][0]
                    y = a*self.latVec[0][1] + b*self.latVec[1][1] + c*self.latVec[2][1]
                    z = a*self.latVec[0][2] + b*self.latVec[1][2] + c*self.latVec[2][2]
                self.data.append([x,y,z])
        self.splitDat = self.splitter(data=self.data)
        infile.close()

    def shift(self,i,j,outDir):

        '''
        shift atoms at top layer by (i/nSecs*latVecA,i/nSecs*latVecB)
        '''

        writePath = outDir+'/'+str(i)+'_'+str(j)
        if not os.path.exists(writePath):
            os.makedirs(writePath)             #create write path if not exist
        out = open(writePath+"/POSCAR", "w")

        out.write(self.Header)              #write the header of POSCAR
        for idx,line in enumerate(self.data):
            delta = []
            if self.splitDat[idx]:  #top layer atoms
                delta = [self.delta[0][0] * i + self.delta[1][0] * j, self.delta[0][1] * i + self.delta[1][1] * j]
            else:
                delta = [0, 0]
            for dim in range(2):
                out.write(str(line[dim] + delta[dim]))
                out.write(' ')
            out.write(str(line[2]))
            out.write('\n')
        out.close()

    def shift4Min(self, i, j, outDir):
        '''
        shift atoms at top layer by (i/nSecs*latVecA,i/nSecs*latVecB)
        '''
        flag = [True, True]
        writePath = outDir+'/'+str(i)+'_'+str(j)
        if not os.path.exists(writePath):
            os.makedirs(writePath)             #create write path if not exist
        out = open(writePath+"/POSCAR", "w")
        out.write(self.minHeader)              #write the header of POSCAR
        for idx,line in enumerate(self.data):
            delta = []
            if self.splitDat[idx]:  #top layer atoms
                delta = [self.delta[0][0] * i + self.delta[1][0] * j, self.delta[0][1] * i + self.delta[1][1] * j, 0]
            else:
                delta = [0, 0, 0]
            for dim in range(3):
                out.write(str(line[dim] + delta[dim]))
                out.write(' ')
            if self.splitDat[idx] and flag[0]:
                out.write('F F T \n')
                flag[0] = False
            elif (not self.splitDat[idx]) and flag[1]:
                out.write('F F T \n')
                flag[1] = False
            else:
                out.write('T T T \n')
        out.close()
    
    def shiftAll(self, outDir):
        for i in range(self.nSecs):
            for j in range(self.nSecs):
                self.shift(i, j, outDir=outDir)

    def splz(self,data):

        '''
        an example of spliter.
        '''

        splitDat = []
        avgZ = 0
        for i in range(len(data)):
            avgZ += data[i][2]
        avgZ /= len(data)
        for i in range(len(data)):
            if data[i][2]>=avgZ:
                splitDat.append(True)
            else:
                splitDat.append(False)
        return splitDat
        
