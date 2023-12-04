class readFile:
    
    def __init__(self,fileName):
        
        # this refers to x - y values
        self.xValues=[]
        self.yValues=[]
        
        f=open(fileName,'r')

        for row in f: 
            if "#" in row:
                pass
            else:
                row = row.split()
                self.xValues.append(float(row[0]))
                self.yValues.append(float(row[1]))
        
    def x(self):
        return self.xValues

    def y(self):
        return self.yValues
import time
class readAutoWrite:
    
    def __init__(self,fileName):
        
        # this refers to x - y values
        self.xValues=[]
        self.totalCellsNumber=1
        
        f=open(fileName,'r')
        counter=0
        i=0

        for row in f:
            counter+=1
            if counter==21:
                row = row.split()
                self.totalCellsNumber=int(row[0])
            elif counter>22:
                i+=1
                row = row.split()
                self.xValues.append(float(row[0]))
            if i==self.totalCellsNumber:
                break
        
    def x(self):
        return self.xValues

import numpy as np

class readVelocity:
    
    def __init__(self,fileName):
        
        # this refers to x - y values
        self.xValues=[]
        self.UxValues=[]
        self.UyValues=[]
        self.UzValues=[]
        
        f=open(fileName,'r')

        for row in f: 
            if "#" in row:
                pass
            else:
                row = row.split()
                c=np.size(row)
                self.xValues.append(float(row[0]))
                
                for i in range(1,c):
                    if "(" or ")" in row:
                        row[i]=row[i].replace("(","")
                        row[i]=row[i].replace(")","")
                    if i==1:
                        self.UxValues.append(float(row[i]))
                    elif i==2:
                        self.UyValues.append(float(row[i]))
                    elif i==3:
                        self.UzValues.append(float(row[i]))
        
    def x(self):
        return self.xValues

    def y(self):
        Umag=[]
        for i in range(0,np.size(readVelocity.x(self))):
            Umag.append((self.UxValues[i]**2+self.UyValues[i]**2+self.UzValues[i]**2)**0.5)
        return Umag
    

class readConvergence:
    
    def __init__(self,fileName):
        
        # this refers to x - y values
        self.xValues=[]
        self.kValues=[]
        self.epsilonValues=[]
        self.M3Values=[]

        f=open(fileName,'r')

        for row in f: 
            if "#" in row:
                pass
            else:
                row = row.split()
                self.xValues.append(float(row[0]))
                self.kValues.append(float(row[4]))
                self.epsilonValues.append(float(row[5]))
                self.M3Values.append(float(row[6]))                
        
    def x(self):
        return self.xValues

    def k(self):
        return self.kValues

    def epsilon(self):
        return self.epsilonValues

    def M3(self):
        return self.M3Values