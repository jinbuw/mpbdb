from math import sqrt
import numpy as Npy
_array = Npy.array
    
class Vector:
    def cross(self, v1, v2):
        "Returns the cross product of two vectors """
        return _array([(v1[1]*v2[2] - v1[2]*v2[1]),
                      (v1[2]*v2[0] - v1[0]*v2[2]),
                      (v1[0]*v2[1] - v1[1]*v2[0])])
    
    
    def prnVec(self, vec, cmt= " "):
        print("# vec %s : %s " % (cmt, self.getVecStr(vec)))

    
    def getVecLen(self, vec):
        return  sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])
    
    
    def getVectNorm(self, vec):
        return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2] * vec[2])
    
    def getVecInnerProduct(self, v1, v2):
        return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])
    
    def getDistSquare(self, ca, cb):
        dx = ca[0] - cb[0]
        dy = ca[1] - cb[1]
        dz = ca[2] - cb[2]
        dr2 = dx * dx + dy*dy + dz*dz 
        return dr2 
    
    
    def getDist(self, ca, cb):
        dx = ca[0] - cb[0]
        dy = ca[1] - cb[1]
        dz = ca[2] - cb[2]
        dr2 = dx * dx + dy*dy + dz*dz 
        dr = sqrt(dr2)
        return dr 
 
     