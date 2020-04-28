from math import * #log, sin, cos, sqrt, pi, fabs, radians,degrees, acos, asin, atan2, pi, exp, ceil, sqrt,fmod, modf 
import numpy as Npy
_array = Npy.array
pi2 = 2 * pi 

def angleFromSineAndCosine(sine, cosine):
    """ return angle in radians """
    sine = min(1., max(-1., sine))
    angle = Npy.arcsin(abs(sine))
    if sine*cosine < 0.:
        angle = -angle
    if cosine < 0.:
        angle = pi + angle
    if angle < -pi:
        angle = pi2 + angle
    elif angle > pi:
        angle = angle - pi2
    return angle

def getVecLen(vec):
    return  sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])

def getVectNorm( vec):
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2] * vec[2])

def getNormVect( vec):
    return vec/getVectNorm(vec)
    # Npy.add.reduce(self.array*other.array)
 
def cross(va, vb):
    "Returns the cross product with vector |other|."
    #print "# va= %s, vb= %s" % (va, vb)
    #print type(va), type(vb)
    return Npy.array([va[1]*vb[2] - va[2]*vb[1], va[2]*vb[0] - va[0]*vb[2],
                va[0]*vb[1] - va[1]*vb[0]])
    
def dihedral_vect(pv1, pv2, pv3, pv4, deg=False):
    """ return dihedral angle from 4 position vectors in radians3"""
    if pv1 is None or pv2 is None or pv3 is None or pv4 is None: return -999
    v1 = pv2 - pv1
    v2 = pv3 - pv2
    v3 = pv3 - pv4
    
    a = cross(v1, v2)
    a = getNormVect(a)
    
    b = cross(v3, v2)
    b = getNormVect(b)
    #print "a= %s b= %s ab= %s " % (a, b, a*b)
    #print Npy.dot(a, b)
    #mcos = Npy.add.reduce(a*b)
    #print "# mcos: ", mcos 
    mcos = Npy.dot(a,b)
    a_b = cross(a,b)
    msin = Npy.dot(a_b, getNormVect(v2))
    #print "msin= ", msin 
    dihd=  angleFromSineAndCosine(msin, mcos)
    #print " psi(rad)= %6.4f  %6.2f(deg)" % (psi, psi*rad2deg)
    if deg:
        return degrees(dihd)
    else:
        return dihd
    




class Vector:
    
    def normVect(self, vect):
        """ return a normalized vector)"""
        return vect/sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2])    
    
    def cross(self, v1, v2):
        "Returns the cross product of two vectors """
        return _array([(v1[1]*v2[2] - v1[2]*v2[1]),
                      (v1[2]*v2[0] - v1[0]*v2[2]),
                      (v1[0]*v2[1] - v1[1]*v2[0])])
    
    def getVecStr(self, vec, cmt=""):
        return  "# %7.2f %7.2f %7.2f " % (vec[0], vec[1], vec[2])
        
    
    def prnVec(self, vec, cmt= " "):
        print("# vec %s : %s " % (cmt, self.getVecStr(vec)))

    
    def getVecLen(self, vec):
        return  sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])
    
    def getVectSub(self, ca1, ca2):
        return (ca2[0]-ca1[0], ca2[1]-ca1[1], ca2[2]-ca1[2])

    def getVectAdd(self, ca1, ca2):
        return (ca2[0]+ca1[0], ca2[1]+ca1[1], ca2[2]+ca1[2])
    
    def getVectMulbyVal(self, vec, val):
        return Npy.array(vec)*val 
    
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

    def getVecAngle(self, vec1, vec2, deg=0, vb=0):
        """ return angle between vec1 and vec2 in radians """
        vec1 = self.normVect(vec1)
        vec2 = self.normVect(vec2)
        ang =  acos(vec1[0] * vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2])
        if vb: print("# ang: %6.1f " % degrees(ang))
        if deg:
            return degrees(ang)
        else:
            return ang 
    
    #deprecated
    def rot_along_X(self, angle, vect):
        vect_axle = [1, 0, 0]
        return self.arb_rot(vect, vect_axle, angle)

    def rot_along_Y(self, angle, vect):
        vect_axle = [0, 1, 0]
        return self.arb_rot(vect, vect_axle, angle)

    def rot_along_Z(self, angle, vect):
        vect_axle = [0, 0, 1]
        return self.arb_rot(vect, vect_axle, angle)
    

    def arb_rot(self, vect, vect_axis, gamma):
        #print "# gamma= %6.1f" % gamma 
        #print "# ve ang: %6.1f "  % degrees(self.getVecAngle(vect, vect_axis))
        gamma = radians(gamma)
        #print "# gamma: ", gamma 
        vect = self.normVect(vect)
        #vect_axis = normVect(vect_axis)
        x,y,z  = vect[0], vect[1],vect[2]
        v0,v1,v2 = vect_axis[0],vect_axis[1],vect_axis[2]
        v1_2 = v1 * v1 
        v2_2 = v2 * v2 
        lyz_2 = v1_2 + v2_2
        cos_gamma = cos(gamma)
        sin_gamma = sin(gamma)
        #Ry = nmatrix([[lyz,    0.,    -vect_axis[0] ],
                        #[0,      1.,         0 ],
                        #[vect_axis[0],  0,  lyz]])
            
        #Rz = nmatrix([[ cos_gamma,     -sin_gamma,         0],
                        #[ sin_gamma,      cos_gamma,         0],
                        #[0,           0,         1   ]])
        if lyz_2 == 0.0:
            #Rgamma = inv(Ry)*Rz*Ry
            rot_vector = _array([x, y*cos_gamma - z * sin_gamma, z * cos_gamma + y * sin_gamma])
        else:
            #Rx = npy.array([[ 1,           0,          0],
                            #[0,       vect_axis[2]/lyz,   -vect_axis[1]/lyz],
                            #[0,       vect_axis[1]/lyz,   vect_axis[2]/lyz ]])
            #Rgamma = inv(Rx)*inv(Ry)*Rz*Ry*Rx;
            v0x = v0 * x 
            v1y = v1 * y
            v2z = v2 * z 
            vv = v0x + v1y + v2z 
            vx = v0*vv -(-lyz_2*x +v0*v1y +v0*v2z)*cos_gamma + \
                 (-v2*y + v1*z) * sin_gamma
            vy = v1 * vv - (v1*v0x + (-1.0 + v1_2)*y + v1*v2z)*cos_gamma + \
                 (v2*x - v0*z) * sin_gamma 
            vz = v2*v0x + v2*v1y + v2_2 * z + \
               (-v2*v0x + (1.0 - v2_2)*z -v2*v1y) * cos_gamma + \
               (-v1*x +v0*y)*sin_gamma 
            rot_vector = _array([vx,vy,vz])
        return  self.normVect(rot_vector)

    
   


    def vector2thePhi(self, axis, vb=1):
        """ transfer vector in Cartesian space to the, phi(in units of degree) in Polar space """
        rAB = self.normVect(axis)
        if vb: self.prnVec(rAB, cmt="# normal axis: ")
        _theta  = acos(rAB[2])
        if fabs(rAB[0]) < 0.001 and fabs(rAB[1]) < 0.001:
            _phi = 0
        else:
            _phi = atan2(rAB[1],rAB[0])
        if vb: print("# bond: theta= %6.1f  phi= %6.1f " % (degrees(_theta), degrees(_phi)))
        return degrees(_theta), degrees(_phi)
    
    
    
def std( set, vb=0, name="", NUM=False):
    """ calculate the standard deviation of one set of data (a list) 
        return av_sum, var_std, v_max, v_min
    """
    num = len(set)
    v_max = max(set)
    v_min = min(set)
    v_rang = v_max - v_min
    _sum  = 0
    sum2 = 0
    for x in set: 
        _sum += x
        sum2 += x*x
    av_sum  = _sum/float(num)
    av_sum2 = sum2/float(num)
    tmp2 = av_sum2 - av_sum*av_sum
    #print '# av_sum=%8.3f  av_2sum = %8.3f av_sum2=%8.3f ' % (av_sum,tmp2,av_sum2)
    if fabs(tmp2)<1.0e-12 :
        var_std = 0
    elif tmp2 > 0:
        var_std = sqrt(tmp2)
    else:
        raise "Std:sqrt(x), x < 0"
    if vb: 
        dv = v_max -v_min
        if name:
            print('# Std(%s): av= %8.2f std= %8.3f max=% 8.2f min=% 8.2f  Diff(max-min)=% 8.2f' % (name,av_sum, var_std, v_max,v_min, dv))
        else:
            print('# Std: av= %8.2f std= %8.3f max=% 8.2f min=% 8.2f  Diff(max-min)=% 8.2f' % (av_sum, var_std, v_max,v_min, dv))
    if NUM:
        return num,  av_sum, var_std, v_min, v_max
    else:
        return av_sum, var_std, v_min, v_max

