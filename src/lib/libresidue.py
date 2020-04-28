from .libmath import cross, dihedral_vect,Vector
from numpy import array as _array


class RESIDUE(Vector):
    def __init__(self, atom, **rest):
        #self.chain_id = " "
        self.atom_objs = {}
        self.atom_list = []
        self.atom_names = []
        self.num_atoms = 0
        self.setFirstAtom(atom)
        self.aabb_atomnames = ["C", "O", "N", "CA"] 
        #self.note = note
        #self.model_id = rest.get("model_id", 0)
        #self.res_sid = rest.get("res_sid", 1)
        #print "# ** residue obj: resid= %3d res_sid= %3d"  % ( self.resid, self.res_sid)
        #print "# protein residue"
        self.DUPATOM=False
        self.dup_atoms = []
    def setChainID(self, chain_id):
        self.chain_id = chain_id
        for _atom_obj in self.atom_list:
            _atom_obj.setChainID(chain_id)
            
            
    def delAtom(self, atom_name):
        del self.atom_objs[atom_name]
        self.atom_names.remove(atom_name)
        for _obj in self.atom_list:
            if _obj.name == atom_name:
                self.atom_list.remove(_obj)
                break 
        delattr(self, atom_name)
        
        
    def setFirstAtom(self, atom):
        #print atom.__dict__.keys()
        self.chain_id = atom.chain_id
        self.resid = atom.resid
        self.name = atom.resname
        #self.addAtom(atom)
    
    def switchAtomName(self, atomn1, atomn2):
        vec1 = self.getAtomVect(atomn1)
        vec2 = self.getAtomVect(atomn2)
        atom1_obj = self.atom_objs[atomn1]
        atom2_obj = self.atom_objs[atomn2]
        atom1_obj.name = atomn2
        atom2_obj.name = atomn1 
        self.setAtomVect(atomn1, vec2)
        self.setAtomVect(atomn2, vec1)
    
    def resetAtomXYZ(self, atom_name, _atom_xyz):
        self.atom_objs[atom_name].setXYZ(xyz=_atom_xyz)
        self.setAtomVect(atom_name, _atom_xyz)
        
    def replaceAtom(self, old_atomn, new_atomn, new_xyz):
        _atom = self.atom_objs[old_atomn]
        _atom_sid = self.atom_names.index(old_atomn)
        _atom.name = new_atomn
        _atom.setXYZ(xyz=new_xyz)
        self.atom_list[_atom_sid] = _atom
        self.atom_names[_atom_sid] = new_atomn
        delattr(self, old_atomn)
        self.setAtomVect(new_atomn, new_xyz)
        
        
    def setResSid(self, res_sid):
        self.res_sid = res_sid
        
    
    def getResid(self):
        return self.resid
    
    def setResid(self, resid):
        self.resid = resid
        for _atom_obj in self.atom_list:
            _atom_obj.setResid(resid=resid)    
        
    
    def getResSid(self):
        return self.res_sid
    
    
    def getUID(self):
        """ return residue unique ID"""
        uid = "%s_%s_%s_%s" % (self.res_sid, self.resid, self.chain_id, self.model_id)
        return uid 
        
    
    def getNumAtom(self):
        return self.num_atoms
    
    
    def addAtom(self, atom):
        if atom.getResName() != "DUM" and atom.name in self.atom_objs:
            self.dup_atoms.append(atom.name)
            if not self.DUPATOM:
                #print "# Error duplicate atom name found: %s" % atom.name
                #atom.prn()
                self.DUPATOM = True
            return 
        self.atom_list.append(atom)
        self.atom_objs[atom.name] = atom 
        self.setAtomVect(atom.name, atom.axyz)
        self.atom_names.append(atom.name)
        self.num_atoms += 1
        
    
    def getAtom(self, atom_name):
        return self.atom_objs[atom_name]
    
    def setAtomVect(self, atom_name, avec):
        setattr(self, atom_name, avec)
        
    def getAtomVect(self, atom_name, vb=1):
        try:
            return getattr(self, atom_name)
        except AttributeError:
            if vb: print("#getAtomVect():  resid= %4d res_sid= %4d  could not find atom: %s  in %s " % (self.resid, self.res_sid, atom_name, self.name))
            #self.prn()
            pass
        
    
    def prn(self):
        print(self.__dict__)
    
    def __str__(self):
        return "chain_id: %s resid= %4d  resname= %4s" % (self.chain_id, self.resid, self.name )
    
        
    def getBondAngle_By_Atoms(self, atom1, atom2, atom3):
        "return bond angle Atom1 -- Atom2 -- Atom3 in deg """
        P1 = self.getAtomVect(atom1)
        P2 = self.getAtomVect(atom2)
        P3 = self.getAtomVect(atom3)
        return degrees(self.getVecAngle(P1-P2, P3-P2))
        
    def getBondVect(self, atom_name1, atom_name2):
        """ input: (atom_name1,  atom_name2) 
            return  vect  atom_name1 --> atom_name2 """
        vec1 = self.getAtomVect(atom_name1)
        vec2 = self.getAtomVect(atom_name2)
        vec = vec2 - vec1
        return vec 
    
    def getBondLen(self, *args):
        """ input: (atom_name1,  atom_name2) 
            return  vector lenght :   atom_name1 --> atom_name2  """
        vec = self.getBondVect(*args)
        vec_len = self.getVecLen(vec)
        return vec_len 
    
    def getCenter(self):
        """ get residue center  """
        if not hasattr(self, "center"):
            self.setCenter()
        return self.center
    
    
    def setCenter(self):
        res_c = _array([0., 0., 0.])
        for _atom in self.atom_list:
            res_c += _atom.axyz
        self.center = res_c / self.num_atoms
        #print "# residue %3d  num_atom= %2d center:  %8.2f, %8.2f, %8.2f" % (self.resid, self.num_atoms, res_c[0], res_c[1], res_c[2])
        
    def getSAXCenter(self):
        """ get residue center without proton  """
        res_c = _array([0., 0., 0.])
        cnt = 0
        for _atom in self.atom_list:
            if _atom.name[0] != "H":
                res_c += _atom.axyz
                cnt += 1.0
        res_c = res_c / cnt
        #print "# residue %3d  num_atom= %2d center:  %8.2f, %8.2f, %8.2f" % (self.resid, self.num_atoms, res_c[0], res_c[1], res_c[2])
        return  res_c 
    
    def getRadius(self):
        """ spherical approximation"""
        res_c = self.getCenter()
        res_sax_c = self.getSAXCenter()
        radius_sax = 0.0
        radius = 0.0 
        for _atom in self.atom_list:
            if _atom.name[0] != "H":
                vec = _atom.axyz - res_sax_c
                vec_len = getVecLen(vec)
                if vec_len > radius_sax: 
                    radius_sax  = vec_len 
                    atom_name_sax = _atom.name 
        
        for _atom in self.atom_list:
            vec = _atom.axyz - res_c
            vec_len = getVecLen(vec)    
            if vec_len > radius: 
                radius  = vec_len 
                atom_name = _atom.name 
                    
        print("# resid= %4d res_name= %4s  atom_name= %4s radius= %8.3f  sax: atom_name= %4s radius= %8.3f dif= %8.4f" %  \
              (self.resid, self.name, atom_name, radius, atom_name_sax, radius_sax, radius-radius_sax))
        return radius, radius_sax
    
    def getResXYZ(self, list=0):
        """ return residue's  cartesian coordinates """
        res_xyz = []
        for _atom in self.atom_list: 
            res_xyz.append(_atom.getXYZ())
        if list: 
            return res_xyz
        else:
            return _array(res_xyz)
    
    
    def getResXYZDt(self):
        """ return residue's  cartesian coordinates """
        res_xyz = {}
        for _atom_name in self.atom_names:
            res_xyz[_atom_name] = getattr(self, _atom_name)
        return res_xyz
    
        
    def setResXYZ(self, res_xyz):
        for i in range(len(self.atom_names)):
            self.atom_list[i].setXYZ(xyz=res_xyz[i,:])
        
 
    def getDihedraAngle_by_Atoms(self, dih_atoms ):
        """ return dihedral angle by 4 atoms of dihedral angle"""
        atom1, atom2,atom3,atom4 = dih_atoms  
        return degrees(dihedral_vect(self.getAtomVect(atom1), self.getAtomVect(atom2), self.getAtomVect(atom3), self.getAtomVect(atom4)))
    
        
    def getNT(self):
        if len(self.name) ==  1:
            return self.name
        else:
            return self.name[0]
        
    def getChainID(self):
        return self.chain_id
    
    def getResSideChainWoHCenter(self):
        """ return the center of sidechain without Hydrogen atoms"""
        res_c = _array([0., 0., 0.])
        cnt = 0
        for _atom in self.atom_list:
            if _atom.name in self.aabb_atomnames: continue
            if _atom.name in ["H","h"]: continue 
            res_c += _atom.axyz
            cnt += 1
        if cnt > 0:
            center = res_c / cnt   
        else:
            center = None
        return center
    
    
    def getResXYZ_Wo_H(self, list=0, names=None, BACKBONE=False, SIDECHAIN=False):
        """ return residue atom  cartesian coordinates  without H"""
        res_xyz = []
        #res_xyz_dt = {}
        if names is None:
            names =  self.atom_names
        for _atom_name in names: 
            if _atom_name[0] in ["H", 'h']: continue
            if BACKBONE and _atom_name not in self.aabb_atomnames: continue
            if SIDECHAIN and _atom_name in self.aabb_atomnames: continue
            _atom = getattr(self, _atom_name)
            
##            if _atom_name == "OP1":
##                try: 
##                    _atom = getattr(self, _atom_name)
##                except:
##                    _atom = getattr(self, "O1P")
##            elif _atom_name == "OP2": 
##                try: 
##                    _atom = getattr(self, _atom_name)
##                except:
##                    _atom = getattr(self, "O2P")
##            else:
##                _atom = getattr(self, _atom_name)
            res_xyz.append(_atom)
        if list: 
            return res_xyz
        else:
            return _array(res_xyz)    
 
    def getName(self):
        return self.name
    
    def setName(self, name ):
        self.name = name 
    
    def getResAtomVect(self, BACKBONE=False, SIDECHAIN=False, WOH=True):
        """ return residue atom  cartesian coordinates  without H"""
        res_xyz = []
        names =  self.atom_names
        for _atom_name in names: 
            if WOH and _atom_name[0] in ["H", 'h']: continue
            if BACKBONE and _atom_name not in self.aabb_atomnames: continue
            if SIDECHAIN and _atom_name in self.aabb_atomnames: continue
            _atom = getattr(self, _atom_name)
            res_xyz.append((_atom_name, _atom))   
        return res_xyz
    
    def getResCharge(self):
        """ return the totoal charge  of all atoms in this residue """
        tot_crg = 0.0
        for _atom_name in self.atom_names: 
            _atom = self.getAtom(_atom_name) #getattr(self, _atom_name)
            tot_crg += _atom.getCharge()  
        return tot_crg

    
    def getResSideAtomVect(self):
        return self.getResAtomVect(SIDECHAIN=True)
        