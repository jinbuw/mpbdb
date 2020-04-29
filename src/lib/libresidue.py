from .libmath import Vector
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
        
    def getNT(self):
        if len(self.name) ==  1:
            return self.name
        else:
            return self.name[0]
        
    def getChainID(self):
        return self.chain_id
    
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
    
    def getResSideAtomVect(self):
        return self.getResAtomVect(SIDECHAIN=True)
        