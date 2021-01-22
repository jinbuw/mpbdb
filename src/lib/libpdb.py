 
 
from numpy import array as _array 
from .libio import MIO

from .libresidue import RESIDUE as RESI
from .libmath import Vector

class ATOM:
    def __init__(self, **rest):
        self.setAtomId(**rest)
        self.setAtomName(**rest)
        self.setResid(**rest)
        self.setResName(**rest)
        self.setXYZ(**rest)
        self.setSegid(**rest)
        self.setMisc(**rest)
        self.chain_id = rest.get("chain_id", None)
        self.setCharge(**rest)
        self.setRadius(**rest)
    
    
    def setChainID(self, chain_id):
        self.chain_id = chain_id
        
        
    def setAtomName(self, **rest):
        self.name = rest.get("atomname", None ) 
    
    def setResName(self, **rest):
        self.resname = rest.get("resname", None)
    
    def setXYZ(self, **rest):
        if 'xyz' in rest:
            xyz = rest.get("xyz")
            self.x = xyz[0]
            self.y = xyz[1]
            self.z = xyz[2]
            #self.xyz = Vector(xyz)
            self.axyz = _array(xyz)
    
    def setAtomXYZ(self, xyz):
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        
    def getChainID(self):
        return self.chain_id
    
    def isH(self):
        if self.name[0] == "H":
            return True
        else:
            return False

    def getXYZ(self):
        #return #self.axyz
        return _array([self.x, self.y, self.z])
    
    def upDateXYZ(self):
        self.new_axyz = _array([self.x, self.y, self.z])
    
    def setX(self, x=None):#**rest):
        self.x = x #rest.get('x',None)
    
    def setY(self, y=None):#**rest):
        self.y = y  #rest.get('y',None)
    
    def setZ(self, z=None):#**restt):
        self.z = z #rest.get('z',None)

    def setAtomId(self, **rest):
        self.atomid = rest.get('atomid', None)
    
    def setCharge(self, **rest):
        self.atom_charge = rest.get('charge', 0)
    
    def getCharge(self):
        return self.atom_charge
    
    def setRadius(self, **rest):
        self.atom_radius = rest.get('radius', 0)
        
        
    def setResid(self, **rest):
        self.resid = rest.get('resid', None)
        
    def setSegid(self, **rest):
        self.segid = rest.get('segid', None)
        
    def setMisc(self, **rest):
        self.misc = rest.get('misc', None)
        
    def name4Pdb(self):
        """ return atom name in pdb format """
        atom_name = self.name 
        if len(atom_name) <= 3: 
            atom_name= " %s" % atom_name
        return atom_name 
    

    
    def getResid(self):
        return self.resid
    
    
    def getResName(self):
        return self.resname
    
    
    def prn(self):
        print("# name: %4s  resid= %4d chain_id= %4s " % (self.name, self.resid, self.chain_id))
        
        
class PDB(MIO, Vector): # deprecated
    #iresid = 0  # handle different chains have same residue number 
    fmt_pdb = '%-6s%5d %-4s %-4s%1s%4d%1s%3s%8.3f%8.3f%8.3f%-22s%2s'
    oneletter = {
                    'ASP':'D','GLU':'E','ASN':'N','GLN':'Q',
                    'ARG':'R','LYS':'K','PRO':'P','GLY':'G',
                    'CYS':'C','THR':'T','SER':'S','MET':'M',
                    'TRP':'W','PHE':'F','TYR':'Y','HIS':'H',
                    "HID":"H", "CYZ":"C","HIZ":"H",
                    'ALA':'A','VAL':'V','LEU':'L','ILE':'I',
                    'CIR': 'C', # mutated
                     }        
    
    letter2name = { 'D':'ASP', 'E':'GLU', 'N':'ASN', 'Q':'GLN',
                    'R':'ARG', 'K':'LYS', 'P':'PRO', 'G':'GLY',
                    'C':'CYS', 'T':'THR', 'S':'SER', 'M':'MET',
                    'W':'TRP', 'F':'PHE', 'Y':'TYR', 'H':'HIS',
                    'A':'ALA', 'V':'VAL', 'L':'LEU', 'I':'ILE'}
                 
    nt_names = {"U":"Uri", "G":"Gua", "A":"Ade", "C":"Cyt"} 
    nt_Cnames = {"U":"URI", "G":"GUA", "A":"ADE", "C":"CYT", "ADE":"ADE"}      
    aa_names = ['ALA', 'LYS', 'GLY', 'GLU', 'PHE', 'ILE', 'ARG', 'THR', 'PRO', 
                     'HIS', 'VAL', 'ASN', 'ASP', 'LEU', 'TYR', 'SER', 'CYS', 'MET', 'GLN', 'TRP']
    aa_names += ['CIR'] # mutated residues
        
    def __init__(self,  loadPDB=0, **rest):
        self.fn_pdb = rest.get('fn_pdb', None)
        self.setPDBInit(**rest)
        self.unvalid_resnames = []
        
        self.fmt_atom_pdb = '%-6s%5d %-4s %-4s%1s%4d%1s%3s%8.3f%8.3f%8.3f%-22s%2s'
          
        if loadPDB:
            self.readPDB(self.fn_pdb)
            self.setResSids()
            self.setTranResids()
            
    def loadPDBFile(self):
        self.setPDBInit()
        self.readPDB(self.fn_pdb) #fn_pdb)
        #res_names = self.getResNames()
        self.setResSids()
        
            
    def seqFasta2AANames(self, seq_fasta):
        name_list = []
        for _letter in seq_fasta:
            name_list.append(self.letter2name[_letter])
        return name_list
    
    
    def resName2Letter(self, resn):
        if len(resn) == 2:
            return resn[1]
        else:
            return self.oneletter[resn]
    
    
    def seq2fasta(self, seq, vb=0):
        seq_fasta = ""
        for _name in seq:
            seq_fasta += self.oneletter[_name]
        if vb:
            print("seq_len=%4d  fasta_len=%4d" % (len(seq), len(seq_fasta)))
            print("# fasta_seq:%s:" % seq_fasta)
        return seq_fasta            


    def getPDBSeqFasta(self):
        _seq = self.getPDBSeq()
        return self.seq2fasta(_seq)
            
            
    def setResSids(self):    
        self.num_res = len(self.res_sids)
        #print "# num_res= %4d" % self.num_res
        self.pdb_first_res_sid = self.res_sids[0]
        self.pdb_last_res_sid = self.res_sids[-1]
                
        
    def getSeqDt(self):
        seq_dt =  {}
        for res_sid in self.res_sids:
            _res =  self.residues[res_sid]
            seq_dt[res_sid] =  _res.name,  _res.resid
            pre_resid = _res.resid
        return seq_dt
    
    
    def getResNames(self):
        res_names  =  []
        for res_sid in self.res_sids:
            _res =  self.residues[res_sid]
            res_names.append(_res.name)
        return res_names
    
    
    def reset(self):
        self.setPDBInit()
        
        
    def setPDBInit(self, **rest):
        self.residues = {}
        self.res_sids = []
        self.res_sid2resid_tbl = {}
        self.resid2res_sid_tbl = {}   
        self.chain_res_sids = {} 
        self.chain_lengths = None
        self.pdb_atom_lns = []
        self.pdb_head_lns = []
        self.res_sid = 1
        self.HETATOM = rest.get("HETATOM", False)
        #print "# self.Het: ", self.HETATOM, rest 
        self.helices = {}
        self.helix_sids = []    
        self.sheets = {}
        self.sheet_sids = []
    
    def isValidResName(self, res_name):
        res_name = res_name.upper()
        if len(res_name)== 2  and res_name[0] in ["C", "N", "D"]:
            return True
        elif res_name == "DUM":
            return True
        elif res_name not in self.aa_names:
            return False
        else:
            return True
    
    def setChainLengths(self):
        if self.chain_lengths is None:
            self.chain_lengths = {}
            for chain_id in self.chain_ids:
                self.chain_lengths[chain_id] = len(self.chain_res_sids[chain_id])
    
    def resetChainLengths(self):
        self.chain_lengths = {}
        for chain_id in self.chain_ids:
            self.chain_lengths[chain_id] = len(self.chain_res_sids[chain_id])

        
    def getChainLengths(self):
        if self.chain_lengths is None:
            self.setChainLengths()
        return self.chain_lengths
    
    
    def getChainIDs(self):
        return self.chain_ids 
    
    
    def getChainLength(self, chain_id):
        if self.chain_lengths is None:
            self.setChainLengths()
        #print("# chain_lengths: {}".format(self.chain_lengths))
        return self.chain_lengths.get(chain_id, 0)
        
    def setTranResids(self, tran_resids ={} ):
        self.TRANS_RESIDS  = tran_resids
                      
    
    def synResid(self):
        new_resids = []
        new_seq_dt = {}
        new_residues = {}
        #print '# selfresis:', self.resids 
        for resid in list(self.residues.keys()):
            res = self.residues[resid]
            new_resids.append(res.resid)
            new_seq_dt[res.resid] = res.nt_name 
            new_residues[res.resid] = self.residues[resid]
        self.resids = new_resids
        self.seq_dict = new_seq_dt
        self.residues = new_residues
        
        
    def getTranResid(self, resid):
        return self.TRANS_RESIDS.get(resid, resid)
    
    
    def readPDB(self, fn_pdb, newpdb=True, RESI=RESI, vb=0):
        if vb: print("# reading pdb: %s " % fn_pdb)
        self.pdb_lns = self.readTxtFile(fn_pdb) 
        self.parsePDB(self.pdb_lns, resi=RESI, newpdb=newpdb)
        self.chain_ids  = list(self.chain_res_sids.keys())
        self.chain_ids.sort()
        self.num_chain = len(self.chain_ids)
        self.num_res = len(self.res_sids)
        if vb: print("# num_chains=%2d  num_res: %d" %  (self.num_chain, self.num_res))            


    def loadPDBFile(self, fn_pdb):
        self.setPDBInit()
        self.readPDB(fn_pdb) #fn_pdb)
        #res_names = self.getResNames()
        self.setResSids()
        

    def prnChainInf(self, prn=1):
        lns = [""]
        ln  =  "%12s, num_chain=%2d" % (self.getBaseName(self.fn_pdb),self.num_chain)
        for _chain_id in self.chain_ids:
            ln += ", %2s:%2d" % (_chain_id, len(self.chain_res_sids[_chain_id]))
        lns.append( ln )
        self.prnLines(lns, prn=prn)
        self.add2log(lns)
        return lns
    
    def isMHC(self):
        """ define MHC by chain lengths which is not right."""
        AC = False
        BC = False
        CC = False
        a_id = None
        b_id = None
        c_id = None
        for chain_id in self.chain_ids:
            chain_len = len(self.chain_res_sids[chain_id])
            if chain_len > 270 and chain_len < 280:
                if not AC:
                    AC= True
                    a_id = chain_id
            if chain_len >90 and chain_len < 110:
                if not BC:
                    BC=True
                    b_id = chain_id
            if chain_len >7 and chain_len < 10:
                if not CC:
                    CC=True
                    c_id = chain_id
        self.mhc_chain_ids = [a_id, b_id, c_id]        
        if AC and BC and CC:
            return True
        else:
            return False
        
        
    def getMHCRessids(self):
        res_sids = []
        for chain_id in self.mhc_chain_ids:
            res_sids += self.chain_res_sids[chain_id]
        res_sids.sort()
        return res_sids
    
    def resSid2Resid(self, res_sid, vb=0): #v2
        resid = self.res_sid2resid_tbl.get(res_sid, res_sid)
        return resid
  
    
    
    def getRes2ResCloseAtoms(self, res1, res2, cutoff2=25.0, SAMEATOM=True):
        close_atoms = []
        for atn1 in res1.atom_names:
            atv1 = res1.getAtomVect(atn1)
            for atn2 in res2.atom_names:
                if SAMEATOM and atn1[0] == atn2[0]: continue # exclude same atom type 
                atv2 = res2.getAtomVect(atn2)
                dr2 = self.getDistSquare(atv1, atv2)
                if dr2 < cutoff2:
                    close_atoms.append((atn1, atn2, sqrt(dr2)))
        return close_atoms    


    def getRes2ResClosestAtom(self, res1, res2, cutoff2=25.0):
        #close_atoms = []
        min_dr2 = 9999.0
        for atn1 in res1.atom_names:
            atv1 = res1.getAtomVect(atn1)
            for atn2 in res2.atom_names:
                if atn1[0] == atn2[0]: continue # exclude same atom type 
                atv2 = res2.getAtomVect(atn2)
                dr2 = self.getDistSquare(atv1, atv2)
                if min_dr2 > dr2:
                    min_dr2 = dr2 
                    min_p = atn1, atn2, sqrt(dr2)
                    #close_atoms.append((atn1, atn2, sqrt(dr2)))
        return min_p #close_atoms    


    def resid2ResSid(self, resid, chain_id=None): #v2
        #print self.resid2res_sid_tbl
        res_sids = self.resid2res_sid_tbl[resid]
        #print self.resid2res_sid_tbl
        if len(res_sids) > 1:
            for _ressid, _chain_id in res_sids:
                if _chain_id == chain_id:
                    return _ressid
            print("# Warning: 2+ residues have same residue name in the pdb file, only return first resid in the pdb file ")    
        return res_sids[0][0]
        

    def getPdbCode(self):
        fn_pdb = os.path.basename(self.fn_pdb)
        pdb_code, fn_ext = fn_pdb.rsplit(".", 1) #split("\.", fn_pdb)
        #print split("\.", fn_pdb)
        return pdb_code 

    
    def resid2Res(self, resid, chain_id=None): #v2
        #print  "*** resid=",  resid 
        res_sid = self.resid2ResSid(resid, chain_id=chain_id)
        return self.resSid2Res(res_sid)
    
    
    def resSid2Res(self, res_sid): #v2 
        return self.residues[res_sid]
    
    
    def getRes(self, res_sid):
        return self.residues[res_sid]

    
    def AtomName4Pdb(self, atom_name):
        """ return atom name in pdb format """
        if len(atom_name) <= 3: 
            atom_name= " %s" % atom_name
        return atom_name 
    
    def resids2segs(self, resids, break_resids=[]):
        resids.sort()
        segs = []
        _seg = [resids[0]]
        num_res = len(resids)
        #print "# tot_res:%4d " % num_res 
        for i in range(1, num_res):
            pre_resid = resids[i-1]
            resid = resids[i]
            if resid in break_resids:
                segs.append(_seg)
                _seg = [resid]
            elif resid - pre_resid == 1:
                _seg.append(resid)
            else:
                segs.append(_seg)
                _seg = [resid]
        segs.append(_seg)
        return segs          
    
    def hasDupRes(self, vb=0):
        has_dup_res = False
        for res_sid in self.res_sids:
            _res = self.resSid2Res(res_sid)
            if _res.DUPATOM:
                if vb: print("# dupres: resid: %s   resn: %s   dupatoms: %s " % (_res.resid, _res.name, _res.dup_atoms))
                has_dup_res = True
                break 
        return has_dup_res
    
    
    def getSimpResidStr(self, resids):
        segs = self.resids2segs(resids)
        ln = ""
        for _seg in segs:
            num_elem = len(_seg)
            if num_elem >1:
                ln += "%s-%s, " % (_seg[0], _seg[-1])
##            elif num_elem > 1:
##                ln +="%s-%s, " % (_seg[0], _seg[1])
            else:
                ln += "%s, " % _seg[0]
        return ln 
    
    def isConnectedRes(self, pre_res, next_res):
        if pre_res is None: return False
        C_vec =  pre_res.getAtomVect("C", vb=0)        
        N_vec = next_res.getAtomVect("N", vb=0)
        if C_vec is None or N_vec is None:
            return True
        dist = self.getAtomDist(C_vec, N_vec)
        #print "# %4d C -- %4d N dist= %8.1f"  % (pre_res.resid, next_res.resid, dist)
        if dist < 1.5: 
            return True
        else:
            print("# chain: %s %4d/%4d C -- chain: %s %4d/%4d N dist= %8.1f"  % (pre_res.chain_id, pre_res.res_sid, pre_res.resid,  
                                                                   next_res.chain_id,next_res.res_sid, next_res.resid, dist))
            print("# pre_C: ", C_vec)
            print("# nextN: ", N_vec)
            return False
    
    
    def prnSeq(self, cols=10, id=None):
        res_sids = list(self.residues.keys())
        res_sids.sort()
        i = 0 
        ln = ""
        cnt = 0
        print("# seq: ")
        if cols > 1:
            while i < self.num_res:
                res = self.residues[res_sids[i]]
                ln += "%4d %4s " % (res.resid, res.name)
                cnt += 1
                i+= 1
                if cnt == cols:
                    cnt = 0 
                    print(ln) 
                    ln = ""
            print(ln)
        else:
            
            while i < self.num_res:
                res = self.residues[res_sids[i]]
                ln += "%s" %  res.name
                i += 1
                cnt += 1
            print("%s:%s" %  (cnt, ln ))
    
                
    def atomTxt2Atom(self, ln):
        """ parse atom text to atom object:
        Cols.  
        1-6    Record name "ATOM  " or "HETATM"
        7-11   Atom serial number                   (see note i)
        13-14   Chemical symbol (right justified)  )
        15      Remoteness indicator               ) (see note ii)
        16      Branch designator                  )
        17      Alternate location indicator         (see note iii)
        18-20   Residue name                         (see note iv)
        21      Reserved                   )
        22      Chain identifier           )
                                           )         (see note v)
        23-26   Residue sequence number    )
        27      Code for inserting residue )
        31-38   X   )
        39-46   Y   ) Orthogonal Angstrom coordinates
        47-54   Z   )
        55-60   Occupancy
        61-66   Isotropic B-factor
        73-76   Segment identifier, left justified (used by XPLOR)
        77-78   Element symbol, right justified )
                                                )    (see note vi)
        79-80   Charge on atom                  )
        """
        try:
            resid = int(ln[22:26])
        except ValueError:
            return 
        ln = ln.lstrip() 
        try:
            atomid = int(ln[6:11])
        except ValueError:
            atomid = -999
            
        atomname= ln[12:16].strip()
        resname = ln[17:20]
        x = float(ln[30:38])
        y = float(ln[38:46])
        z = float(ln[46:54])
        misc = ln[54:72]
        segid =  ln[72:76].rstrip()
        chain_id = ln[21]
        
        try:
            atom_radius = float(ln[54:60])
        except ValueError:
            atom_radius = 0
            
            
        try:
            atom_charge = float(ln[60:67])
        except ValueError:
            atom_charge = 0
            
        atom = ATOM(atomname=atomname.strip(), atomid=atomid, xyz=(x,y,z), resid=resid,  resname=resname.strip(), misc=misc, \
                    segid=segid, chain_id=chain_id, charge=atom_charge, raius=atom_radius)
        
        #atom = ATOM(atomname=atomname.strip(), atomid=atomid, xyz=(x,y,z), resid=resid,  resname=resname.strip(), misc=misc, segid=segid, chain_id=chain_id)
        atom.ln = ln 
        #atom.prn()
        return atom         
    
    
##    def getAminoNames(self):
##        """ use 1TTT.pdb generate 20 aminoacid names """
##        aa_names  = []
##        for resid in range(229, 1440):
##            if self.residues[resid].name not in aa_names:
##                aa_names.append(self.residues[resid].name)
##        print aa_names
##        print "# num_aa: %3d " % len(aa_names)
        
        
    def parsePDB(self, pdb_lns, resi=None, newpdb=True): #v2
        """ parse pdb datafile to get residues, sequence, cartesian coordinates
        """
        residue  = None
        pre_resid = -9999
        pre_resname = None
        body = 0
        num_ter = 0
        num_res = 0
        #print "# num_lns: %6d"  % len(pdb_lns)
        #print "# het: ", self.HETATOM
        for ln in pdb_lns:
            #if ln.startswith("REMA") and body == 0:
            #print("pdb: %s" % ln) 
            if ln.startswith("ATOM") or (self.HETATOM and ln.startswith("HETATM")):
                body = 1
                self.pdb_atom_lns.append(ln)
                atom = self.atomTxt2Atom(ln)
                if atom is None: 
                    print("#*** Warning: non recognized atom: %s" % ln) 
                    continue
                #print ln 
                resname = atom.getResName()
                if not self.isValidResName(resname): 
                    if resname not in self.unvalid_resnames:
                        print("# Error: not a valid residue name: %s : %s " % (resname, ln.strip()))
                        self.unvalid_resnames.append(resname)
                    continue
                chain_id = atom.getChainID()
                #if chain_id == "C":
                #    print ln 
                resid = atom.getResid()
                if pre_resid != resid or (pre_resname != resname): 
                    self.addResidue(residue)
                    num_res += 1
                    pre_resid = resid
                    pre_resname = resname                    
                    residue = resi(atom) #RESI(atom)
                residue.addAtom(atom)
                self.pdb_atom_lns.append(ln)
                #print ln 
            elif ln.startswith("TER") and num_ter:
                num_ter += 1
            elif ln.startswith('ENDMDL'): 
                #print "# find endmod "
                break
            elif body == 0:
                self.pdb_head_lns.append(ln)
        self.addResidue(residue)
        num_res += 1
        
        #print "#num_res= %4d" %  num_res
        return num_res
                  
    
    def getTotEnergy(self):
        if hasattr(self, "totE"):
            return self.totE
        else:
            for ln in self.pdb_head_lns:
                if ln.startswith("REMARK summary total"):
                    tmp = ln.split()
                    self.totE = float(tmp[3])
                    return self.totE
                
                
    def parseXPRDCOut(self):
        for ln in self.pdb_lns:
            if ln.startswith("REMARK RDC RDC"):
                tmp = ln.split()
                self.rdc_xp_rms = float(tmp[3])
                break
        
    
    def getXpRdcRms(self, vb=0):
        if not hasattr(self, "rdc_xp_rms"):
            self.parseXPRDCOut()
        if vb: 
            print("# RDCV rms: %7.2f " % self.rdc_xp_rms)
        return self.rdc_xp_rms
    
        
    def addResSid2Chain(self, res_sid, chain_id):
        if chain_id in self.chain_res_sids:
            self.chain_res_sids[chain_id].append(res_sid)
        else:
            self.chain_res_sids[chain_id] = [res_sid]        
    
    
    def addResidue(self, residue): #v2
        if residue is not None:
            if residue.name == "ANI": return 
            #if residue.getNumAtom() > 4:
            residue.setResSid(self.res_sid)
            self.res_sids.append(self.res_sid)
            res_chain_id = residue.getChainID()
            self.res_sid2resid_tbl[self.res_sid] = residue.getResid()
            self.addResid2Tbl(residue.getResid(), self.res_sid,res_chain_id)
            self.addResSid2Chain(self.res_sid, res_chain_id)
            self.residues[self.res_sid] = residue
            self.res_sid += 1
    

    def addResid2Tbl(self, resid, res_sid, chain_id): #v2
        if resid in self.resid2res_sid_tbl:
            self.resid2res_sid_tbl[resid].append((res_sid, chain_id))
        else:
            self.resid2res_sid_tbl[resid] = [(res_sid, chain_id)]
    
    
    def getPDBCode(self, vb=0):
        fn_pdb = os.path.basename(self.fn_pdb)
        pdb_code, tt = fn_pdb.split(".")
        if vb: print("# pdb_code(or fn_pre): %s "  % pdb_code)
        return pdb_code 
    
    
    def parseRNA_in_ComboPDB(self):
        """ parse protein or/and RNA pdb files to get residues, sequence, cartesian coordinates
        """
        self.residues = {}
        self.pdb_atom_lns = []
        residue  = None
        pre_resid = -9999
        body = 0
        self.pdb_misc_lns = []
        atoms = []
        for ln in self.pdb_lns:
            if ln[0:4] == "REMA" and body == 0:
                self.pdb_head_lns.append(ln)
            elif ln[0:4] == "ATOM" or ln[0:6] in 'HETATM':
                body = 1
                atom = self.atomTxt2Atom(ln)
                atoms.append(atom)
                resid = atom.resid
                if pre_resid != resid: 
                    if len(atoms) > 3:  # new residue
                        self.addResidue_By_Atoms(atoms)
                    else: # metal ion atoms or HOH
                        for _atom in atoms:
                            self.pdb_misc_lns.append(_atom.ln)
                    atoms = []
                    pre_resid = resid
            else:
                self.pdb_misc_lns.append(ln)
        
        if len(atoms) > 3:   # last residue
            self.addResidue_By_Atoms(atoms)
        self.setPDBresids()
    
    
    
    def getRNASeq(self):
        rna_resids = []
        rna_seqs = ""
        for res_sid in self.res_sids:
            
            _resname = self.residues[res_sid].name
            if _resname not in self.aa_names and _resname[-1] in ['A','U','C','G']:
                resid =  self.resSid2Resid(res_sid)
                rna_resids.append(resid)
                nt_name = _resname[-1]
                #rna_seqs.append(nt_name)
                rna_seqs += nt_name
        #print "# rna_seq: ", rna_seqs
        #print "# rna_resids: ", rna_resids 
        return rna_seqs, rna_resids
            
            
##    def setPDBresSids(self):
##        """ set pdb resids from residue objects """
##        self.resSids = self.residues.keys()
##        self.resids.sort()
##        self.num_res = len(self.residues)
##        #self.last_atomid = atomid
        
        
    def getNumRes(self):
        """ return the number of residues in the pdb file """
        return self.num_res 
                  
    
    def prnResName(self):
        print("# pdb: cutinase ") 
        for resid in self.resids:
            print("%4d %8s " % (resid, self.residues[resid].name))
    
    def prnResSids(self):
        cnt = 0 
        for res_sid in self.res_sids:
            cnt += 1
            if cnt > 10  and  cnt < 300: continue
            _res = self.resSid2Res(res_sid)
            print("# sid= %4s    resn: %s  resid= %s  chainID: %s" %  (res_sid, _res.name,  _res.resid, _res.chain_id))


    def prnAAtomXYZ(self, res_first, res_end, vb=1):
        # array index starts from 0, but the amino residues from 1.
        # frist resudue's N atom is not included in the pdb file generated by gen_pdb function
        resxyz, atom_names = self.getResXYZ(res_first, res_end)
        xyz_lns = []
        for i in range(len(resxyz)):
            #print "%3d %3s %8.4f %8.4f %8.4f" % (i+1, atom_names[i],xyz[i].x(), xyz[i].y(), xyz[i].z())
            ln =  "%3d %3s %8.3f %8.3f %8.3f" % (atom_names[i][0],atom_names[i][1], resxyz[i,0], resxyz[i,1], resxyz[i,2])
            xyz_lns.append(ln)
            if vb: print(ln) 
        return resxyz, xyz_lns
    
  
    
    def getAtomType(self, atom_name):
        if atom_name[0] in digits:
            return atom_name[1]
        else:
            return atom_name[0]
        
        
    
    
    def __addDumAtomLns(self):
        if hasattr(self, "dum_atom_xyzs"):
            lns = []
            dummy_atoms = ["C","CA"]
            _xyzs = self.dum_atom_xyzs
            for i in range(len(_xyzs)):
                ln =  self.fmt_atom_pdb  % ('ATOM', 0, dummy_atoms[i], "DUM",  " ", -1, "","", _xyzs[i][0], _xyzs[i][1], _xyzs[i][2],  "", "C")
                lns.append(ln)
            self.pdb_lns = lns + self.pdb_lns
        #return lns 
                       
    def __setInitWritePDB(self, lns):
        self.fmt_pdb = '%-6s%5d %-4s %-4s%1s%4d%1s%3s%8.3f%8.3f%8.3f%-22s%2s'
        if lns is None: 
            self.pdb_lns = []
        else:
            self.pdb_lns = lns 
        self.serial_atom_id = 1

        
    def writePDB(self, fn_pdb, lns=None, vb=1, resids=None, shift_resid=0, RESID=False, res_sid=False, \
                 atom_names=None, new_atomid=False, tran_resid=False,uniq_resid=False, uniq_res=False, \
                 TER=False,chain_id=None, res_sids=None,REMARK=True, WITHEND=True, CN=False,DMA=False):
        """ DMA: dummy major atom for overcoming the grid size problem in delphi computation """
        
        self.__setInitWritePDB(lns)
        res_sids, num_res = self.__getPDBResSids(res_sids)
        if vb: print("# saving: num_res= %4d " % num_res)
        for i in range(num_res):
            res_sid = res_sids[i]
            resid = self.resSid2Resid(res_sid)
            if resids is not None and resid not in resids: continue
            _amino = self.resSid2Res(res_sid)
            _resid  = self.__getPDBResid(res_sid, uniq_resid, tran_resid, _amino, RESID, resid, shift_resid)
            resn = self.__getPDBResn(res_sid, _amino, CN)
            self._res_atom_names = []
            for _atom in _amino.atom_list:
                self.__addAtomPDBLn(_atom, atom_names, uniq_res, _resid, new_atomid, chain_id, resn)
                
            self.__addTer2ResBreak(TER, i)
        self.__addDumAtomLns()
        self.__addTER(TER)
        self.__addREMARK( REMARK)
        self.__addEND( WITHEND)
        self.__writePDBLns(fn_pdb, vb=1)
        return self.pdb_lns 
    
    def __addAtomPDBLn(self, _atom, atom_names, uniq_res, _resid, new_atomid, chain_id, resn):
        atom_name = _atom.name 
        
        if atom_name in ["H3T", "H5T"]: return   
        if uniq_res and atom_name in self._res_atom_names: return 
        
        self._res_atom_names.append(atom_name)
        if atom_names is not None and atom_name not in atom_names: return 
        atom_type = self.getAtomType(atom_name)
        if len(atom_name) <= 3: atom_name= " %s" % atom_name
     
        if _atom.misc is None: 
            misc = "  1.00  0.00"
        else:
            misc  = _atom.misc
        #misc = "  1.00  0.00"
     
        if _atom.segid is None:
            segid = ""
        else:
            segid = _atom.segid 
      
        
        if new_atomid:
            _atom_id = self.serial_atom_id
            self.serial_atom_id += 1
        else:
            _atom_id = _atom.atomid
     
        if chain_id is None:
            _chain_id = _atom.chain_id
        else:
            _chain_id = chain_id
            
        ln =  self.fmt_pdb  % ('ATOM', _atom_id, atom_name, resn, _chain_id, _resid, "","", _atom.x, _atom.y, _atom.z,  misc, atom_type)
        self.pdb_lns.append(ln)
        #return ln 
                
    
    def __getPDBResSids(self, res_sids):
        if res_sids is None: res_sids = self.res_sids
        self.uniq_resid_tbl = []
        self.first_res_sid = res_sids[0]
        self.last_res_sid = res_sids[-1]
        num_res = len(res_sids)
        self.num_res_1 = num_res - 1
        return res_sids, num_res
    
    
    def __addTer2ResBreak(self, TER, i):
        if TER:
            if i >= self.num_res_1: 
                next_res_sid = res_sid
            else:
                next_res_sid = res_sids[i+1]
            next_resid = self.resSid2Resid(next_res_sid)     
                #print "# next_resid= %4d  resid= %4d " % (next_resid, resid)
            if next_resid - resid > 1:
                self.pdb_lns.append( '%-6s%5d %-4s %-4s%1s%4d' % ("TER", _atom_id+1, "", resn, chain_id, _resid))    
        
        
    def __addTER(self, TER):
        if TER and  not lns[-1].startswith("TER"):
            self.pdb_lns.append( '%-6s%5d %-4s %-4s%1s%4d' % ("TER", _atom_id+1, "", resn, chain_id, _resid))    
        
        
    def __addREMARK(self, REMARK):
        if REMARK:
            self.pdb_lns =  self.pdb_head_lns + self.pdb_lns + ["END"]
        
        
    def __addEND(self, WITHEND):
        if WITHEND:
            self.pdb_lns =  self.pdb_lns + ["END"]
        
    def __writePDBLns(self, fn_pdb, vb=0):
        if fn_pdb is not None:
            if vb: print("# pdb was saved in < %s >: " %  fn_pdb)
            self.saveLines(self.pdb_lns, fn_pdb, vb=vb)
            if self.uniq_resid_tbl:
                self.saveLines(self.uniq_resid_tbl, "%s.tbl" % fn_pdb[:-4], vb=1)
        else:
            self.err_exit( " no output pdb filename: %s" % fn_pdb)
        
        
    def __getPDBResid(self, res_sid, uniq_resid, tran_resid, _amino, RESID, resid, shift_resid):
        if uniq_resid:
            _resid = res_sid
            self.uniq_resid_tbl.append("%8d %s%d" %  (_resid, _amino.chain_id, _amino.resid))
        elif tran_resid:
            _resid = self.getTranResid(_amino.resid)
        elif RESID:
            _resid = _amino.resid
        else:
            _resid = shift_resid + resid  
        return _resid
    
            
    def getNormResn(self, resn):
        resn_len = len(resn)
        if resn_len == 3:
            return resn
        elif resn_len == 2:
            return self.letter2name[resn[1]]
        elif resn_len == 4:
            return resn[1:]
        else:
            print("# Error in getNormResn(): resn: %s"  % resn) 
            
    def __getPDBResn(self, res_sid, _amino, CN):
        #print "# PARSEFPM:", self.PARSEFPM
        if hasattr(self, "PARSEFPM") and self.PARSEFPM:
            resn = self.getNormResn(_amino.name)
        elif CN:
            if hasattr(self, "null_resids") and res_sid in self.null_resids:
                resn = " D%s" % self.resName2Letter(_amino.name)
            elif res_sid == self.first_res_sid:
                resn = " N%s" % self.resName2Letter(_amino.name)
            elif res_sid == self.last_res_sid:
                resn = " C%s" % self.resName2Letter(_amino.name)
            else:
                resn = _amino.name 
        else:
            resn = _amino.name 
        return resn 
        
    ### new function 01/29/2020 
        
    def writeDelphiPDB(self, fn_pdb, lns=None, vb=1, resids=None, shift_resid=0, RESID=False, res_sid=False, \
                 atom_names=None, new_atomid=False, tran_resid=False,uniq_resid=False, uniq_res=False, \
                 TER=False,chain_id=None, res_sids=None,REMARK=True, WITHEND=True, CN=False,DMA=False):
        """ write delphi pdb for delphi computation  """

      
        self.__setInitWritePDB(lns)
        res_sids, num_res = self.__getPDBResSids(res_sids)
        for i in range(num_res):
            res_sid = res_sids[i]
            sresid = self.resSid2Resid(res_sid)
            if resids is not None and resid not in resids: continue
            _amino = self.resSid2Res(res_sid)
            _resid  = self.__getPDBResid(res_sid, uniq_resid, tran_resid, _amino, RESID, sresid, shift_resid)

            if CN:
                resn = self.pdb_dp_resns[res_sid]
            else:
                resn = _amino.name 
                
            self._res_atom_names = []
            for _atom in _amino.atom_list:
                self.__addAtomPDBLn(_atom, atom_names, uniq_res, _resid, new_atomid, chain_id, resn)
            self.__addTer2ResBreak(TER, i)
        self.__addDumAtomLns()
        self.__addTER(TER)
        self.__addREMARK( REMARK)
        self.__addEND( WITHEND)
        self.__writePDBLns(fn_pdb, vb=1)
        return self.pdb_lns 
        
 
        
    def savePDB_Chains(self, fn_pdb, chainids,  lns=None, vb=1, resids=None, shift_resid=0, res_sid=False, \
                 atom_names=None, new_atomid=False, tran_resid=False,uniq_resid=False,TER=False,  res_sids=None):
        """ save pdb residues by their chain ids """
       
        if chainids is None: 
            print("# Error:  No input chain ids, no pdb will be saved")
            return 
        if lns is None: lns = [] 
        lns += self.pdb_head_lns
        serial_atom_id = 1
        for chain_id in chainids:
            res_sids = self.chain_res_sids[chain_id]
            #print "# chain: %s  res_sids: %s"  % (chain_id, res_sids)
            for _res_sid in res_sids:
                #print "# _res_sid=%s"  %  _res_sid
                _amino = self.resSid2Res(_res_sid)
                resid = _amino.resid
                for _atom in _amino.atom_list:
                    atom_name = _atom.name 
                    if atom_name in ["H3T", "H5T"]: continue  
                    atom_type = self.getAtomType(atom_name)
                    if len(atom_name) <= 3: atom_name= " %s" % atom_name
                    if _atom.misc is None: 
                        misc = "  1.00  0.00"
                    else:
                        misc  = _atom.misc
                    misc = "  1.00  0.00"
                    if _atom.segid is None:
                        segid = ""
                    else:
                        segid = _atom.segid 
                    #segid = atom_name[0]
                    if uniq_resid:
                        _resid = res_sid
                    elif tran_resid:
                        _resid = self.getTranResid(_atom.resid)
                    else:
                        _resid = shift_resid + resid 
                    if new_atomid:
                        _atom_id = serial_atom_id
                        serial_atom_id += 1
                    else:
                        _atom_id = _atom.atomid
                    if len(_amino.name) == 3:
                        resn = _amino.name
                    else:
                        resn = "  " + _amino.name
                    if chain_id is None:
                        _chain_id = _atom.chain_id
                    else:
                        _chain_id = chain_id
                        
                    ln =  self.fmt_pdb  % ('ATOM', _atom_id, atom_name, resn, _chain_id, _resid, "","", _atom.x, _atom.y, _atom.z,  misc, atom_type)
                    lns.append(ln)
        lns =  self.pdb_head_lns + lns + ["END"]
        self.saveLines(lns, fn_pdb, vb=vb)

        
    def getPDBTXT(self, Ter=False, vb=1, res_sids=None, resids=None, shift_resid=0, end=True, res_sid=False, atom_names=None, new_atomid=False, tran_resid=False,uniq_resid=False):
        fmt_pdb = '%-6s%5d %-4s %-4s%1s%4d%1s%3s%8.3f%8.3f%8.3f%-22s%2s'
        cnt = 0
        lns = []
        serial_atom_id = 0
        #print "# num_resid= %4d " % len(self.res_sids)
        if res_sids is None:
            res_sids = self. res_sids
        for res_sid in res_sids:
            resid = self.resSid2Resid(res_sid)
            #print "# res_sid= %4d resid= %4d "% (res_sid, resid)
            if resids is not None and resid not in resids: continue
            _amino = self.resSid2Res(res_sid)
           
            for _atom in _amino.atom_list:
                atom_name = _atom.name 
                if atom_name in ["H3T", "H5T"]: continue  
             
                #if atom_name == "'HO2": continue # DNA
                atom_type = self.getAtomType(atom_name)
                if len(atom_name) <= 3: atom_name= " %s" % atom_name
                if _atom.misc is None: 
                    misc = "  1.00  0.00"
                else:
                    misc  = _atom.misc
                misc = "  1.00  0.00"
                if _atom.segid is None:
                    segid = ""
                else:
                    segid = _atom.segid 
                #segid = atom_name[0]
                if uniq_resid:
                    _resid = res_sid
                elif tran_resid:
                    _resid = self.getTranResid(_atom.resid)
                else:
                    _resid = shift_resid + resid 
                
                #ln =  fmt_pdb  % ('ATOM', _atom.atomid, atom_name, _atom.resname," ", _atom.resid, "","", _atom.x, _atom.y, _atom.z,  misc, segid)
                if new_atomid:
                    _atom_id = serial_atom_id
                    serial_atom_id += 1
                else:
                    _atom_id = _atom.atomid
                    #print "# atomid in atom: ", _atom_id
                #print '# segid: ', segid, "misc: ", misc  
                if len(_amino.name) == 3:
                    resn = _amino.name
                else:
                    resn = "  " + _amino.name
                ln =  fmt_pdb  % ('ATOM', _atom_id, atom_name, resn," ", _resid, "","", _atom.x, _atom.y, _atom.z,  misc, atom_type)
                lns.append(ln)
                
        if Ter:
            lns.append( '%-6s%5d %-4s %-4s%1s%4d' % ("TER", _atom_id+1, "", resn, "", _resid))
        elif end:
            lns =  lns + ["END"]
        else:
            pass
        
        return lns    
            
            
    def addPDB_REMARK(self, lns):
        new_lns = []
        for _ln in lns:
            new_lns.append("REMARK " + _ln)
        return new_lns 
        
    def setPDBXYZ(self): 
        pdb_res_iatom_ids = []
        pdb_res_eatom_ids = []
        cnt  = 0    
        pdb_xyz = []
        for resid in self.res_sids:
            pdb_res_iatom_ids.append(cnt)
            for _atom in self.residues[resid].atom_list:
                pdb_xyz.append(copy(_atom.axyz))
                cnt += 1
            pdb_res_eatom_ids.append(cnt-1)
        #print pdb_res_eatom_ids
        #print pdb_res_iatom_ids
        self.pdb_xyz_arr = _array(pdb_xyz)
        self.pres_ia_ids_arr = _array(pdb_res_iatom_ids, dtype=Npy.int32)
        self.pres_ie_ids_arr = _array(pdb_res_eatom_ids, dtype=Npy.int32)
    
    def setPDBXYZwoH(self): 
        #pdb_res_iatom_ids = []
        #pdb_res_eatom_ids = []
        cnt  = 0    
        pdb_xyz = []
        for resid in self.res_sids:
            #pdb_res_iatom_ids.append(cnt)
            for _atom in self.residues[resid].atom_list:
                if _atom.isH():
                    continue
                pdb_xyz.append(copy(_atom.axyz))
                cnt += 1
            #pdb_res_eatom_ids.append(cnt-1)
        #print pdb_res_eatom_ids
        #print pdb_res_iatom_ids
        self.pdb_xyz_arr_woH = _array(pdb_xyz)
        #self.pres_ia_ids_arr = _array(pdb_res_iatom_ids, dtype=Npy.int32)
        #self.pres_ie_ids_arr = _array(pdb_res_eatom_ids, dtype=Npy.int32)
        
        
    def getMajorDim(self):
        """ find the major direction of pdb and return minimum and maximum coordinates of this direction """
        pdb_xyz = self.getPDBXYZwoH()
        xyz_max = Npy.amax(pdb_xyz, axis=0)
        xyz_min = Npy.amin(pdb_xyz, axis=0)
        dxyz = xyz_max - xyz_min
        md = Npy.argmax(dxyz)
        #print md, dxyz, xyz_min, xyz_max 
        return md, xyz_min[md],xyz_max[md]
    
        
    def getPDBXYZ(self):
        if not hasattr(self, "pdb_xyz_arr"):
            self.setPDBXYZ()
        return self.pdb_xyz_arr
    
    def getPDBXYZwoH(self):
        if not hasattr(self, "pdb_xyz_arr"):
            self.setPDBXYZwoH()
        return self.pdb_xyz_arr_woH
    

    def getPDBCenter(self):
        pdb_xyz_arr = self.getPDBXYZ()
        return Npy.sum(pdb_xyz_arr, axis=0)/len(pdb_xyz_arr)
    
        
    def upDatePDBXYZ(self, pdbxyz):
        cnt = 0
        for resid in self.res_sids:
            for _atom in self.residues[resid].atom_list:
                _atom.setAtomXYZ(pdbxyz[cnt])
                cnt += 1
    
    
    def setResCenter(self): # pass 
        """ set Monomer A SAX atoms residue centers, Monomer A geometry center as the Cartesian origin point(0, 0, 0) """
        self.res_centers = {}
        for resid in list(self.residues.keys()):
            self.res_centers[resid] = self.residues[resid].getCenter()
            
    
    def getMatchedRes(self, chain_id, ref_vec):
        res_sids = self.getChainResSids(chain_id)
        min_dr = 999.0 
        min_res_sid = -1 
        for res_sid in res_sids:
            vec = self.getResAtomVect(res_sid, "CA")
            dr = self.getAtomDist(vec, ref_vec)                    
            if dr < min_dr:
                min_dr = dr 
                min_res_sid = res_sid
        return min_dr, self.residues[min_res_sid]
    
        
    def getAtomDist(self, ca, cb):
        dx = ca[0] - cb[0]
        dy = ca[1] - cb[1]
        dz = ca[2] - cb[2]
        dr2 = dx * dx + dy*dy + dz*dz 
        dr = sqrt(dr2)
        return dr 
    
    def getResAtomVect(self, res_sid, atom, vb=0):
        _res = self.residues[res_sid]
        vec =  _res.getAtomVect(atom, vb=0)
        if vb:
            if vec is None:
                print("#Warning:  res: %4d has not atom: %s " % (res_sid, atom ))
            else:
                print("# res_sid= %4d  atom: %4s   xyz: %s " % (res_sid, atom, vec))
        return vec
    
    
    def getAnyAtomDist(self, res_sid1, atom1, res_sid2, atom2, vb=0):
        vec1 = self.getResAtomVect(res_sid1, atom1, vb=vb)
        vec2 = self.getResAtomVect(res_sid2, atom2, vb=vb)
        return self.getAtomDist(vec1, vec2)
    
    
    def getDomainMinDist(self, res_sids1, res_sids2, atom="CA"):
        min_dr = 999999.0
        min_rs1 = -1
        min_rs2 = -1 
        for res_sid1 in res_sids1:
            vec1 = self.getResAtomVect(res_sid1, atom)
            if vec1 is None: continue 
            for res_sid2 in res_sids2:
                vec2 = self.getResAtomVect(res_sid2, atom)
                if vec2 is None: continue
                dr = self.getAtomDist(vec1, vec2)
                if dr < min_dr:
                    min_dr = dr 
                    min_rs1 = res_sid1
                    min_rs2 = res_sid2 
        return min_dr, min_rs1, min_rs2
    
                
    def getAtomsCenter(self, atoms, vb=0):
        num_atom = len(atoms)
        xyz= _array([0.,0.,0.])
        for resid, atom in atoms:
            res_sid = self.resid2ResSid(resid)
            if vb: print("# resid= %4d res_sid=%4d  atom=%4s" % (resid, res_sid, atom))
            xyz += self.getResAtomVect(res_sid, atom)
        xyz = xyz/num_atom
        return xyz 
        
    def getAtomsCenterDist(self, atoms1, atoms2, vb=0):
        ct1 = self.getAtomsCenter(atoms1, vb=vb)
        ct2 = self.getAtomsCenter(atoms2, vb=vb)
        return self.getAtomDist(ct1, ct2)
    
    
    def getDihe(self, atoms, vb=0):
        atom1, atom2, atom3, atom4 = atoms 
        vec1 = self.getResAtomVect(self.resid2ResSid(atom1[0]), atom1[1], vb=vb)
        vec2 = self.getResAtomVect(self.resid2ResSid(atom2[0]), atom2[1], vb=vb)
        vec3 = self.getResAtomVect(self.resid2ResSid(atom3[0]), atom3[1], vb=vb)
        vec4 = self.getResAtomVect(self.resid2ResSid(atom4[0]), atom4[1], vb=vb)
        dihe  = dihedral_vect(vec1, vec2, vec3, vec4,deg=True)
        if vb:
            print("# atoms: ", atoms) 
            print("# dihe: %7.2f" % dihe)
        return dihe
     
   
    def isInSameChain(self, res_sids):
        chain_id = self.resSid2Res(res_sids[0]).getChainID()
        for res_sid in res_sids[1:]:
            _cid = self.resSid2Res(res_sid).getChainID()
            if _cid != chain_id: return False
        else:
            return True
        
    
    def atomsXYZDt2Res(self, atoms_xyz_dt, resid, resn, atomid=1):
        _res = None
        for atom_name in list(atoms_xyz_dt.keys()):
            xyz = atoms_xyz_dt[atom_name]
            #print "# new atomid: ", atomid, atom_name
            atom = ATOM(atomname=atom_name, atomid=atomid, xyz=xyz, resid=resid,  resname=resn,  misc="", segid="", chain_id="")
            if _res is None:
                _res = RESIDUE(atom)
                #print "new_res: ",  _res.atom_names
            else:
                _res.addAtom(atom)
            atomid += 1
        #print "# res", _res.atom_names
        self.addResidue(_res)    
        
    
    def prnResSidsAndResids(self, res_sids, cmt=""):
        print(cmt)
        print("# res_sids:", res_sids)
        print("# resids: ",  self.segResSid2Resid(res_sids))
        
        
    
    def segResSid2Resid(self, res_sids):
        resids = []
        for res_sid in res_sids:
            resids.append(self.resSid2Resid(res_sid))
        return resids

    def getPDBResSids(self):
        return self.res_sids
    
    def getPDBResids(self, vb=0):
        pdb_resids = []
        for res_sid in self.res_sids:
            resid = self.resSid2Resid(res_sid)
            if resid in pdb_resids:
                print("# Waring  resid alreaid existed in the PDB file")
                print("# res_sid= %4d  resid= %4d  " % (res_sid, resid))
            else:
                pdb_resids.append(resid)
        if vb: 
            print("# pdb resids: ", pdb_resids)
        return pdb_resids
    
    def getResCaDistSquare(self, res1, res2):
        vec1 = res1.getAtomVect("CA")
        vec2 = res2.getAtomVect("CA")
        return self.getDistSquare(vec1, vec2)
    
    def getResCentDist(self, resid1, resid2, vb=0):
        ca1 = self.res_centers[resid1]
        ca2 = self.res_centers[resid2]
        dr2 = self.getDistSquare(ca1, ca2)
        dr = sqrt(dr2)
        if vb: 
            print("# ca1: ", ca1)
            print("# ca2: ", ca2) 
            print("# rescenter dist: %4d -- %4d : dr= %8.4f  dr2= %8.4f " % (resid1, resid2, dr, dr2))
        return dr 
    
    def resid2name(self, resid, chain_id):
        _res = self.resid2Res(resid, chain_id=chain_id)
        return _res.getName()
    
    def ressid2name(self, ressid):
        _res = self.resSid2Res(ressid)
        return _res.getName()
    
    ### New Functs 11/07/2014
    def getSegSeq(self, seg, chain_id=None, vb=0):
        if chain_id is None:
            chain_id = self.chains[0]
        seg_seq = []
        for _resid in range(seg[0], seg[1]+1):
            seg_seq.append(self.resid2name(_resid, chain_id))
        if vb: 
            print("# seg: %4d -- %4d" % (seg[0], seg[1]))
            print(seg_seq)               
        return seg_seq
    
    def getSegSids(self, seg):
        """ segment in the same chain"""
        seg_start, seg_end = seg 
        chain_id = seg_start[0]
        assert chain_id == seg_end[0]
        seg_start_resid = seg_start[1]
        seg_start_resn = seg_start[2]
        seg_end_resid = seg_end[1]
        seg_end_resn = seg_end[2]
        seg_sids = []
        for resid in range(seg_start_resid, seg_end_resid+1):
            res_sid = self.resid2ResSid(resid, chain_id=chain_id)
            seg_sids.append(res_sid)
        return seg_sids        
        
    
    def getChainSeq(self, chain_id=None, vb=0, INFO=False):
        if chain_id is None:
            return "# Error: no chain_id input"
        _chain_ressids = self.chain_res_sids[chain_id]
        chain_seq = []
        cnt = 0 
        for _ressid in _chain_ressids:
            cnt += 1
            #print "%4d: %4s   %4d - %4d"  % (cnt, chain_id, _ressid, self.resSid2Resid(_ressid))
            chain_seq.append(self.ressid2name(_ressid))
        num_seq = len(chain_seq)
        num_res = len(_chain_ressids)
        if vb:  print("# chain: %4s num_res: %4d  seq_len: %4d" % (chain_id, num_res, num_seq))
        if INFO:
            chain_first_resid = self.resSid2Resid(_chain_ressids[0])
            chain_last_resid = self.resSid2Resid(_chain_ressids[-1])
            return chain_seq, chain_first_resid, chain_last_resid, num_res 
        else:
            return chain_seq
        
    
    
    def getPDBSeq(self):
        """ return pdb sequence in full amino acid name format:
        eg. ['Pro', ...
        """
        pdb_seq = []
        for _ressid in self.res_sids:
            pdb_seq.append(self.ressid2name(_ressid))
        return pdb_seq
                        
    
    def getPDBSeqLabel(self,  chain_id=None):
        seq_dt =  {}
        if chain_id is None:
            for res_sid in self.res_sids:
                _res =  self.residues[res_sid]
                seq_dt[res_sid] =  _res.name,  _res.resid, _res.chain_id
        else:
            for res_sid in self.res_sids:
                _res =  self.residues[res_sid]
                if chain_id == _res.chain_id:
                    seq_dt[res_sid] =  _res.name,  _res.resid, _res.chain_id
        return seq_dt
    
    
    def getChainResidRange(self, chain_id):
        _chain_ressids = self.chain_res_sids[chain_id]
        chain_start_resid = self.resSid2Resid(_chain_ressids[0])
        chain_end_resid = self.resSid2Resid(_chain_ressids[-1])
        num_res = len(_chain_ressids)
        return num_res, chain_start_resid, chain_end_resid
    
    
    def getChainResSids(self, chain_id):
        return self.chain_res_sids[chain_id]
    
    
    
    def seq2resSid(self, seq, chain_id=None):
        """ match sequence(fasta) to residues, no more or less residue in sequence"""
        num_seq = len(seq)
        i = 0 
        j = 0 
        matched = {}
        excluded = ["-", "*"]
        while i < num_seq: 
            if seq[i] in excluded:
                i+=1 
                continue
            res_sid = self.res_sids[j]
            res = self.resSid2Res(res_sid)
            if chain_id is not None and res.chain_id != chain_id:
                j+= 1
                continue
            nt = self.resName2Letter(res.name)
            if nt!= seq[i]: 
                j+=1 
                continue
            else:
                matched[i] = res_sid, nt, res.name, res.resid
                i+= 1
                j+=1
        return matched
    
                
    ### possible deprecated functions
    def setUniqResids(self, bs_resid=None, vb=0):
        """ adjust residue resids in chains so that each residue has a uniq resid
            bs_resid: chain b start resid 
        """
        #print "# test"
        if hasattr(self, 'UNIQRES') and self.UNIQRES:
            print("# uniqres, nothing done")
            return 
        else:
            self.UNIQRES = 1
        if len(self.chain_ids) == 1:
            resids = list(self.residues.keys())
            resids.sort()
            self.resids = resids 
        elif len(self.chain_ids)==2:
            chain_A_residues = self.all_residues[self.chain_ids[0]]
            chain_A_resids = list(chain_A_residues.keys())
            chain_A_resids.sort()
            chain_B_residues = self.all_residues[self.chain_ids[1]]
            chain_B_resids = list(chain_B_residues.keys())
            chain_B_resids.sort()
            residues = chain_A_residues
            last_a_resid = chain_A_resids[-1]
            if bs_resid is None:
                b_resid = int(ceil(last_a_resid/100.0)*100.0)
                #if (chain_B_resids[0] - last_a_resid) < 5:
                    
                    #if chain_B_resids[0] - b_resid < 5:
                    #    b_resid += 100 
                b_start_resid = b_resid + 1
                    #self.b_start_resid = b_start_resid
            else:
                b_start_resid = bs_resid
            print("# Chain B: start_resid: ", b_start_resid)    
            self.b_start_resid = b_start_resid 
            for resid in chain_B_resids:
                residues[b_start_resid] = chain_B_residues[resid]
                b_start_resid += 1
            self.residues = residues
            resids = list(residues.keys())
            resids.sort()
            self.resids = resids
        else:
            print("# --------------------------------------")
            print("# only allow two chains in the pdb file.")
            print("#_______________________________________")
            sys.exit(0)
        if vb:
            print("# num_res= %4d " % len(self.resids))
            print("# resids: ", self.resids) 
            

    def getSegXYZ_Wo_H(self, seg, list=1, names=None, BACKBONE=False, SIDECHAIN=False):
        ''' seg includes only the res_sids, not resids to accomdate multiple chains with same residue numbers
        '''        
        #print self.residues.keys()
        seg_xyz = []
        for res_sid in seg:
            _seg_xyz = self.residues[res_sid].getResXYZ_Wo_H(list=1, names=names, BACKBONE=BACKBONE, SIDECHAIN=SIDECHAIN)
            seg_xyz += _seg_xyz
        if list:
            return seg_xyz
        else:
            return _array(seg_xyz)
    
    
    
    def getSegResXYZ(self, seg, list=1, fmt="a", resid_type="resid"):
        #print self.residues.keys()
        if resid_type ==  "resid":
            _id2res =  self.resid2Res
        else:
            _id2res =  self.resSid2Res
        if fmt=="d":
            seg_xyz = {}
            for resid in seg:
                seg_xyz[resid] = _id2res(resid).getResXYZ(list=1)            
            return seg_xyz
        elif fmt == 'a':
            seg_xyz = []
            for resid in seg:
                _seg_xyz = _id2res(resid).getResXYZ(list=1)
                seg_xyz += _seg_xyz
            return _array(seg_xyz)        
        else:
            raise "unknow format for atom xyz"
    
    
    def getSegResXYZ_WoH(self,  seg,  list=1,  fmt="a"):
        if fmt=="d":
            seg_xyz = {}
            for resid in seg:
                seg_xyz[resid] = self.resid2Res(resid).getResXYZ_Wo_H(list=1)            
            return seg_xyz
        elif fmt == 'a':
            seg_xyz = []
            for resid in seg:
                _seg_xyz = self.resid2Res(resid).getResXYZ_Wo_H(list=1)
                seg_xyz += _seg_xyz
            return _array(seg_xyz)        
        else:
            raise "unknow format for atom xyz"
    
    
    def getSegAtomVecs(self, atom_name, chain_id, start_resid, end_resid):
        seg_resids = list(range(start_resid , end_resid+1))
        seg_xyz =[]
        for resid in seg_resids:
            seg_xyz.append(self.resid2Res(resid, chain_id=chain_id).getAtomVect(atom_name))
        return seg_xyz    
    
    def getResAtomsVecs(self, ressid, atom_names):
        _res = self.resSid2Res(ressid)
        vecs = []
        for _atom_name in atom_names:
            vecs.append(_res.getAtomVect(_atom_name))
        return _array(vecs)
    
                      
    #def getResAtomVect(self, resid, atom_name, chain_id=None): 
    #    _res = self.resid2Res(resid, chain_id=chain_id)
    #    return _res.getAtomVect(atom_name)
    
    
    def getSegSidAtomVecs(self, ressids, atom_name):
        vecs = []
        for _ressid in ressids:
            _res = self.resSid2Res(_ressid)
            vecs.append(_res.getAtomVect(atom_name))
        return _array(vecs)
            
            
    def getSegResXYZ_WoH_atom(self, seg, list=1):
        seg_res_xyz = {}
        for res_sid in seg:
            seg_res_xyz[res_sid] = self.residues[res_sid].getResXYZ_Wo_H(list=1)
            #seg_res_xyz[resid] = self.residues[resid].getResXYZ_WoH_atom()
        return seg_res_xyz


    def getSegResXYZ_WoH_dict(self, seg):
        seg_res_xyz = {}
        for resid in seg:
            #seg_res_xyz[resid] = self.residues[resid].getResXYZ_Wo_H(list=1)
            _res =  self.resid2Res(resid)
            seg_res_xyz[resid] = _res.getResXYZ_Wo_H()
        return seg_res_xyz            

    
    def getDmax(self, num_block=20, vb=0):
        """ Calcuate Dmax by looping all atoms in molecule """
        dr2_max = 0
        _getDistSq = self.getDistSquare
        resids = self.res_sids
        #num_res = len(resids) + 1
        #resids = range(7, num_res)
        #atom_XYZs = self.getSegResXYZ_WoH_atom(resids)
        atom_XYZs = self.getSegResXYZ(resids, fmt="d")
        #print "atom_xyzs: ", atom_XYZs
        num_res = len(resids)
        dr2_max = 0.0
        dr_cut = 0.0
        for i in range(num_res-num_block-1):
            resid1 = resids[i]
            res1_xyz = atom_XYZs[resid1]
            for j in range(i+num_block, num_res):
                resid2 = resids[j]
                res2_xyz = atom_XYZs[resid2]
                #print "# res1: xyz: ", res1_xyz
                #print "# res2: xyz: ", res2_xyz
                dr_res = sqrt(_getDistSq(res1_xyz[0], res2_xyz[1]))
                if dr_res < dr_cut: 
                    #print "# resid1= %4d resid2= %4d  dr_res= %7.2f  dr_cut= %7.2f" % (resid1, resid2, dr_res, dr_cut)
                    continue
                for v1 in res1_xyz:
                    for v2 in res2_xyz:
                        dr2_12 = _getDistSq(v1, v2)
                        if dr2_12 > dr2_max:
                            dr2_max = dr2_12
                            dr_max = sqrt(dr2_max)
                            dr_cut = dr_max - 20.0
                            #print "#new: resid1= %4d  resid2=%4d dr= %7.2f dr_cut= %7.2f" % (resid1, resid2, dr2_max, dr_cut)
        Rmax = sqrt(dr2_max)
        if vb: 
            print("# num_res= %4d " % num_res)
            print("# full search(over all atoms): Dmax= %8.4f" % Rmax)
        return Rmax             
    
    
    def cGetDmax(self, vb=0, **rest): # pass
        """ get Dmax by using C functions: c_getDmax() """
        resids = self.resids
        #print "# **** ressids: ",  ressids
        atom_XYZs = self.getSegResXYZ(resids, list=0, **rest)        
        dmax = c_getDmax(atom_XYZs, len(atom_XYZs))
        if vb: print("# cget_Dmax= %8.2f " % dmax) 
        return  dmax     
  
    
    def delRes(self, res_sid):
        del self.residues[res_sid]
        if res_sid in self.res_sids:
            self.res_sids.remove(res_sid)
        self.num_res -= 1
        
    
  
 # HELIX parsing
    def parseHelix(self):
        """ HELIX records are used to identify the position of helices in the molecule. Helices are named and numbered. The residues where the helix begins and ends are noted, as well as the total length.
Record Format


Details

Additional HELIX records with different serial numbers and identifiers occur if more than one helix is present.
The initial residue is the N-terminal residue of the helix.
Helices are classified as follows:

                TYPE OF HELIX          CLASS NUMBER 
                                       (COLUMNS 39 - 40)
      ---------------------------------------------------
      Right-handed alpha (default)       1
      Right-handed omega                 2
      Right-handed pi                    3
      Right-handed gamma                 4
      Right-handed 310                   5
      Left-handed alpha                  6
      Left-handed omega                  7
      Left-handed gamma                  8
      27 ribbon/helix                    9
      Polyproline                       10

Relationships to Other Record Types

There may be related information in the REMARKs.
Example

         1         2         3         4         5         6         7
1234567890123456789012345678901234567890123456789012345678901234567890123456
HELIX    1  HA GLY A   86  GLY A   94  1                                   9
HELIX    2  HB GLY B   86  GLY B   94  1                                   9

"""
        for ln in self.pdb_head_lns:
            if ln.startswith("HELIX"):
                self.helix2attr(ln)
                #print ln 
                
    def helix2attr(self, ln):
        """
        COLUMNS       DATA TYPE        FIELD        DEFINITION
        --------------------------------------------------------------------
         1 -  6       Record name      "HELIX "
         8 - 10       Integer          serNum       Serial number of the helix.
                                                    This starts at 1 and increases
                                                    incrementally.
        12 - 14       LString(3)       helixID      Helix identifier. In addition
                                                    to a serial number, each helix is
                                                    given an alphanumeric character
                                                    helix identifier.
        16 - 18       Residue name     initResName  Name of the initial residue.
        20            Character        initChainID  Chain identifier for the chain
                                                    containing this helix.
        22 - 25       Integer          initSeqNum   Sequence number of the initial
                                                    residue.
        26            AChar            initICode    Insertion code of the initial
                                                    residue.
        28 - 30       Residue name     endResName   Name of the terminal residue of
                                                    the helix.
        32            Character        endChainID   Chain identifier for the chain
                                                    containing this helix.
        34 - 37       Integer          endSeqNum    Sequence number of the terminal
                                                    residue.
        38            AChar            endICode     Insertion code of the terminal
                                                    residue.
        39 - 40       Integer          helixClass   Helix class (see below).
        41 - 70       String           comment      Comment about this helix.
        72 - 76       Integer          length       Length of this helix.
        """
        ser_num = int(ln[7:10])
        helix_id =  ln[11:14].strip()
        start_resname = ln[15:18]
        start_chain_id = ln[19]
        start_resid = int(ln[21:25])
        end_resname = ln[27:30]
        end_chain_id = ln[31]
        end_resid = int(ln[33:37])
        helix_class = int(ln[38:40])
        helix_lengh = int(ln[71:76])
        if ser_num in self.helices:
            print("#***  Error in parse pdb HELIX: duplicate serNum. ")
            
        else:
            self.helices[ser_num] = {"id":helix_id,
                                     "start_resname":start_resname,
                                     "start_chain_id":start_chain_id,
                                     "start_resid":start_resid,
                                     "end_resname":end_resname,
                                     "end_chain_id":end_chain_id,
                                     "end_resid":end_resid,
                                     "class":helix_class,
                                     "length":helix_lengh
                                     }
            self.helix_sids.append(ser_num)
            
    
    def getChainHelix(self, chain_id):
        chain_helices = {}
        chain_hids = []
        for helix_sid in self.helix_sids:
            _helix = self.helices[helix_sid]
            if _helix["start_chain_id"] == chain_id:
                #print _helix
                chain_helices[helix_sid] = _helix
                chain_hids.append(helix_sid)
        return chain_helices, chain_hids
    
    def parseSheet(self):
        """
        SHEET
        
        Overview
        
        SHEET records are used to identify the position of sheets in the molecule. Sheets are both named and numbered. The residues where the sheet begins and ends are noted.
        
        Record Format
        
        COLUMNS        DATA TYPE       FIELD           DEFINITION
        ----------------------------------------------------------------------------------
         1 -  6        Record name     "SHEET "
        
         8 - 10        Integer         strand          Strand number which starts at 1 for
                                                       each strand within a sheet and
                                                       increases by one.
        
        12 - 14        LString(3)      sheetID         Sheet identifier.
        
        15 - 16        Integer         numStrands      Number of strands in sheet.
        
        18 - 20        Residue name    initResName     Residue name of initial residue.
        
        22             Character       initChainID     Chain identifier of initial residue
                                                       in strand.
        
        23 - 26        Integer         initSeqNum      Sequence number of initial residue
                                                       in strand.
        
        27             AChar           initICode       Insertion code of initial residue
                                                       in strand.
        
        29 - 31        Residue name    endResName      Residue name of terminal residue.
        
        33             Character       endChainID      Chain identifier of terminal
                                                       residue.
        
        34 - 37        Integer         endSeqNum       Sequence number of terminal residue.
        
        38             AChar           endICode        Insertion code of terminal residue.
        
        39 - 40        Integer         sense           Sense of strand with respect to
                                                       previous strand in the sheet. 0
                                                       if first strand, 1 if parallel,
                                                       -1 if anti-parallel.
        
        42 - 45        Atom            curAtom         Registration. Atom name in current
                                                       strand.
        
        46 - 48        Residue name    curResName      Registration. Residue name in
                                                       current strand.
        
        50             Character       curChainId      Registration. Chain identifier in
                                                       current strand.
        
        51 - 54        Integer         curResSeq       Registration. Residue sequence
                                                       number in current strand.
        
        55             AChar           curICode        Registration. Insertion code in
                                                       current strand.
        
        57 - 60        Atom            prevAtom        Registration. Atom name in
                                                       previous strand.
        
        61 - 63        Residue name    prevResName     Registration. Residue name in
                                                       previous strand.
        
        65             Character       prevChainId     Registration. Chain identifier in
                                                       previous strand.
        
        66 - 69        Integer         prevResSeq      Registration. Residue sequence
                                                       number in previous strand.
        
        70             AChar           prevICode       Registration. Insertion code in
                                                       previous strand.
        
        Details
        
        * The initial residue for a strand is its N-terminus. Strand registration information is provided in columns 39 - 70. Strands are listed starting with one edge of the sheet and continuing to the spatially adjacent strand.
        
        * The sense in columns 39 - 40 indicates whether strand n is parallel (sense = 1) or anti-parallel (sense = -1) to strand n-1. Sense is equal to zero (0) for the first strand of a sheet.
        
        * The registration (columns 42 - 70) of strand n to strand n-1 may be specified by one hydrogen bond between each such pair of strands. This is done by providing the hydrogen bonding between the current and previous strands. No registration information should be provided for the first strand.
        
        * For structures which form a closed sheet (beta-barrel), the first strand is repeated as the last strand. An explanatory remark is included in the REMARK section.
        
        * Split strands, or strands with two or more runs of residues from discontinuous parts of the amino acid sequence, are explicitly listed. Provide a description to be included in the REMARK section.
        
        Verification/Validation/Value Authority Control
        
        SHEET records are now being generated automatically by PDB using the Kabsch and Sander algorithm [Kabsch and Sander, Biopolymers 22: 2577-2637 (1983)], although they may be provided by the depositor instead. PDB verifies that named residues exist in the ATOM records.
        
        Relationships to Other Record Types
        
        If the entry contains bifurcated sheets or beta-barrels, the relevant REMARK records must be provided. See the REMARK section for details.
        
        Example
        
                 1         2         3         4         5         6         7
        1234567890123456789012345678901234567890123456789012345678901234567890
        SHEET    1   A 5 THR A 107  ARG A 110  0
        SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
        SHEET    3   A 5 ARG A  87  SER A  91 -1  N  LEU A  89   O  TYR A  97
        SHEET    4   A 5 TRP A  71  ASP A  75 -1  N  ALA A  74   O  ILE A  88
        SHEET    5   A 5 GLY A  52  PHE A  56 -1  N  PHE A  56   O  TRP A  71
        SHEET    1   B 5 THR B 107  ARG B 110  0
        SHEET    2   B 5 ILE B  96  THR B  99 -1  N  LYS B  98   O  THR B 107
        SHEET    3   B 5 ARG B  87  SER B  91 -1  N  LEU B  89   O  TYR B  97
        SHEET    4   B 5 TRP B  71  ASP B  75 -1  N  ALA B  74   O  ILE B  88
        SHEET    5   B 5 GLY B  52  ILE B  55 -1  N  ASP B  54   O  GLU B  73
        
        The sheet presented as BS1 below is an eight-stranded beta-barrel. This is represented by a nine-stranded sheet in which the first and last strands are identical.         
        """
        for ln in self.pdb_head_lns:
            if ln.startswith("SHEET"):
                self.sheet2attr(ln)          
    
    def sheet2attr(self, ln):
        """
            COLUMNS        DATA TYPE       FIELD           DEFINITION
        ----------------------------------------------------------------------------------
         1 -  6        Record name     "SHEET "
         8 - 10        Integer         strand          Strand number which starts at 1 for
                                                       each strand within a sheet and
                                                       increases by one.
        12 - 14        LString(3)      sheetID         Sheet identifier.
        15 - 16        Integer         numStrands      Number of strands in sheet.
        18 - 20        Residue name    initResName     Residue name of initial residue.
        22             Character       initChainID     Chain identifier of initial residue
        23 - 26        Integer         initSeqNum      Sequence number of initial residue
        27             AChar           initICode       Insertion code of initial residue
                                                       in strand.
        29 - 31        Residue name    endResName      Residue name of terminal residue.
        33             Character       endChainID      Chain identifier of terminal
                                                       residue.
        34 - 37        Integer         endSeqNum       Sequence number of terminal residue.
        38             AChar           endICode        Insertion code of terminal residue.
        39 - 40        Integer         sense           Sense of strand with respect to
                                                       previous strand in the sheet. 0
                                                       if first strand, 1 if parallel,
                                                       -1 if anti-parallel.
        42 - 45        Atom            curAtom         Registration. Atom name in current
                                                       strand.
        46 - 48        Residue name    curResName      Registration. Residue name in
                                                       current strand.
        50             Character       curChainId      Registration. Chain identifier in
                                                       current strand.
        51 - 54        Integer         curResSeq       Registration. Residue sequence
                                                       number in current strand.
        55             AChar           curICode        Registration. Insertion code in
                                                       current strand.
        57 - 60        Atom            prevAtom        Registration. Atom name in
                                                       previous strand.
        61 - 63        Residue name    prevResName     Registration. Residue name in
                                                       previous strand.
        65             Character       prevChainId     Registration. Chain identifier in
                                                       previous strand.
        66 - 69        Integer         prevResSeq      Registration. Residue sequence
                                                       number in previous strand.
        70             AChar           prevICode       Registration. Insertion code in
                                                       previous strand.
    Example
           
                    1         2         3         4         5         6         7
           1234567890123456789012345678901234567890123456789012345678901234567890
           SHEET    1   A 5 THR A 107  ARG A 110  0
           SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
           SHEET    3   A 5 ARG A  87  SER A  91 -1  N  LEU A  89   O  TYR A  97
           SHEET    4   A 5 TRP A  71  ASP A  75 -1  N  ALA A  74   O  ILE A  88
           SHEET    5   A 5 GLY A  52  PHE A  56 -1  N  PHE A  56   O  TRP A  71
           SHEET    1   B 5 THR B 107  ARG B 110  0
           SHEET    2   B 5 ILE B  96  THR B  99 -1  N  LYS B  98   O  THR B 107
           SHEET    3   B 5 ARG B  87  SER B  91 -1  N  LEU B  89   O  TYR B  97
           SHEET    4   B 5 TRP B  71  ASP B  75 -1  N  ALA B  74   O  ILE B  88
           SHEET    5   B 5 GLY B  52  ILE B  55 -1  N  ASP B  54   O  GLU B  73                                                       
                                                       
        """
        sheet_strand = int(ln[7:10])
        sheet_id = ln[11:14].strip()
        start_resname = ln[17:20]
        start_chain_id = ln[21]
        start_resid = int(ln[22:26])
        end_resname = ln[28:31]
        end_chain_id = ln[32]
        end_resid = int(ln[33:37])
        #helix_class = int(ln[38:40])
        #helix_lengh = int(ln[71:76])
        #print ln 
        #print "# sheet id: %s" % sheet_id
        if sheet_id in self.sheets:
            if sheet_strand in self.sheets[sheet_id]:
                print("# Error in parse pdb SHEET: duplicate sheet")
                
                
        else:
            self.sheets[sheet_id] = {}
            self.sheet_sids.append(sheet_id)
        self.sheets[sheet_id][sheet_strand] = {"id":sheet_id,
                        "strand":sheet_strand,
                        "start_resname":start_resname,
                        "start_chain_id":start_chain_id,
                        "start_resid":start_resid,
                        "end_resname":end_resname,
                        "end_chain_id":start_chain_id,
                        "end_resid":end_resid,
                        #"class":helix_class,
                        #"length":helix_lengh
                        }
            
    
    def getChainSheet(self, chain_id):
        chain_sheets = {}
        chain_sids = []
        #print self.sheet_sids
        for _sid in self.sheet_sids:
            _sheet = self.sheets[_sid]
            for _strand in list(_sheet.keys()):
                #print '# strand id: ', _strand
                beta = _sheet[_strand]
                #print "# beta: %s " % beta
                #print beta["start_chain_id"]
                #print "###"
                if beta["start_chain_id"] == chain_id:
                    #print "# chain: %s sheet: " % chain_id
                    chain_sheets[_sid] = _sheet
                    chain_sids.append(_sid)
                    break
        #print "# chain_id: %s" % chain_id
        #print "# chain_sheets: %s " %  chain_sheets
        #print "# chain_sids: %s " % chain_sids
        return chain_sheets, chain_sids                
                    
    
class MHCIPDB(PDB):
    def __tmp(self):
        #self.new_chain_ids = []
        #self.ligand_chain_ids = [] # any chain lengh < 30 and > 1 is treated as ligand 
        #self.other_chain_ids = []
        #self.possible_alpha_chain_ids = []
        #self.MIN_NUM_STRAND = 6 
        #self.MIN_HALF_NUM_STRAND = 3 
        pass 
    
    def findLigand(self, vb=0):
        if self.mono_chain_ids is not None:
            chain_ids = self.mono_chain_ids
        else:
            chain_ids = self.chain_ids 
        self.ligand_chain_ids = []
        self.other_chain_ids = []
        self.possible_alpha_chain_ids =[]
        #print("chain_ids: {}".format(chain_ids))
        #print("# mono_chain_ids: {}".format(self.mono_chain_ids))
        for chain_id in chain_ids:
            chain_length = self.getChainLength(chain_id)
            #lns.append( "# chain_id=%4s  chain_length= %4s"  % (chain_id, chain_length))
            if  chain_length > 1:
                if chain_length < 30:
                    #lns.append( "# ligand chain_id: %s" % chain_id)
                    self.ligand_chain_ids.append((chain_id, chain_length))
                elif chain_length < 160:
                    self.other_chain_ids.append((chain_id, chain_length))
                elif chain_length < 369:
                    self.possible_alpha_chain_ids.append((chain_id, chain_length))
                else:
                    self.other_chain_ids.append((chain_id, chain_length))
        self.num_ligand = len(self.ligand_chain_ids)
        if vb: 
            print("# num_ligand= %4d" % self.num_ligand)
            print("# ligand chains: %s" % self.ligand_chain_ids)
    
    
    def mergeHelices(self, helices, hids):
        """ treat helices with a connection less than 10 resiudes as one long helix """
        num_helix = len(hids)
        long_helices = []
        pre_end_resid = -9999
        min_long_helix_length = 5 
        long_helix_start = None
        #for hid in helices: #i in range(num_helix-1):
        for i in range( num_helix):
            hid = hids[i]
            helix = helices[hid]
            helix_start_resid = helix['start_resid']
            helix_end_resid = helix['end_resid']
            helix_chain_id = helix['start_chain_id']
            assert helix_chain_id == helix['end_chain_id']
            #print "#helix: %s %s -- %s " % (hid, helix_start_resid, helix_end_resid)
            resid_gap = helix_start_resid - pre_end_resid
            if resid_gap > 10: # new long helix 
                if long_helix_start is not None:
                    long_helix_length = pre_end_resid - long_helix_start 
                    if long_helix_length > min_long_helix_length:
                        long_helices.append((long_helix_start, pre_end_resid,long_helix_length, pre_chain_id ))
                long_helix_start = helix_start_resid    
            pre_end_resid = helix_end_resid
            pre_chain_id = helix_chain_id
        if long_helix_start is not None:
            long_helix_length = pre_end_resid - long_helix_start 
            if long_helix_length > min_long_helix_length:
                long_helices.append((long_helix_start, pre_end_resid,long_helix_length, pre_chain_id ))        
        return long_helices
        
    def filterLongHelices(self, helices, min_len=30):
        """ filiter out short helices, in MHCI, one helix is 36+-, the other is 43+-"""
        long_helices = []
        for _helix in helices:
            if _helix[2] > min_len:
                long_helices.append(_helix)
        return long_helices 
    
    
    def findMHCIAlphaChain3D(self, vb=0):
        """ find MHC I alpha chain by 2 long  helices and 3/4 beta sheets before the first long helix,
        3/4 beta sheets between first and second long helices """
        #print "# possible alpha chains: ", self.possible_alpha_chain_ids
        self.mhcI_alphas = []
        self.tpt_drs = []
        self.dist_range = 20
        
        for chain_id, chain_length in self.possible_alpha_chain_ids:
            helices, hids = self.getChainHelix(chain_id)
            long_helices = self.mergeHelices(helices, hids)
            if vb: print("# long_helices: {}".format(long_helices))
            long_helices = self.filterLongHelices(long_helices, min_len=20)
            if len(long_helices) > 1: 
                sheets, sheet_ids  = self.getChainSheet(chain_id)
                bind_sheet_ids = self.findBindingSheet(sheets, sheet_ids)
                #print("# bind_sheetids: ", bind_sheet_ids)
                if len(bind_sheet_ids)> 0:
                    mhcI_bind_sheet_ids = self.findMHCIBindSheet(bind_sheet_ids, sheets, long_helices)
                    num_bsheet = len(mhcI_bind_sheet_ids)
                    if num_bsheet == 1:
                        self.mhcI_alphas.append((chain_id, mhcI_bind_sheet_ids[0]))
                        sheet_range_ressids = self.getBetaSheetResSids(sheets[mhcI_bind_sheet_ids[0]])
                        self.getBindingGrooveTemplateDrs(long_helices, sheet_range_ressids)  
                    elif num_bsheet > 1:
                        self.add2log("# Error in find mhcI alpha for chain: %s" % chain_id)
                        self.add2log("# mhcI_bind_sheet: %s" %  mhcI_bind_sheet_ids)     
                    else:
                        pass 
            else:
                if vb: 
                    print("# num_long_helix < 2: %s" % long_helices)
                else:
                    pass 
                
        #return mhcI_alphas
    
    def getBetaSheetResSids(self, sheet):
        range_ressids = []
        for strand_sid in sheet:
            strand = sheet[strand_sid]
            ss_resid = strand["start_resid"]
            se_resid = strand["end_resid"]
            s_chain = strand["start_chain_id"]
            ss_res_sid = self.resid2ResSid(ss_resid, chain_id=s_chain)
            se_res_sid = self.resid2ResSid(se_resid, chain_id=s_chain)
            range_ressids += range(ss_res_sid, se_res_sid)
        return range_ressids
       

    def getSeg2SegMinDr(self, seg_ressids, range_ressids, atom="CA"):
        min_dr = 100
        min_res_sids = -1, -1
        for src_ressid in seg_ressids:
            for res_sid in range_ressids:
                dr = self.getAnyAtomDist(src_ressid, atom, res_sid, atom)
                if dr < min_dr:
                    min_dr = dr 
                    min_res_sids = src_ressid, res_sid
        return min_dr, min_res_sids
  
  
    def getBindingGrooveTemplateDrs(self, long_helices, sheet_range_ressids):
        tpt_drs = []
        #self.getAnyAtomDist(res_sid1, atom1, res_sid2, atom2)
        helix1_sressid, helix1_eressid = self.getHelixTermRessids(long_helices[0])
        helix2_sressid, helix2_eressid  = self.getHelixTermRessids(long_helices[1])
        
        t1_src_ressid = helix1_sressid + 9  
        tgt_ressids = range(helix2_eressid, helix2_eressid - self.dist_range, -1)
        t1_min_dr, t1_min_ressids = self.getSeg2SegMinDr(range(t1_src_ressid, t1_src_ressid +5), tgt_ressids)
        res1 = self.resSid2Res(t1_min_ressids[0])
        res2 = self.resSid2Res(t1_min_ressids[1])
        print("#t1:         res: {:4d}/{}/{}  -- res: {:4d}/{}/{}: dist= {:.4f}".format(res1.resid, res1.name, res1.chain_id, \
                                                                                   res2.resid, res2.name, res2.chain_id, t1_min_dr))
        tpt_drs.append((t1_min_dr, res1.resid, res1.name, res1.chain_id, res2.resid, res2.name, res2.chain_id))
        
        t2_src_ressid = helix1_eressid #- 10
        tgt_ressids = range( helix2_sressid, helix2_sressid+ self.dist_range)
        t2_min_dr, t2_min_ressids = self.getSeg2SegMinDr(range(t2_src_ressid, t2_src_ressid+5),tgt_ressids)
        res1 = self.resSid2Res(t2_min_ressids[0])
        res2 = self.resSid2Res(t2_min_ressids[1])
        #print("#t2: ressid: {:4d}  -- ressid {:4d}: dist= {:.4f}".format(helix1_eressid, t2_min_ressid, t2_min_dr))
        print("#t2:         res: {:4d}/{}/{}  -- res: {:4d}/{}/{}: dist= {:.4f}".format(res1.resid, res1.name, res1.chain_id, \
                                                                                res2.resid, res2.name, res2.chain_id, t2_min_dr))
        tpt_drs.append((t2_min_dr, res1.resid, res1.name, res1.chain_id, res2.resid, res2.name, res2.chain_id))

        
        helix1_mid_ressid = helix1_eressid - 14 #(helix1_sressid + helix1_eressid)/2 
        mid_range = 8 
        tgt_ressids = range(helix2_sressid+mid_range, helix2_eressid-mid_range+1)
        mid_min_dr, mid_min_ressids = self.getSeg2SegMinDr(range(helix1_mid_ressid, helix1_mid_ressid+5), tgt_ressids)
        #print("#mid ressid: {:4d}  -- ressid {:4d}: dist= {:.4f}".format(helix1_mid_ressid, mid_min_ressid, mid_min_dr))
        res1 = self.resSid2Res(mid_min_ressids[0])
        res2 = self.resSid2Res(mid_min_ressids[1])
        #print("#t2: ressid: {:4d}  -- ressid {:4d}: dist= {:.4f}".format(helix1_eressid, t2_min_ressid, t2_min_dr))
        print("#mid:        res: {:4d}/{}/{}  -- res: {:4d}/{}/{}: dist= {:.4f}".format(res1.resid, res1.name, res1.chain_id, \
                                                                                res2.resid, res2.name, res2.chain_id, mid_min_dr))
        tpt_drs.append((mid_min_dr, res1.resid, res1.name, res1.chain_id, res2.resid, res2.name, res2.chain_id))

        
        helix1_mid2sheet_min_dr, h1m_ressids = self.getSeg2SegMinDr(range(helix1_mid_ressid, helix1_mid_ressid+5), sheet_range_ressids)
        #print("#h1 mid2sheet ressid: {:4d}  -- ressid {:4d}: dist= {:.4f}".format(helix1_mid_ressid, h1m_ressid, helix1_mid2sheet_min_dr))
        res1 = self.resSid2Res(h1m_ressids[0])
        res2 = self.resSid2Res(h1m_ressids[1])
        #print("#t2: ressid: {:4d}  -- ressid {:4d}: dist= {:.4f}".format(helix1_eressid, t2_min_ressid, t2_min_dr))
        print("#h1mid2sheet res: {:4d}/{}/{}  -- res: {:4d}/{}/{}: dist= {:.4f}".format(res1.resid, res1.name, res1.chain_id, \
                                                               res2.resid, res2.name, res2.chain_id, helix1_mid2sheet_min_dr))
        tpt_drs.append((helix1_mid2sheet_min_dr, res1.resid, res1.name, res1.chain_id, res2.resid, res2.name, res2.chain_id))

        
        helix2_mid_ressid = (helix2_sressid + helix2_eressid)/2
        helix2_mid2sheet_min_dr, h2m_ressids = self.getSeg2SegMinDr(range(helix2_mid_ressid, helix2_mid_ressid+5), sheet_range_ressids)
        #print("#h2 mid2sheet res: {:4d}  -- res {:4d}: dist= {:.4f}".format(helix2_mid_ressid, h2m_ressid, helix2_mid2sheet_min_dr))
        res1 = self.resSid2Res(h2m_ressids[0])
        res2 = self.resSid2Res(h2m_ressids[1])
        #print("#t2: ressid: {:4d}  -- ressid {:4d}: dist= {:.4f}".format(helix1_eressid, t2_min_ressid, t2_min_dr))
        print("#h2mid2sheet res: {:4d}/{}/{}  -- res: {:4d}/{}/{}: dist= {:.4f}".format(res1.resid, res1.name, res1.chain_id, \
                                                               res2.resid, res2.name, res2.chain_id, helix2_mid2sheet_min_dr))
        tpt_drs.append((helix2_mid2sheet_min_dr, res1.resid, res1.name, res1.chain_id, res2.resid, res2.name, res2.chain_id))
        #print(tpt_drs)
        self.tpt_drs.append(tpt_drs)
    
    def findMHCIBindSheet(self, bind_sheet_ids, sheets, long_helices, vb=0):
        #print bind_sheet_ids
        first_helix_start_resid = long_helices[0][0]
        second_helix_start_resid = long_helices[1][0]
        mhcI_bind_sheet_ids = []
        for sheet_id in bind_sheet_ids:
            sids_before_helix = []
            sids_middle_helix = []            
            sheet = sheets[sheet_id]
            strand_ids = list(sheet.keys())
            for sid in strand_ids:
                strand = sheet[sid]
                strand_start_resid = strand['start_resid']
                if strand_start_resid < first_helix_start_resid:
                    sids_before_helix.append(sid)
                elif strand_start_resid < second_helix_start_resid:
                    sids_middle_helix.append(sid)
                else:
                    pass 
            num_strand_before = len(sids_before_helix)
            num_strand_middle = len(sids_middle_helix)
            if num_strand_before >= self.MIN_HALF_NUM_STRAND and num_strand_middle >= self.MIN_HALF_NUM_STRAND:
                mhcI_bind_sheet_ids.append(sheet_id)
            if vb: print("num_strand: before= %s  middle= %s" % (num_strand_before, num_strand_middle))
        return mhcI_bind_sheet_ids
        
        
        
    def findBindingSheet(self, sheets, sheet_ids):
        bind_sheet_ids = []
        for sheet_id in sheet_ids:
            num_strand_in_sheet = len(sheets[sheet_id])
            #print("# sheet: %s num_strand: %s" % (sheet_id, num_strand_in_sheet))
            if num_strand_in_sheet >= self.MIN_NUM_STRAND:
                bind_sheet_ids.append(sheet_id)
                #print sheets[sheet_id]
        return bind_sheet_ids        
         
            
        
    def getChain2ChainMinDist(self, src_chain_id, dst_chain_id):
        src_chain_res_sids = self.chain_res_sids[src_chain_id]
        dst_chain_res_sids = self.chain_res_sids[dst_chain_id]
        min_dist = 99999
        min_src_res_sid = None
        min_dst_res_sid = None
        for _src_res_sid in src_chain_res_sids:
            for _dst_res_sid in dst_chain_res_sids:
                try:
                    _dist = self.getAnyAtomDist(_src_res_sid, "CA", _dst_res_sid, "CA")
                except TypeError:
                    continue
                if _dist < min_dist:
                    min_dist = _dist
                    min_src_res_sid = _src_res_sid
                    min_dst_res_sid = _dst_res_sid
        ln = "# min_dist: %4s -- %4s : %6.2f"  % ( src_chain_id, dst_chain_id, min_dist)
        print(ln) 
        return min_dist
    
    
    def getligand2ChainMinDist(self, src_chain_id, dst_chain_id):
        src_chain_res_sids = self.chain_res_sids[src_chain_id]
        dst_chain_res_sids = self.chain_res_sids[dst_chain_id]
        min_src_res_sid = None
        min_dst_res_sid = None
        num_res = len(src_chain_res_sids)
        lns = ["# num_res= %4d" % num_res ]
        mid_res_sid = src_chain_res_sids[num_res/2 -1]
        pre_res_sid = src_chain_res_sids[int(num_res * 0.25)]
        tail_res_sid = src_chain_res_sids[int(num_res * 0.75)]
        ligand_res_sids = [ pre_res_sid, mid_res_sid, tail_res_sid] 
        num_ligand_res = len(ligand_res_sids)
        sum_dist = 0.0
        sum_cnt = 0
        for _src_res_sid in ligand_res_sids:
            min_dist = 99999
            for _dst_res_sid in dst_chain_res_sids:
                try:
                    _dist = self.getAnyAtomDist(_src_res_sid, "CA", _dst_res_sid, "CA")
                except TypeError:
                    continue
                if _dist < min_dist:
                    min_dist = _dist
                    #min_src_res_sid = _src_res_sid
                    min_dst_res_sid = _dst_res_sid
            if min_dist < 999:
                _src_res = self.residues[_src_res_sid]
                _dst_res = self.residues[min_dst_res_sid]
                sum_dist += min_dist
                sum_cnt += 1
                lns.append( "# min_dist: %4s(%4d) -- %4s(%4d) : %6.2f"  % ( src_chain_id, _src_res.resid,  dst_chain_id, _dst_res.resid, min_dist) )
        min_dist = sum_dist/float(sum_cnt)
        self.prnLines(lns)
        self.add2log(lns)
        return min_dist        


    def getNewChainID(self):
        for cid in self.possible_chain_ids:
            if cid not in self.chain_ids and cid not in self.new_chain_ids:
                return cid 
        raise "# Error: all possible Chain IDs(%s) are used up "  % self.possible_chain_ids

        
    def checkNegResid(self):
        neg_segs = []
        for chain_id in self.chain_ids:
            chain_res_sids = self.chain_res_sids[chain_id]
            neg_res_sids = []
            for _res_sid in chain_res_sids:
                _res = self.residues[_res_sid]
                if _res.resid < 0:
                    neg_res_sids.append(_res_sid)
            num_neg_res = len(neg_res_sids)
            if num_neg_res > 5:
                print("# found a possible segment(%4d)  in the chain: %s" % (num_neg_res, chain_id))
                neg_segs.append((chain_id, num_neg_res, neg_res_sids))
        return neg_segs

    
    def reformNegChain(self, neg_segs):
        for chain_id, num_e, neg_res_sids in neg_segs:
            new_chain_id = self.getNewChainID()
            new_resid = 0
            for _res_sid in neg_res_sids:  
                new_resid += 1
                _res = self.residues[_res_sid]
                _res.setChainID(new_chain_id)
                _res.setResid(new_resid)
                self.chain_res_sids[chain_id].remove(_res_sid)
            self.chain_res_sids[new_chain_id] = neg_res_sids    
            print("# new_chain_id=%4s  num_res= %4d "  % (new_chain_id, new_resid))
            self.chain_ids.append(new_chain_id)
        self.resetChainLengths()    

    
    def getMHCIBindCore(self, vb=0):    
        if vb: 
            print("ligands:" ,  self.ligand_chain_ids)
            print("alphas:", self.mhcI_alphas)
        return self.ligand_chain_ids, self.mhcI_alphas
    
    
    def prnMHCIBindCore(self):
        lns = []
        lns += self.pdb_chain_inf
        lns.append("# ligands: %s" % self.ligand_chain_ids)
        lns.append("# alphas: %s" %  self.mhcI_alphas)
        return lns 
    
    
    def isMHCI3D(self):
        """ check the binding region by two long helices and 6-8 beta strand sheet """
        if len(self.mhcI_alphas)> 0:
            return True
        else:
            return False
        
    def start(self, mono_chain_ids=None, vb=0):
        self.mono_chain_ids = mono_chain_ids
        self.MIN_NUM_STRAND = 5
        self.MIN_HALF_NUM_STRAND = 2 
        
        self.parseHelix()
        self.parseSheet()
        self.pdb_chain_inf = self.prnChainInf(prn=0)
        
        self.findLigand()
        self.findMHCIAlphaChain3D(vb=vb)
        if self.isMHCI3D():
            print('# < %s > is a mhcI pdb' % self.fn_pdb)
            #return self.getMHCIBindCore()
            return True, self.prnMHCIBindCore()
        else:
            #print('# < %s > is not a mhcI pdb' % self.fn_pdb)
            return False, self.pdb_chain_inf
    
    def getMHCIPDBInf(self):
        alpha_chain_id = self.mhcI_alphas[0][0]
        alpha_seq = self.getChainSeq(alpha_chain_id)
        #print("# alpha_chain: {}".format(self.mhcI_alphas))
        #print("# ligands: {}".format(self.ligand_chain_ids))
        return alpha_chain_id, alpha_seq, self.ligand_chain_ids      
         
      

if (__name__ == "__main__"):
    pdb_obj = PDB(fn_pdb='riboA_wH.pdb')
    pdb_obj.start()
