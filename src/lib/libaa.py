
aa_seq_names = { "R": "Arginine",  #"Arg"        Positive
                "K": "Lysine" , # Lys - K        positive
                "D": "Aspartic acid", # - Asp    Negative 
                "E": "Glutamic acid", #  - Glu   Negative 
                
                #Polar (may participate in hydrogen bonds):
                "Q": "Glutamine",  # Gln
                "N": "Asparagine", # Asn
                "H": "Histidine",  #  His        #positive 
                "S": "Serine",                  # - Ser - S
                "T": "Threonine",               # - Thr - T
                "Y": "Tyrosine" , # - Tyr - Y
                "C": "Cysteine",  #  - Cys - C
                "M": "Methionine", # - Met - M
                "W": "Tryptophan", # - Trp - W

                # Hydrophobic (normally buried inside the protein core):
                "A": "Alanine"    ,              #  - Ala - A
                "I": "Isoleucine",               # - Ile - I
                "L": "Leucine",                  # - Leu - L
                "F": "Phenylalanine" ,           # - Phe - F
                "V": "Valine",  # - Val - V
                "P": "Proline", #  - Pro - P
                "G": "Glycine" # - Gly - G
                }
                
def aaseq2name(seq):
    return aa_seq_names[seq]

def AAHBond(object):
    """ Donors (associated hydrogens are deduced from topology)
        N of the main chain, water OH2/OW, ARG NE/NH1/NH2, ASN ND2, 
        HIS ND1/NE2, SER OG, TYR OH, CYS SG, THR OG1, GLN NE2, LYS NZ, TRP NE1
   
    Acceptors
        O of the main chain, water OH2/OW, ASN OD1, ASP OD1/OD2, CYH SG, GLN OE1,
        GLU OE1/OE2, HIS ND1/NE2, MET SD, SER OG, THR OG1, TYR OH

"""
    
    
    aa_hb_dornors = {"ALA":["N"], "ARG": ["N", "NE", "NH1", "NH2"], "ASN":["N", "ND2"], "ASP":["N"], 
                    "CYS":["N", "SG"], "GLN":["N", "NE2"], "GLU":["N"], "GLY":["N"], 
                    "HIS":["N", "ND1", "NE2"], "ILE":["N"], "LEU":["N"], "LYS":["N", "NZ"],
                    "MET":["N"], "PHE":["N"], "PRO":["N"], "SER":["N", "OG"], 
                    "THR":["N", "OG1"], "TRP":["N", "NE1"], "TYR":["N", "OH"], "VAL":["N"], 
                     }
    aa_hb_acceptors = {"ALA":["O"], "ARG": ["O", "NE", "NH1", "NH2"], "ASN":["O", "OD1"], "ASP":["O"], 
                    "CYH":["O", "SG"], "GLN":["O", "OE1", "OE2"], "GLU":["O"], "GLY":["O"], 
                    "HIS":["O", "ND1", "NE2"], "ILE":["O"], "LEU":["O"], "LYS":["O"],
                    "MET":["O", "SD"], "PHE":["O"], "PRO":["O"], "SER":["O", "OG1"], 
                    "THR":["O"], "TRP":["O"], "TYR":["O", "OH"], "VAL":["O"], 
                     }
   
    
class AACHARGE:
    """ 
        
    letter2name = { 'D':'ASP', 'E':'GLU', 'N':'ASN', 'Q':'GLN',
                    'R':'ARG', 'K':'LYS', 'P':'PRO', 'G':'GLY',
                    'C':'CYS', 'T':'THR', 'S':'SER', 'M':'MET',
                    'W':'TRP', 'F':'PHE', 'Y':'TYR', 'H':'HIS',
                    'A':'ALA', 'V':'VAL', 'L':'LEU', 'I':'ILE'}
    large hydrophobic residue: 
    Phe, Trp, Tyr, Leu, Ile, Val or Met
    
    a short and/or hydrophobic residue at the P6 position (including Ser, Thr, Cys, Ala, Pro, Val, Ile, Leu or Met).
    """ 
    aa_charged =   ['R', 'K', 'D', 'E']
    aa_polar =     ['Q', 'N', 'H', 'S', 'T', 'Y', 'C','M', 'W']
    aa_hydrophic = ['A', 'I', 'L', 'F', 'V', 'P', 'G']
    
    
    charge_types = {"C": "Charged", 'P':'Polar', "H":'Hydrophoic'}

    aa_positive = ["R", "K"]
    aa_negative = ["D", "E"]
    aa_ext_negative = aa_negative +  ["C", "F", "H", "W", "Y" ]
    
    def catAACharge(self, aa):
        if aa in self.aa_charged:
            return "C"
        elif aa in self.aa_polar:
            return "P"
        elif aa in self.aa_hydrophic:
            return "H"
        else:
            print("# Error: unknown charge in 20 AAs")
    
    
    def resNT2HPC(self, res_nt):
        """ return residue's Hydrophobicity, polarity and charge"""
        return self.catAACharge(res_nt)
        
    
    def resNT2PN(self, res_nt):
        """ return positive charge or negative charge or neutral """
        if res_nt in self.aa_positive:
            return 1
        elif res_nt in self.aa_negative:
            return -1
        else:
            return 0 
        
    def resNT2EPN(self, res_nt):
        if res_nt in self.aa_positive:
            return 1
        elif res_nt in self.aa_ext_negative:
            return -1
        else:
            return 0 
  
    def resNT2PNCP(self, res_nt):
        if res_nt in self.aa_positive:
            return "+"
        elif res_nt in self.aa_ext_negative:
            return "-"
        else:
            return "O" 

        
        
    def grpAAbyCharge(self, aas):
        charges = {}
        for aa in aas:
            _charge_type = self.catAACharge(aa)
            if _charge_type in charges:
                charges[_charge_type].append(aa)
            else:
                charges[_charge_type] = [aa]
        return charges
    
    def prnGrpAACharge(self, aa_grp):
        ln = ""
        num_aa = 0
        for _ct in ["C", "P", "H"]:
            if _ct in aa_grp:
                ln += " %s:[%s]" % (self.charge_types[_ct], ",".join(aa_grp[_ct]))
                num_aa += len(aa_grp[_ct])
        ln = "# num_aa= %d %s"  % (num_aa, ln)
        print(ln) 
    
    def grpReduLigand(self, aa_pdbids):
        aa_grp_pdbids = {}
        for pdbid in aa_pdbids:
            grpid = self.pdbid_grp[pdbid]
            if grpid in aa_grp_pdbids:
                aa_grp_pdbids[grpid].append(pdbid)
            else:
                aa_grp_pdbids[grpid] = [pdbid]
        return aa_grp_pdbids
    