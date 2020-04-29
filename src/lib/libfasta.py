# transform one letter nucleotide name to 3-letter Xplor format 
# deprecated module !!!

from .libio import MIO

class FASTA(MIO):
    oneletter = {
            'ASP':'D','GLU':'E','ASN':'N','GLN':'Q',
            'ARG':'R','LYS':'K','PRO':'P','GLY':'G',
            'CYS':'C','THR':'T','SER':'S','MET':'M',
            'TRP':'W','PHE':'F','TYR':'Y','HIS':'H',
            'ALA':'A','VAL':'V','LEU':'L','ILE':'I',
        }        
    
    nt_names = {"U":"Uri", "G":"Gua", "A":"Ade", "C":"Cyt"} 
    nt_Cnames = {"U":"URI", "G":"GUA", "A":"ADE", "C":"CYT", "ADE":"ADE"}  
    def seq2fasta(self, seq, vb=0):
        seq_fasta = ""
        for _name in seq:
            seq_fasta += self.oneletter[_name]
        if vb:
            print("seq_len=%4d  fasta_len=%4d" % (len(seq), len(seq_fasta)))
            print("# fasta_seq:%s:" % seq_fasta)
        return seq_fasta
        
    def resn2fnt(self, resn):
        """ convert resn(3 letter) to 1 letter """
        return self.oneletter[resn]
    
     
    
    def readPDBFasta(self, fp):
        ccs = self.readTxtFile(fp)
        seq = ""
        pdb_seqs = {}
        pdbid = "NA"
        chain_id = "NA"
        for ln in ccs:
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue 
            if ln.startswith(">") or ln.startswith("<"): # new sequence #1AQD:G|PDBID|CHAIN|SEQUENCE
                if seq:
                    if pdbid not in pdb_seqs:
                        pdb_seqs[pdbid] = {}
                    pdb_seqs[pdbid][chain_id] = seq     
                    seq = ""
                tmp = ln[1:].split("|")
                pdb_code = tmp[0]
                if ":" in pdb_code:
                    pdbid, chain_id = pdb_code.split(":")
                else:
                    pdbid = pdb_code
                    chain_id = "NA"
            else:
                seq += ln 
        #print "# pdbid=%4s chain_id=%2s"  %(pdbid, chain_id)        
        #print len(ccs)
        if seq:
            #print pdbid, chain_id, seq 
            if pdbid not in pdb_seqs:
                pdb_seqs[pdbid] = {}            
            pdb_seqs[pdbid][chain_id] = seq 
        return pdb_seqs            
    
    def write_PDB_fasta(self, pdb_seqs, fn_fasta, pdbids=None, vb=0):
        """ Write PDB  chain sequence(s) into a fasta file """
        if pdbids is None:
            pdbids = list(seqs.keys())
        pdbids.sort()
        lns = []
        for pdbid in pdbids:
            chain_ids = list(pdb_seqs[pdbid].keys())
            chain_ids.sort()
            for chain_id in chain_ids:
                lns.append(">{}:{}|PDBID|CHAIN:SEQUENCE".format(pdbid, chain_id))
                lns.append(pdb_seqs[pdbid][chain_id])
        self.saveLines(lns, fn_fasta, vb=vb)
        
        
    def readAlleleFasta(self, fp):
        """ read allele's fasta file and parse allele name, sequence length:
        >HLA:HLA08460 A*01:118 181 bp
        SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWD
        QETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGK
        DYIALNEDLRSWTAADMAAQITKRKWEAVHAEEQRRVYLEGRCVDGLRRYLENGKETLQR
        T
        """
        ccs = self.readTxtFile(fp)
        self.allele_seqs =  {}
        self.allele_name_seqlen = {}
        seq = ""
        for ln in ccs:
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue 
            if ln.startswith(">"):
                if seq:
                    self.__addAllele(allele_name, seq_len, seq)
                    seq = ""
                elems = ln[1:].split()
                allele_name = elems[1]
                seq_len = int(elems[2])
            else:
                seq += ln 
        if seq: self.__addAllele(allele_name, seq_len, seq)
        return self.allele_seqs, self.allele_name_seqlen
       
    def __addAllele(self, allele_name, seq_len, seq):
        assert allele_name not in self.allele_seqs
        self.allele_seqs[allele_name] = seq     
        assert len(seq) == seq_len
        self.allele_name_seqlen[allele_name] = seq_len
    
    
    def readPdbSeqs(self, fn_inp):
        print("# fasta file for all downloaded pdbs: {}".format(fn_inp))
        self.pdb_seqs = self.readPDBFasta(fn_inp)
        self.pdbids = list(self.pdb_seqs.keys())
        self.num_pdb = len(self.pdbids)
        print("# num_pdb= {:4d} in file:{}".format(self.num_pdb, fn_inp))
            
        
class FilterMHCbySeq(FASTA):
    """ filter MHC class II pdbs by their sequences """
    def rmMultimers(self):
        #pdb_rm_chains = {}
        self.pdb_mono_chains = {}
        for pdbid in self.pdbids:
            _pdbseq = self.pdb_seqs[pdbid]
            chain_ids = list(_pdbseq.keys())
            chain_ids.sort()
            _rm_chains = []
            existed_seqs = []
            chain_inf = []
            mono_chains = []
            for chain_id in chain_ids:
                _seq = _pdbseq[chain_id]
                chain_len = len(_seq)
                if _seq in existed_seqs:
                    _rm_chains.append(chain_id)
                else:
                    existed_seqs.append(_seq)
                    mono_chains.append(chain_id)
                chain_inf.append((chain_id, chain_len))
            self.add2log("# {} : mono {}  {}".format(pdbid, mono_chains, chain_inf))
            self.add2log("# pdbid: {}  rm_chains: {}".format(pdbid, _rm_chains))
            self.pdb_mono_chains[pdbid] = mono_chains
                
    def filterPDBbySeqLen(self, min_chain_len=160, min_pep_len= 3,  max_pep_len=30):
        """ filter out longest chain lenth < 160, maximum binding peptide length= 30  """
        not_valid_pdbs = []
        no_ligand_pdbs = []
        valid_pdbs = []
        num_not_valid = 0
       
        stat_ligand = {}
        for pdbid in self.pdbids:
            num_ligand  = 0
            num_valid_chain = 0
            _pdbseq = self.pdb_seqs[pdbid]
            ln = "%s: " % pdbid
            LIGSTAT = True
            for chain_id in self.pdb_mono_chains[pdbid]: 
                chain_len = len(_pdbseq[chain_id])
                if chain_len > min_chain_len:
                    num_valid_chain += 1
                elif chain_len < max_pep_len and chain_len > min_pep_len:
                    num_ligand += 1
                    if LIGSTAT:
                        stat_ligand.setdefault(chain_len,[]).append(pdbid)
                        LIGSTAT = False
                ln += "%s:%-5s"  % (chain_id, chain_len)
            self.add2log(ln)
            if num_valid_chain  < 1:
                self.add2log("# pdb: %4s  is not valide MHCI " % pdbid)
                num_not_valid += 1
                not_valid_pdbs.append((pdbid, num_valid_chain))
            elif num_ligand == 0:
                no_ligand_pdbs.append(pdbid)
            else:
                valid_pdbs.append(pdbid)
        self.saveFilterPDB(valid_pdbs, no_ligand_pdbs, not_valid_pdbs,stat_ligand)        
    
    
    def prnFilteredbySeq(self):
        self.num_no_ligand = len(self.no_ligand_pdbs)
        self.num_valid_pdb = len(self.valid_pdbs)
        self.num_not_valid = len(self.not_valid_pdbs)        
        ln =  "# num_pdb=%4d  num_not_valid_pdb = %4d num_no_ligand=%4d num_ligand= %4d " % \
              (self.num_pdb, self.num_not_valid, self.num_no_ligand, self.num_valid_pdb)
        print(ln)
        self.add2log(ln)
        
                
    def saveFilterPDB(self,valid_pdbs, no_ligand_pdbs, not_valid_pdbs, stat_ligand):
        fp_bin = self.getSeqFilteredFpBin()
        self.dumpObj((self.num_pdb, valid_pdbs, no_ligand_pdbs, not_valid_pdbs, self.pdb_mono_chains), fp_bin, vb=1)
        self.valid_pdbs = valid_pdbs
        self.no_ligand_pdbs = no_ligand_pdbs
        self.not_valid_pdbs = not_valid_pdbs
        #self.prnFilteredbySeq()
        lig_lengths = list(stat_ligand.keys())
        lig_lengths.sort()
        for lig_len in lig_lengths:
            self.add2log("# ligand len: %2d  num_pdb: %4d" % (lig_len, len(stat_ligand[lig_len])))
        self.saveLines(self.log, self.fp2FpLog(fp_bin), vb=1)
                
    