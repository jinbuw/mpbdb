# Mapping protein seuence to HLA allele protein sequence 
 
from lib.libtree import HLATree
from lib.libfasta import FASTA

class PdbSeq2Allele(HLATree, FASTA):
    """ map aligned protein sequence to allele protein sequence """
    def loadMHCIPDBSeq(self):
        """ mhci_pdbs[pdbid] = [(chain_ids, chain_seq(xxx,...), [ligand chains]), ...] """
        fp = self.getMHCIPDBFpBin()
        self.mhcI_pdbs = self.loadObj(fp)
        self.mhcI_pdbids = self.mhcI_pdbs.keys() 
        self.mhcI_pdb_alpha_seqs = {}
        for pdbid in self.mhcI_pdbids:
            alpha_seq = self.mhcI_pdbs[pdbid][1]
            self.mhcI_pdb_alpha_seqs[pdbid] = self.seq2fasta(alpha_seq)
    
        
    def loadMHCIPDB2Allele(self):
        fp_bin  = self.getPDB2ALLeleFpBin()
        if not self.isNew(fp_bin):
            self.matched_pdbs, self.non_matched_pdbs = self.loadObj(fp_bin)
        else:
            self.matched_pdbs = {}
            self.non_matched_pdbs = {}
   
    def saveMHCIPDB2Allele(self):
        fp_bin  = self.getPDB2ALLeleFpBin()
        self.dumpObj((self.matched_pdbs, self.non_matched_pdbs), fp_bin, vb=1)
    
    
    def getPDBAlphaSeqDt(self, pdbid):
        pdb_alpha_seq = self.mhcI_pdb_alpha_seqs[pdbid]
        #aseq_dt = self.alignPDBSeq2RefSeq(pdb_alpha_seq)        
        return self.alignPDBSeq2RefSeq(pdb_alpha_seq)        
        
    
    def mhcIPdb2Allele(self):
        pdbids = self.mhcI_pdbids
        #pdbids = ["1A1M"]
        #pdbids =  ["4QRQ", "1M05", "4U1J", "4U1M", "3FT2"]
        UPDATE = False #True
        self.loadMHCIPDB2Allele()
        cnt = 0 
        pdb_cnt = len(self.matched_pdbs)
        for pdbid in pdbids:
            if not UPDATE:
                if pdbid in self.matched_pdbs or  pdbid in self.non_matched_pdbs: 
                    #print(self.matched_pdbs[pdbid])
                    continue 
            #self.matchProtAseq(aseq_dt)
            #seq_len = len(aseq_dt)
            #self.add2log('# pdb: {}  len: {}'.format(pdbid, seq_len))
            aseq_dt = self.getPDBAlphaSeqDt(pdbid)
            com_aseq_dt = self.getCommRegSeqdt(aseq_dt)
            match_names, mns, mdifs = self.matchProtAseq(com_aseq_dt)
            if match_names:
                self.matched_pdbs[pdbid] = match_names, mns, mdifs
                pdb_cnt += 1
                max_num = min(5, len(match_names))
                ln = "# {:4d} pdb: {} -> {}".format(pdb_cnt, pdbid, match_names[:max_num])
                self.add2log(ln)
                print(ln)
            else:
                self.non_matched_pdbs[pdbid] = mns, mdifs 
            cnt += 1
            #if cnt > 10: break 
        if cnt > 0:
            self.saveMHCIPDB2Allele()
            self.saveLog("mhcIPDB2Allele.log")
        print("# num_pdb2allele: {}".format(len(self.matched_pdbs)))
        
            
    def start(self):
        self.loadHLAAllele() 
        self.createDTree()
        self.loadMHCIPDBSeq() 
        self.setCommonSeqReg((25,206))
        #self.prnComonRefSeq()
        self.mhcIPdb2Allele()
        #self.testMapPdb2Allele()
        
    def testMapPdb2Allele(self):
        check_pdbids = ["4QRQ", "1M05", "4U1J", "4U1M", "3FT2"]
        for pdbid in check_pdbids:
            if pdbid not in self.mhcI_pdbids:
                print("# {} not in MHCI_pdb".format(pdbid))
            #aseq_dt = self.getPDBAlphaSeqDt(pdbid)
            #com_seq_dt = self.getCommRegSeqdt(aseq_dt)
            #seq_dif = self.getSeqDiff("B*08:01", com_seq_dt)
            #print(pdbid, seq_dif)
            #break 
        #return        
        #pdbids = self.mhcI_pdbids
        self.matched_pdbs = {}
        self.non_matched_pdbs = {}
        cnt = 0 
        pdb_cnt = len(self.matched_pdbs)
        for pdbid in check_pdbids: #pdbids:
            if pdbid in self.matched_pdbs or pdbid in self.non_matched_pdbs:
                continue 
            aseq_dt = self.getPDBAlphaSeqDt(pdbid) #self.alignPDBSeq2RefSeq(pdb_alpha_seq)
            #self.matchProtAseq(aseq_dt)
            seq_len = len(aseq_dt)
            self.add2log('# pdb: {}  len: {}'.format(pdbid, seq_len))
            com_aseq_dt = self.getCommRegSeqdt(aseq_dt)
            match_names, mns, mdifs = self.matchProtAseq(com_aseq_dt)
            print(match_names)
            
            #self.chkNodeChildren("B", com_aseq_dt)
            #self.chkNodeAnchor("B*08", com_aseq_dt)
            #self.chkNodeAnchor("B*08:01", com_aseq_dt)
            #self.chkNodeChildren("B*08", com_aseq_dt)
            #print(com_aseq_dt)
            #print(self.root.getMatchChild(com_aseq_dt))
            
            
            
    def testMapSeq2Allele(self):
        #self.testSeq2Allel()
        #self.testGrpMatch("C*14")
        #self.testHLAProtMatch()
        self.testMatchProt("B*27:07")
        #self.testGrpNode("C*14", "C*14:92")
        #self.testNode("C", "C*14:92")
        
    def testNode(self, grp, prot):
        prot_seq = self.getProtSeq(prot)
        node = self.getNode(node_name=grp) #"C*14")
        node.getPossMChildren(prot_seq, vb=1)        
    
        
if  __name__ == "__main__":
    import argparse, sys 
    parser = argparse.ArgumentParser(prog="%s" % __file__, description="cluster pdbs by their sequences")
    parser.add_argument("-v","--version", action='version', version='%(prog)s 1.0')
     
    parser.add_argument("-d", "--data_dir",dest="data_dir", default="example_data",
                        help="input existed  data directory")   
    
    parser.add_argument("-u","--update", dest="UPDATE",  action="store_true",  default=False,
                        help="Update all existing binary files")       
    
    args =parser.parse_args()    
    tobj = PdbSeq2Allele(args)
    tobj.start()
    
      