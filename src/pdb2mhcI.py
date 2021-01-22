# find mhc II in pdb file by the geometry of MHC class II (two helices and 8 beta sheets)
""" select MHCI protein 3D structures with the peptide binding groove """

from lib.libpdb import MHCIPDB
from lib.libmhcidb import MHCIDBData
from lib.libfasta import FilterMHCbySeq

class PDB2MHCI(MHCIDBData, FilterMHCbySeq):
    """ find MHC  MHC class II  in pdb file """
     
    
    def loadFilteredPDBbySeq(self):
        fp = self.getSeqFilteredFpBin()
        if self.isNew(fp):
            fp_fasta = self.getDLFastaFp()
            self.readPdbSeqs(fp_fasta)
            self.rmMultimers()
            self.filterPDBbySeqLen()
        else:
            self.num_pdb, self.valid_pdbs, self.no_ligand_pdbs, self.not_valid_pdbs, self.pdb_mono_chains = self.loadObj(fp)
        self.prnFilteredbySeq()
    
    
    def getMHCIPDBS(self):
        """ return existing MHCI pdbids. if not existed, return empty set """
        fp = self.getMHCIPDBFpBin()
        if not self.isNew(fp):
            mhcI_pdbs = self.loadObj(fp)  
        else:
            mhcI_pdbs = {}
        return mhcI_pdbs, fp  
    
    
    def getCheckPDBIDs(self):
        fn = "checked_pdbids.bin"
        fp = self.joinPath(self.mhcIdb_pdb_path, fn)
        if not self.isNew(fp): 
            checked_pdbids, mhcI_pdbids  = self.loadObj(fp)
            if isinstance(checked_pdbids, list):
                checked_pdbids = set(checked_pdbids)
                mhcI_pdbids = set(mhcI_pdbids)
        else:
            checked_pdbids = set([])
            mhcI_pdbids = set([])
        return checked_pdbids, mhcI_pdbids, fp
    
        
    def filteredBy3D(self):
        mhcI_pdbs, fp_mhcI_pdbs  = self.getMHCIPDBS()
        checked_pdbids, mhcI_pdbids, fn = self.getCheckPDBIDs()
        pdbids = self.valid_pdbs + self.no_ligand_pdbs
        num_new_pdb = 0 
        cnt = 0 
        for pdbid in pdbids: 
            if pdbid in checked_pdbids or pdbid in mhcI_pdbs:
                 continue 
            pdbid = pdbid.upper()
            fp = self.pdbid2Fp(pdbid)
            t_obj = MHCIPDB(fn_pdb=fp, loadPDB=1)
            is_mhcI, lns = t_obj.start(mono_chain_ids=self.pdb_mono_chains[pdbid]) #print(fp) 
            if is_mhcI:
                mhcI_pdbids.add(pdbid)
                num_new_pdb += 1
                mhcI_pdbs[pdbid] = t_obj.getMHCIPDBInf()
                self.add2log(lns)
            else:
                ln = "# {} not a MHCI PDB with peptide binding groove".format(fp)
                self.add2log(ln, prn=1)
            cnt += 1
            checked_pdbids.add(pdbid)
            if cnt == 50:
                self.dumpObj((checked_pdbids,mhcI_pdbids), fn)
                cnt = 0 
        if cnt > 0:
            self.dumpObj((checked_pdbids, mhcI_pdbids), fn)
            #break
        ln = "# num_mhcI_pdb:{}".format(len(mhcI_pdbids))
        self.add2log(ln, prn=1)
        if num_new_pdb > 0:
            self.dumpObj(mhcI_pdbs, fp_mhcI_pdbs, vb=1)
            self.add2log(ln)
            self.saveLog("checked_pdbs_{}.log".format(num_new_pdb))
    
    
    def testPDB(self, pdbid, vb=0):
        pdbid = pdbid.upper()
        fp = self.pdbid2Fp(pdbid)
        t_obj = MHCIPDB(fn_pdb=fp, loadPDB=1)
        is_mhcI, lns = t_obj.start(mono_chain_ids=self.pdb_mono_chains[pdbid], vb=vb)          
        self.prnLines(lns)
        print("# isMHCI: {}".format(is_mhcI))
        
        
    def statMHCIPDB(self):
        fn = "checked_pdbids.bin"
        checked_pdbids, mhcI_pdbids  = self.loadObj(fn)
        num_mhcI_pdb = len(mhcI_pdbids)
        print("# num_mhcI_pdb: {}".format(num_mhcI_pdb))
        for pdbid in mhcI_pdbids:
            self.testPDB(pdbid, vb=1)
            break 
        
            
    def start(self):
        #self.UPDATE = True
        self.loadFilteredPDBbySeq() 
        #self.filteredBy3D() 
        #self.statMHCIPDB()
        #self.testPDB("2c7u", vb=1)
        #self.testPDB("4QRQ", vb=1)
        #self.testPDB("1M05", vb=1)

        
if  __name__ == "__main__":
    import sys, argparse
    parser = argparse.ArgumentParser(prog=__file__, description=__doc__)
    parser.add_argument("-v","--version", action='version', version='%(prog)s 1.0')
   
    parser.add_argument("-d", "--data_dir",dest="data_dir", default="example_data",
                        help="input existed  data directory")         
    
    parser.add_argument("-u","--update", dest="UPDATE",  action="store_true",  default=False,
                        help="Update all existing binary files")       
      
    args =parser.parse_args()
    r_obj = PDB2MHCI(args)
    r_obj.start()        
