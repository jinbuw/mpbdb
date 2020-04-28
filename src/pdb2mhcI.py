# find mhc II in pdb file by the geometry of MHC class II (two helices and 8 beta sheets)
""" select MHCI protein 3D structures with the peptide binding groove """

from lib.libpdb import MHCIPDB
from lib.libmhcidb import MHCIDBData
from lib.libfasta import FilterMHCbySeq

class PDB2MHCI(FilterMHCbySeq, MHCIDBData):
    """ find MHC  MHC class II  in pdb file """
    def __init__(self, args):
        """Stores the parameters for an alignment scoring function"""
        self.args = args
        self.fn_inp = args.fn_inp
        self.setMHCIDBPath()

    def setPDBID(self):
        fn_pdb = self.getBaseName(self.fn_inp)
        fn_pre = self.getFnPre(fn_pdb)
        pdbid = self.args.pdbid     
    
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
        fp = self.getMHCIPDBFpBin()
        if self.isExistFile(fp):
            if hasattr(self, "UPDATE") and self.UPDATE:
                pass 
            else:
                return False, None 
            mhcI_pdbs = self.loadObj(fp)  
        else:
            mhcI_pdbs = set({})
        return True, mhcI_pdbs, fp  
    
    
    def getCheckPDBIDs(self):
        fn = "checked_pdbids.bin"
        if not self.isNew(fn): 
            checked_pdbids, mhcI_pdbids  = self.loadObj(fn)
            if isinstance(checked_pdbids, list):
                checked_pdbids = set(checked_pdbids)
                mhcI_pdbids = set(mhcI_pdbids)
        else:
            checked_pdbids = set([])
            mhcI_pdbids = set([])
        return checked_pdbids, mhcI_pdbids, fn
        
    def filteredBy3D(self):
        NEW_MHCI_PDBS, mhcI_pdbs, fp_mhcI_pdbs  = self.getMHCIPDBS()
        if not NEW_MHCI_PDBS:
            return 
        checked_pdbids, mhcI_pdbids, fn = self.getCheckPDBIDs()
        pdbids = self.valid_pdbs + self.no_ligand_pdbs
        #print(mhcI_pdbids)
        num_new_pdb = 0 
        cnt = 0 
        for pdbid in pdbids: 
            #if pdbid in checked_pdbids:
            #     continue 
            if pdbid not in mhcI_pdbids:
                continue
            pdbid = pdbid.upper()
            #print("# checked pdb: {}".format(pdbid))
            fp = self.pdbid2Fp(pdbid)
            t_obj = MHCIPDB(fn_pdb=fp, loadPDB=1)
            is_mhcI, lns = t_obj.start(mono_chain_ids=self.pdb_mono_chains[pdbid]) #print(fp) 
            if is_mhcI:
                #mhcI_pdbids.append(pdbid)
                mhcI_pdbids.add(pdbid)
                num_new_pdb += 1
                mhcI_pdbs[pdbid] = t_obj.getMHCIPDBInf()
                #print(mhcI_pdbs)
                self.add2log(lns)
                #break 
            cnt += 1
            checked_pdbids.add(pdbid)
            #checked_pdbids.append(pdbid)
            if cnt == 50:
                self.dumpObj((checked_pdbids,mhcI_pdbids), fn)
                cnt = 0 
        if cnt > 0:
            self.dumpObj((checked_pdbids, mhcI_pdbids), fn)
            #break
        ln = "# num_mhcI_pdb:{}".format(len(mhcI_pdbids))
        print(ln)
        if num_new_pdb > 0:
            self.dumpObj(mhcI_pdbs, fp_mhcI_pdbs, vb=1)
            self.add2log(ln)
            self.saveLog("checed_pdbs_{}.log".format(num_new_pdb))
    
    
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
        self.testPDB("4QRQ", vb=1)
        self.testPDB("1M05", vb=1)

        
if  __name__ == "__main__":
    import sys, argparse
    parser = argparse.ArgumentParser(prog=__file__, description=__doc__)
    parser.add_argument("-v","--version", action='version', version='%(prog)s 1.0')
    
     
    parser.add_argument("-f", "--filename", dest="fn_inp", default=None, 
                                help="input fasta filename")     
    
    parser.add_argument("-P", "--pdbid", dest="pdbid", default=None, help="PDB ID")
    
    parser.add_argument("-d", "--dir", dest="pdb_dir", default=None, help="PDB directory")

    args =parser.parse_args()
    r_obj = PDB2MHCI(args)
    r_obj.start()        
