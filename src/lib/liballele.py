# class and functions for loading allele related parameters 

#print "# num_pdb2allele=%4d  num_allele= %s" % (num_pdb, num_allele)
    #self.dumpObj(self.mhcI_allele2pdb, "mhcI_allele2pdb.bin")
from .libio import MIO
class Allele(MIO):
    bin_path = r'/home/jiw314/rcse/mhcI/mhcI_pdb_by_seq'
    
    def getBinFp(self, fn):
        return self.joinPath(self.bin_path, fn)
    
    def allele2pdb(self, allele):
        if not hasattr(self, "mhcI_allele2pdb"):
            self.loadAllele2PDB()
        pdbids =  self.mhcI_allele2pdb[allele]
        print("# allele: %s num_pdb: %4d" % (allele, len(pdbids)))
        return pdbids 
            
    def loadAllele2PDB(self):
        self.mhcI_allele2pdb = self.loadObj(self.getBinFp("mhcI_allele2pdb.bin"))

  
    def filterPdb(self, fns_pdb, pdbids):
        low_pdbids = []
        for pdbid in pdbids:
            low_pdbids.append(pdbid.lower()) 
        tgt_fns_pdb = []
        for fn_pdb in fns_pdb:
            for pdbid in low_pdbids:
                #print pdbid, fn_pdb 
                if pdbid in fn_pdb:
                    tgt_fns_pdb.append(fn_pdb)
        return tgt_fns_pdb    