from .libio import MIO
class MHCIDBData(MIO):
    """ existing data path and file names 
    """
    def setMHCIDBPath(self, mhcidb_dirname="MHCIDB"):
        """ set MHCIDB working directory """
        cwd = self.getCWD() 
        pre_dir, after = cwd.split(mhcidb_dirname)
        self.mhcIdb_path = self.joinPath(pre_dir,mhcidb_dirname )
        print(("# MHCIDB workding path: {}".format(self.mhcIdb_path)))
        
        self.mhcIdb_existing_data_path = self.joinPath(self.mhcIdb_path, "existing_data")
        self.mhcIdb_hla_path = self.joinPath(self.mhcIdb_existing_data_path, "hla" )
        self.mhcIdb_pdb_path = self.joinPath(self.mhcIdb_existing_data_path, "pdb" )
        self.mhcIdb_ba_path = self.joinPath(self.mhcIdb_existing_data_path, "ba" ) 
        self.mhcIdb_pdb3d_path = self.joinPath(self.mhcIdb_pdb_path, "raw_pdbs")
        
    def get_hla_aligned_seq_fp(self):
        """ return the path of the file contains the alinged protein seqeuences of 
        HLA gene A, B and C
        """
        fn_aln = "ClassI_prot.txt"
        return self.joinPath(self.mhcIdb_hla_path, fn_aln)    
 
    def getAnchorMajorSeqFp(self):
        fn = "HLA_amseq.bin"
        return self.joinPath(self.mhcIdb_hla_path, fn)
        
    
    def getHLAAlleleGrpFp(self):
        fn = "hla_allele_grps.bin"
        return self.joinPath(self.mhcIdb_hla_path, fn)
    
    
    def getDLFastaFp(self, fn="fasta.txt"):
        return self.joinPath(self.mhcIdb_pdb_path, fn)
    
    
    def pdbid2Fp(self, pdbid):
        return self.joinPath(self.mhcIdb_pdb3d_path, "{}.pdb".format(pdbid.lower()))
    
    
    def getSeqFilteredFpBin(self):
        fn_out_pre = "mhcI_filter_by_seq" 
        fn_out_pre = self.joinPath(self.mhcIdb_pdb_path,fn_out_pre)
        fp_out_bin = "%s.bin" % fn_out_pre
        return fp_out_bin
    
    def getMHCIPDBFpBin(self):
        fn = "mhcI_pdbs.bin"
        return self.joinPath(self.mhcIdb_pdb_path, fn)
    
    def getPDB2ALLeleFpBin(self):
        fn = "mhcI_pdb_to_allele.bin"
        return self.joinPath(self.mhcIdb_path, fn)
        
    def getIEDBFp(self):
        fn_bind_data = 'bdata.20130222.mhci.txt'
        return self.joinPath(self.mhcIdb_ba_path, fn_bind_data)
    
    def getBindDataFp(self):
        fn_bind_bin = "mhcI_bdata.bin"
        return self.joinPath(self.mhcIdb_ba_path, fn_bind_bin)
    
    def loadMHCIBindData(self):
        fp_bind_data = self.getBindDataFp()
        if not self.isNew(fp_bind_data):
            self.mhcI_bind_data =  self.loadObj(fp_bind_data)
        else:
            self.readBindData(self.getIEDBFp())
            self.dumpObj(self.mhcI_bind_data, fp_bind_data)
    
    
    def readBindData(self, fp):
        txt = self.readTxtFile(fp)
        self.mhcI_bind_data = {}
        cnt = 0
        for ln in txt[1:]:
            ln = ln.strip()
            elems = ln.split()
            if not ln: continue
            if len(elems) != 6:
                print("# Error in parsing: %s" % ln) 
            else:
                allele_name = elems[1]
                if allele_name.startswith("HLA-"):
                    cnt += 1
                    lig_len  = elems[2]
                    lig = elems[3]
                    inq = elems[4]
                    bd = elems[5]
                    if allele_name in self.mhcI_bind_data:
                        self.mhcI_bind_data[allele_name].append((lig_len, lig, inq, bd))
                    else:
                        self.mhcI_bind_data[allele_name]= [(lig_len, lig, inq, bd)]
        ln = "# num HLA binding data: %s" % cnt
        self.add2log(ln, vb=1) 
    
    
    def loadMHCIPDB2AlleleData(self):
        fp_bin  = self.getPDB2ALLeleFpBin()
        self.matched_pdbs, self.non_matched_pdbs = self.loadObj(fp_bin)
                    
     

    def loadMHCIPDBInf(self):
        """ mhci_pdbs[pdbid] = [(chain_ids, chain_seq(xxx,...), [ligand chains]), ...] """
        fp = self.getMHCIPDBFpBin()
        self.mhcI_pdbs = self.loadObj(fp)
        self.mhcI_pdbids = self.mhcI_pdbs.keys() 
        
    def getMHCILigandFpBin(self):
        fn = "mhcI_pdb_ligand.bin"
        fp = self.joinPath(self.mhcIdb_pdb_path, fn)
        return fp 
    