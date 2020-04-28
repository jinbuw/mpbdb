# map MHCI pdb with allele name to binding affinity 
from lib.libmhcidb import MHCIDBData 
from lib.libfasta import FASTA

class MHCIPDB2BA(MHCIDBData, FASTA):
    def statMHCIBindAllele(self):
        alleles = list(self.mhcI_bind_data.keys())
        allele_genes = {}
        for allele_name in alleles:
            if "-" not in allele_name:
                print("# wrong allele name format: %s " %  allele_name)
            else:
                _allele_gene = allele_name[:5]
                allele_genes.setdefault(_allele_gene, []).append(allele_name)
        num_allele_proteins = len(alleles)
        lns = self.getCommandLines()
        lns.append( "#num_allele_protein_in_bdata: %4d" % num_allele_proteins)
        genes = list(allele_genes.keys())
        genes.sort()
        tot_bd = 0 
        for _gene in genes:
            _alleles = allele_genes[_gene]
            num_bd_in_gene = 0
            for _allele in _alleles:
                num_bd_in_gene += len(self.mhcI_bind_data[_allele])
            lns.append("gene: %s  num_allele_protein: %4d  num_bd: %4d "  % (_gene, len(_alleles), num_bd_in_gene))
            tot_bd += num_bd_in_gene
        lns.append( "#total_number_bd: %8d" %  tot_bd )
        self.prnLines(lns)
        self.saveLines(lns,"mhcI_bind_gene.txt", vb=1)
        
    
    def statMHCIBind(self, lig_len='9'):
        """ statistically analyze MHCI binding data """
        #print self.mbl.keys()
        lig_bind = self.mbl[lig_len]
        alleles= list(lig_bind.keys())
        alleles.sort()
        print("# ligand lenght=%s" % lig_len) 
        tot = 0 
        for allele in alleles:
            num_bd = len(lig_bind[allele])
            tot += num_bd            
            print("# %s  : %s"  % (allele, num_bd))
        print("#tot= %6d" % tot) 
    
    
    def statMHCBindLiglen(self):
        if not hasattr(self, "mbl"):
            self.setLiglenBD() 
            
        lig_lens = list(self.mbl.keys())
        lig_lens = self.intli(lig_lens)
        lig_lens.sort()
        lns = ["# MHCI binding ligand length vs number of ligands"]
        tot_num_lig = 0
        for lig_len in lig_lens:
            lig_bind = self.mbl[str(lig_len)]
            alleles = list(lig_bind.keys())
            tot = 0 
            for allele in alleles:
                num_bd = len(lig_bind[allele])
                tot += num_bd
            tot_num_lig += tot 
            lns.append("lig_len= %4d  num_bind= %4d" % (lig_len, tot))
        lns.append("# tot_num_lig(human)= %8d" %  tot_num_lig )
        out_lns = self.getCommandLines() + lns 
        self.prnLines(out_lns)
        self.saveLines(out_lns, "mhcI_binding_ligand_lens.txt", vb=1)
   

    def setLiglenBD(self):
        """ set ligand length's binding data """
        self.mbl = {}
        for allele_name in self.mhcI_bind_data:
            #print( self.mhcI_bind_data[allele])
            for lig_len, lig, inq, bd in self.mhcI_bind_data[allele_name]:
                #lig_len, lig, inq, bd = self.mhcI_bind_data[allele]
                if lig_len in self.mbl:
                    if allele_name in self.mbl[lig_len]:
                        self.mbl[lig_len][allele_name].append((lig, inq, bd))
                    else:
                        self.mbl[lig_len][allele_name] = [(lig, inq, bd)]
                else:
                    self.mbl[lig_len] = {}
                    self.mbl[lig_len][allele_name] = [(lig, inq, bd)]
        
        
    def transBD2Sub(self):
        self.mbl_sub = {}
        for lig_len in list(self.mbl.keys()):
            self.mbl_sub[lig_len] = {}
            for allele_name in list(self.mbl[lig_len].keys()):
                asub = self.alleleName2Sub(allele_name)
                for lig, inq, bd in self.mbl[lig_len][allele_name]:
                    if asub in self.mbl_sub[lig_len]:
                        self.mbl_sub[lig_len][asub].append((allele_name, lig, inq, bd))
                    else:
                        self.mbl_sub[lig_len][asub] = [(allele_name, lig, inq, bd)]
                
    def loadMHCIPDBLigandLength(self):
        fn = "mhcI_pdb_ligand_length.bin"
        self.mhcI_pdb_ligand = self.loadObj(fn)
        #print self.mhcI_pdb_ligand.keys()
      
 
    def alleleProt2BindAllele(self, name):
        return "HLA-{}".format(name)
    
    
    def MapPDB2Bind(self):
        """ map MHCI PDB to binding affinity """
        self.mhcI_pdb_2_bd = {}
        num_w_lig = 0 
        self.mhcI_pdb_bd_liglen  = {}
        for pdbid in self.matched_pdbs:
            if len(self.matched_pdbs[pdbid]) == 2:
                 continue 
            ma_names, mns, mdifs = self.matched_pdbs[pdbid]
            #print("#pdbid: {} ma_names: {}".format(pdbid, ma_names))
            try:
                lig_seq, lig_chain_id, lig_len = self.getLigandSeq(pdbid)
            except TypeError:
                lig_seq = None
            #if lig_seq:
            #    print("# lig_seq: {} lig_len: {}".format(lig_seq, lig_len))    
            pdb_bds = []
            cnt = 0 
            for name in ma_names:
                ba_name = self.alleleProt2BindAllele(name)
                if ba_name not in self.mhcI_bind_data:
                    continue 
                if lig_seq:
                    #print("# ba_name: {}".format(ba_name))
                    mbds = self.mapBindLigand(ba_name, lig_seq, lig_len)
                    num_mbds = len(mbds)
                    cnt += num_mbds
                    if num_mbds > 0:
                        for bd in mbds:
                            lig_len = int(bd[0])
                            if lig_len in self.mhcI_pdb_bd_liglen:
                                if pdbid in self.mhcI_pdb_bd_liglen[lig_len]:
                                    self.mhcI_pdb_bd_liglen[lig_len][pdbid].append(bd)
                                else:
                                    self.mhcI_pdb_bd_liglen[lig_len][pdbid] = [bd]
                            else:
                                self.mhcI_pdb_bd_liglen[lig_len] = {pdbid:[bd]}
                    pdb_bds.append((ba_name, mbds))
                else:
                    pdb_bds.append((ba_name, []))
            
            if pdb_bds:
                self.mhcI_pdb_2_bd[pdbid] = pdb_bds, cnt 
                #print(pdbid, pdb_bds)
                if cnt > 0:
                    print("#*** matched mhcI-pep: {} lig_len:{}  lig_seq:{}".format(pdbid, lig_len, lig_seq))
                    num_w_lig += 1
                    #print(pdb_bds)
                    #break 
                    if cnt > 1:
                        print("#multi_bind: {} -> {}".format(pdbid,pdb_bds))
        num_m_pdb = len(self.mhcI_pdb_2_bd)
        ln = "# num_matched_pdb: {}  num_match_mhcI-pep: {}".format(num_m_pdb, num_w_lig)
        #print(ln)
        self.add2log(ln, prn=1)
        fn = "mhcI_pdb_2_bind.bin"
        self.dumpObj(self.mhcI_pdb_2_bd, fn, vb=1)
        lig_lens = list(self.mhcI_pdb_bd_liglen.keys())
        lig_lens.sort()
        for lig_len in lig_lens:
            ln = "# pdb2bd: lig_len={}  num_pdb: {}".format(lig_len, len(self.mhcI_pdb_bd_liglen[lig_len]))
            self.add2log(ln, prn=1)
        fn = "mhcI_pdb_bd_liglen.bin"
        self.dumpObj(self.mhcI_pdb_bd_liglen, fn, vb=1)
      
    def testMap(self):
        """ test map MHCI PDB to binding affinity """
        self.mhcI_pdb_2_bd = {}
        num_w_lig = 0 
        self.mhcI_pdb_bd_liglen  = {}
        #pdbids = ["4QRQ", "1M05", "4U1J", "4U1M", "3FT2"]
        pdbids = ["3FT2"]
        for pdbid in pdbids:
            if len(self.matched_pdbs[pdbid]) == 2:
                print("# unmatched pdb: {}".format(pdbid))
                continue 
            ma_names, mns, mdifs = self.matched_pdbs[pdbid]
            #print("#pdbid: {} ma_names: {}".format(pdbid, ma_names))
            try:
                lig_seq, lig_chain_id, lig_len = self.getLigandSeq(pdbid)
            except TypeError:
                lig_seq = None
            if lig_seq:
                lig_len_in_bd_txt = int(lig_len)
                lig_len = len(lig_seq)
                if lig_len_in_bd_txt != lig_len:
                    print("#*** Error in bd data:  lig_seq: {} lig_len_in_bdtxt: {}  lig_len(seq): {}".format(lig_seq, lig_len, len(lig_seq)))    
            pdb_bds = []
            cnt = 0 
            for name in ma_names:
                ba_name = self.alleleProt2BindAllele(name)
                if ba_name not in self.mhcI_bind_data:
                    continue 
                if lig_seq:
                    #print("# ba_name: {}".format(ba_name))
                    mbds = self.mapBindLigand(ba_name, lig_seq, lig_len)
                    num_mbds = len(mbds)
                    cnt += num_mbds
                    if num_mbds > 0:
                        for bd in mbds:
                            lig_len = int(bd[0])
                            if lig_len in self.mhcI_pdb_bd_liglen:
                                if pdbid in self.mhcI_pdb_bd_liglen[lig_len]:
                                    self.mhcI_pdb_bd_liglen[lig_len][pdbid].append(bd)
                                else:
                                    self.mhcI_pdb_bd_liglen[lig_len][pdbid] = [bd]
                            else:
                                self.mhcI_pdb_bd_liglen[lig_len] = {pdbid:[bd]}
                    pdb_bds.append((ba_name, mbds))
                else:
                    pdb_bds.append((ba_name, []))
            
            if pdb_bds:
                self.mhcI_pdb_2_bd[pdbid] = pdb_bds, cnt 
                #print(pdbid, pdb_bds)
                if cnt > 0:
                    print("#*** matched mhcI-pep: {} lig_len:{}  lig_seq:{}".format(pdbid, lig_len, lig_seq))
                    num_w_lig += 1
                    #print(pdb_bds)
                    #break 
                    if cnt > 1:
                        print("#multi_bind: {} -> {}".format(pdbid,pdb_bds))
     
        
    def mapBindLigand(self, allele_bd, lig_seq, lig_len):
        #print("# pdb: {}, lig_seq:{}".format(pdbid, lig_seq))
        bds = self.mhcI_bind_data[allele_bd]
        mbds = []
        for bd in self.mhcI_bind_data[allele_bd]:
            #print("# bd: {}".format(bd))
            if int(bd[0]) == lig_len and bd[1] ==lig_seq:
                    mbds.append(bd)
        #print(lig_seq, mbds)
        return mbds
    
    
    def getLigandSeq(self, pdbid):
        if pdbid in self.mhcI_pdb_ligands_seqs:
            return self.mhcI_pdb_ligands_seqs[pdbid]
        
        
    def loadMHCIPDBLigandSeq(self):
        fp_bin = self.getMHCILigandFpBin()
        if self.isNew(fp_bin):
            fp_fasta = self.getDLFastaFp()
            self.readPdbSeqs(fp_fasta)
            #print( self.pdb_seqs)
            self.mhcI_pdb_ligands_seqs = {}
            for pdbid in self.mhcI_pdbids:
                ligand_inf = self.mhcI_pdbs[pdbid][-1]
                if ligand_inf:
                    ligand_chain_id = ligand_inf[0][0]
                    ligand_length = ligand_inf[0][1]
                    self.mhcI_pdb_ligands_seqs[pdbid] = self.pdb_seqs[pdbid][ligand_chain_id], ligand_chain_id, ligand_length
            self.dumpObj(self.mhcI_pdb_ligands_seqs, fp_bin, vb=1)
        else:
            self.mhcI_pdb_ligands_seqs = self.loadObj(fp_bin)  
            
    
    
    def listPDBIDS(self):
        fn = "mhcI_pdb_2_bind.bin"
        self.mhcI_pdb_2_bd = self.loadObj(fn)
        pdbids = list(self.mhcI_pdb_2_bd.keys())
        pdbids.sort()
        cnt = 0
        num_pdb = len(pdbids)
        lns = ["# num_pdb= {}".format(num_pdb)]
        
        num_col = 12
        for i in range(0, num_pdb, num_col):
            end_sid = i + num_col 
            ln = " ".join(pdbids[i:end_sid])
            lns.append(ln)
        lns.append(" ".join(pdbids[end_sid:]))
        self.prnLines(lns)
    
    def prnBlockList(self, pdbids,  num_col=15):
        lns = []
        #num_col = 12
        num_pdb =len(pdbids)
        for i in range(0, num_pdb, num_col):
            end_sid = i + num_col 
            ln = " ".join(pdbids[i:end_sid])
            lns.append(ln)
        lns.append(" ".join(pdbids[end_sid:]))
        #self.prnLines(lns)
        return lns 
        
    def listPdbPep(self):
        """ list pdbids with binding peptides """
        fn = "mhcI_pdb_bd_liglen.bin"
        self.mhcI_pdb_bd_liglen = self.loadObj(fn, vb=1)
        pdbids = []
        for liglen in self.mhcI_pdb_bd_liglen:
            pdbids += list(self.mhcI_pdb_bd_liglen[liglen].keys())
        pdbids.sort()
        lns = ["# pdbids with peptides: num_pdb={}".format(len(pdbids))]
        lns += self.prnBlockList(pdbids)
        self.prnLines(lns)
        
    def start(self):
        self.setMHCIDBPath()
        self.loadMHCIBindData()
        self.statMHCIBindAllele()
        self.statMHCBindLiglen()
        self.loadMHCIPDB2AlleleData()
        self.loadMHCIPDBInf()
        self.loadMHCIPDBLigandSeq()
        #self.MapPDB2Bind()
        #self.testMap()
        #self.listPDBIDS() 
        self.listPdbPep()
        
if  __name__ == "__main__":
    import argparse, sys 
    cmd = " ".join(sys.argv)
    cmd_lns = [ "\n---------------------------"]
    cmd_lns += [ "Runing: python %s" % cmd, "#    aligned two sequences with affine gap penalty" ]
    cmd_lns += ["---------------------------\n"]
        
    parser = argparse.ArgumentParser(prog="%s" % __file__, description="cluster pdbs by their sequences")
    parser.add_argument("-v","--version", action='version', version='%(prog)s 1.0')
    
    parser.add_argument("-f", "--filename",dest="fn_fasta", default=None,
                        help="input fasta filename")             

    parser.add_argument("-O","--printout", dest="prnout",  action="store_true",  default=False,
                        help="print aligned sequence")    
    args =parser.parse_args()    
    
    tobj = MHCIPDB2BA(args)
    tobj.start()
          