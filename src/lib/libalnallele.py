# Building HLA allele tree from the aligned MHCI sequence 
#  based on CLASSI_prot.txt

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62
gap_open = -10
gap_extend = -0.5

from string import digits 
from .libmhcidb import MHCIDBData

class HLAName:
    def alleleName2Protein(self, allele_name):
        """ return allele_protein and allele_gene_group """
        elems = allele_name.split(":")
        allele_grp = elems[0]
        allele_protein = ":".join(elems[:2])
        return allele_protein, allele_grp 
    
    
    def alleleGrp2Gene(self, allele_grp):
        if "-" in allele_grp:
            tmp, _inf = allele_grp.split("-")
        elif "_" in allele_grp:
            tmp, _inf = allele_grp.split("_")
        else:
            _inf = allele_grp
        if "*" in _inf:
            gene, grp = _inf.split("*")
        else:
            gene = _inf 
        return gene 
    
    def getAlleleProtVal(self, allele_name):
        """ B*08:103 """
        elems = allele_name.split(":")
        return int(elems[1])
    
    def isChangedExpression(self, allele_name):
        """ if the last charater of name is a letter, it is a protein sequence with changes in expression """
        if allele_name[-1] in digits:
            return False
        else:
            return True
    

class HLAAllele(MHCIDBData, HLAName):
    MAXLEN = 369  # some of alignments in B only with 369, missing last 3 
    hla_genes = ["A","B","C"]
    def __setAlnSeq(self, fp):
        lns = self.readTxtFile(fp)
        self.aln_seqs = {}
        for ln in lns:
            ln  = ln.strip()
            if not ln: continue 
            if len(ln) < 2:
                 continue
            if ln[1] == "*":
                elems = ln.split()
                allele_name = elems[0]
                seq = "".join(elems[1:])
                if allele_name in self.aln_seqs:
                    self.aln_seqs[allele_name] += seq 
                else:
                    self.aln_seqs[allele_name] = seq
                    
    def loadAlignedSeq(self):
        self.fn_aln = self.get_hla_aligned_seq_fp()
        fp_bin = self.fp2FpBin(self.fn_aln)
        self.add2log("# aligned_seq_bin: {}".format(fp_bin))
        if self.isNew(fp_bin):
            self.__setAlnSeq(self.fn_aln)
            self.dumpObj(self.aln_seqs, fp_bin, vb=1)
        else:
            self.aln_seqs = self.loadObj(fp_bin)
        self.num_allele = len(list(self.aln_seqs.keys()))
        self.add2log("# num_allele in aligned file: %s (%s) " % (self.num_allele, self.fn_aln))
    
    
    def __setAlleleGrps(self):
        """ set allele gene groups from allele name
            duplicates in protein names 
        """
        allele_names = list(self.aln_seqs.keys())
        
        self.allele_grp_uniq_proteins = {}
        self.allele_prot_names = {} # allele gene group to allele names  
        self.gene_allele_grp = {} 
        for allele_name in allele_names:
            if self.isChangedExpression(allele_name):
                continue 
            allele_protein, allele_grp = self.alleleName2Protein(allele_name)
            allele_gene = self.alleleGrp2Gene(allele_grp)
            self.gene_allele_grp.setdefault(allele_gene, set([])).add(allele_grp)
            self.allele_grp_uniq_proteins.setdefault(allele_grp, set([])).add(allele_protein)
            self.allele_prot_names.setdefault(allele_protein, []).append(allele_name) 
        for gene in self.hla_genes:
            self.gene_allele_grp[gene] = self.sortNamebyLastSetDigit(list(self.gene_allele_grp[gene])) 
        for grp in self.allele_grp_uniq_proteins:
            self.allele_grp_uniq_proteins[grp] = self.sortNamebyLastSetDigit(list(self.allele_grp_uniq_proteins[grp]), sep=":")
            
    
    def loadAlleleGrps(self):
        fp_allele_grp = self.getHLAAlleleGrpFp()
        self.add2log(["# allele_grp: {}".format(fp_allele_grp)])
        if self.isNew(fp_allele_grp):
            self.__setAlleleGrps()
            self.dumpObj((self.gene_allele_grp, self.allele_grp_uniq_proteins, self.allele_prot_names), fp_allele_grp, vb=1)
        else:
            (self.gene_allele_grp, self.allele_grp_uniq_proteins, self.allele_prot_names) = self.loadObj(fp_allele_grp, vb=1)
            
        
    def loadHLAAllele(self, comm_reg=(25,205)):
        """  common region in allele protein sequences """
        self.loadAlignedSeq()
        self.setHLAAlignedSeqDt()
        self.__setAlleleGrps() 
        self.setRefSeq()
        self.setCommonSeqReg(comm_reg)
        self.loadAlleleGrps()
        self.loadHLAAMSeq()
        
        
    def setHLAAlignedSeqDt(self):
        self.aligned_alleles_seqdt = {}
        
                
    def getProtAlleles(self, prot):
        """ specific protein has multiple alleles in hla allels """
        return self.allele_prot_names[prot]
    
    def getHLAProts(self):
        return list(self.allele_prot_names.keys()) 
    
    
        
    def getGeneGrps(self,gene):
        return self.gene_allele_grp[gene]
    
    def getGrpProts(self, grp):
        """ return group's uniq proteins """
        return self.allele_grp_uniq_proteins[grp]
        
 
    def getProtAllele(self, prot):
        return self.allele_prot_names[prot][0]
    
    
    def getProtCommSeq(self, prot):
        try:
            return self.hla_prot_seq[prot]
        except KeyError:
            allele = self.getProtAllele(prot)
            prot_seq = self.getAlleleCommRealSeqDt(allele)
            self.hla_prot_seq[prot] = prot_seq
            return prot_seq
    
        
    def getGrpSeq(self, grp):
        prot = self.allele_grp_uniq_proteins[grp][0]
        return self.getProtSeq(prot)
    
    def getGrpSeqProt(self, grp):
        prot = self.allele_grp_uniq_proteins[grp][0]
        return prot, self.getProtSeq(prot)
    
    def getAlleleCommRealSeq(self, allele):
        new_seq = ""
        for i in self.getCommSeqRegion():
            new_seq += self.getRealPosSeq(allele, i)
        return new_seq
    
    def getAlleleCommRealSeqDt(self, allele):
        """ return sequence in common region with dictionary format """
        new_seq = {}
        for i in self.getCommSeqRegion():
            new_seq[i] = self.getRealPosSeq(allele, i)
        return new_seq
    
    
    def getAlleleAlignedSeq(self, allele_name):
        return self.aln_seqs[allele_name]


    def setRefSeq(self, ref_allele="A*01:01:01:01"):
        if not hasattr(self, "aln_seqs") or self.aln_seqs is None:
            self.loadAlignedSeq()
        self.ref_aseq = self.aln_seqs[ref_allele]
        
    def getRefSeq(self):
        return self.ref_aseq 
    
    
    def getRealPosSeq(self, allele_name, posid):
        pos_seq = self.aln_seqs[allele_name][posid]
        if pos_seq in [ "-", "*"]:
            return self.ref_aseq[posid]
        elif pos_seq in  ["."]: 
            return " "
        else:
            return pos_seq
    
    def __getAlignedAlleleSeqDt(self, allele):
        """ set allele's seqeunce dt with posid as key """
        seq_dt = {}
        aln_seq = self.aln_seqs[allele]
        for i in range(self.MAXLEN):
            aln_pos_seq = aln_seq[i]
            nt = self.__getAlignedPosSeq(aln_pos_seq, i)
            if nt:
                seq_dt[i] = nt 
        return seq_dt
       
       
    def __getAlignedPosSeq(self, aln_pos_seq, posid):
        """ return pos's sequence in the aligned seqeunce """
        if aln_pos_seq in [ "-"]: 
            return self.ref_aseq[posid]
        elif aln_pos_seq in  [".", "*"]: 
            return 
        else:
            return aln_pos_seq
        
        
    def getAlignedAlleleSeqDt(self, allele):
        try:
            return self.aligned_alleles_seqdt[allele]
        except KeyError:
            self.aligned_alleles_seqdt[allele] = self.__getAlignedAlleleSeqDt(allele)
            return self.aligned_alleles_seqdt[allele]
        
    def statProtAlleleLength(self, prot):
        """ stat prot's allele protein sequence lenghts in hla allele database """
        alleles = self.getProtAlleles(prot)
        print(("# prot: {} -> alleles:{}".format(prot, alleles)))
        for allele in alleles:
            allele_seq_len = len(self.hla_aligend_seq_dt(allele))
            print(("# allele: {} len: {}".format(allele, allele_seq_len)))
    
    
    def __getDescProtAlleles(self, prot, MAXNUM=5):
        """ return alleles with the same specific protein name, descending in the length of the allele protein sequence"""
        alleles = self.getProtAlleles(prot)
        allele_seqs = []
        for allele in alleles:
            allele_seqs.append((len(self.getAlignedAlleleSeqDt(allele)),allele))
        allele_seqs.sort(reverse=True)
        num_allele = len(alleles)
        min_num = min(num_allele, MAXNUM)
        desc_alleles = []
        for seq_len, allele in allele_seqs[:min_num]:
            desc_alleles.append(allele) 
        return desc_alleles     

        
    def statHLAUniqProtSeq(self, prot):
        self.__setProtSeqStat()
        desc_prot_alleles  = self.__getDescProtAlleles(prot)
        for allele in desc_prot_alleles:  
            self.__addSeq2SeqStat(self.getAlignedAlleleSeqDt(allele))
        #print "# seq_stat: ", len(self.seq_stat)
        aseq,  mseq, am_stat  =  self.getCommSeq(COMMREG=False)
        #print("# prot: {}  num_anchor: {}".format(prot, len(aseq)))
        return aseq 
    
    
    def aseq2Str(self, aseq):
        sids = list(aseq.keys())
        sids.sort()
        seq_str = ""
        for sid in sids:
            seq_str += aseq[sid]
        return seq_str
    
    
    def __addSeq2SeqStat(self, prot_seq):
        #cnt = 0 
        for pid, _seq in list(prot_seq.items()): #pid in self.getCommSeqRegion(): 
            _seq = prot_seq[pid] 
            if pid in self.seq_stat:
                try:
                    self.seq_stat[pid][_seq] +=1
                except KeyError:
                    self.seq_stat[pid][_seq] = 1
            else:
                self.seq_stat[pid] = {_seq:1}
            #cnt += 1 
        #print "### cnt ", cnt 
        self.num_seq_in_stat += 1
    
    
    def statProtLengh(self):
        cnt = 0 
        allele_names = list(self.aln_seqs.keys())
        allele_names.sort()
        for allele_name in allele_names:
            if self.isChangedExpression(allele_name):
                continue 
            aln_seq = self.aln_seqs[allele_name]
            aln_seq_len =  len(aln_seq)
            if aln_seq_len < self.MAXLEN:
                print("  # ", aln_seq_len,  allele_name)
             
    
    def mergeContinuosSeg(self, sids):
        sids.sort()
        segs = []
        num_sid = len(sids)
        next_start_sid = sids[0]
        #print sids 
        for i in range(num_sid-1):
            if sids[i+1] - sids[i] != 1:
                end_sid = sids[i]
                seg_len = end_sid - next_start_sid
                segs.append((seg_len, next_start_sid, end_sid))
                next_start_sid = sids[i+1]
        last_seg_len = sids[-1] - next_start_sid
        segs.append((last_seg_len, next_start_sid,sids[-1]))
        #print segs 
        #if len(segs)==1:
        segs.sort(reverse=True)
        #print segs[0] 
        return segs[0] #if next_start_sid == sids[-1]:
    
                
                      
    
    def getCommSeqRegion(self): 
        return self.reg
    
    
    def findAlleleSeqDif(self, allele_name1, allele_name2, region=None):
        new_cmids = {}
        if region is None:
            reg = list(range(self.MAXLEN))
        else:
            _s, _e = region
            reg = list(range(_s, _e+1))
        for pid in reg: #range(self.MAXLEN):
            _seq1 = self.getRealPosSeq(allele_name1, pid)
            _seq2 = self.getRealPosSeq(allele_name2, pid)
            if _seq1 != _seq2:
                new_cmids[pid] = _seq1
        return new_cmids 
    
    def prnComonRefSeq(self):
        print("# HLA abc common reference sequence[25-205](181 aas)")
        print(("### common ref seq: %s " %  self.ref_aseq[25:206]))
    
    def getAlleleGroupProtein(self, grp):
        """ any allele with the same specific protein """
        for prot in self.allele_grp_uniq_proteins[grp]:
            break
        return prot 
    
    
    def getGeneGroupAllele(self, grp):
        allele_prot = self.getAlleleGroupProtein(grp) 
        return self.allele_prot_names[allele_prot][0]
    
     
    def addProt2SeqStat(self, prot):
        prot_seq = self.getProtCommSeq(prot)
        for pid in self.getCommSeqRegion(): 
            _seq = prot_seq[pid] 
            if pid in self.seq_stat:
                try:
                    self.seq_stat[pid][_seq] +=1
                except KeyError:
                    self.seq_stat[pid][_seq] = 1
            else:
                self.seq_stat[pid] = {_seq:1}
        self.num_seq_in_stat += 1
        
    
    def setCommonSeqReg(self, region):
        """ set common sequence region with ids in allele sequences """
        _s, _e = region
        self.reg = list(range(_s, _e+1))        
    
    def __setProtSeqStat(self):
        self.seq_stat = {}
        self.num_seq_in_stat = 0 
    
    
    def getCommSeq(self, cutoff=0.9, vb=0, num_min=None, COMMREG=True):
        """ 90% percent in all sequences """
        if num_min is None:
            num_min = self.num_seq_in_stat * 0.9 
        anchor_seq = {} # 100% in all seqs 
        major_seq = {}  #  > 90% except anchor_seq 
        am_stat = {}
        if not COMMREG:
            sids = list(self.seq_stat.keys())
        else:
            sids = self.getCommSeqRegion()
        for sid in sids:
            for nt, val in list(self.seq_stat[sid].items()): # in self.seq_stat[sid].keys():
                if val == self.num_seq_in_stat:
                    anchor_seq[sid] = nt 
                    am_stat[sid] = val,nt 
                    break
                elif val > num_min:
                    major_seq[sid] = nt
                    am_stat[sid] = val,nt 
                    break 
        if vb:
            print("# anchor_seq: ", len(anchor_seq),  anchor_seq)
            print("# major_seq: ", len(major_seq), major_seq) 
        return anchor_seq, major_seq, am_stat
    
    
    def addGrpSeqStat2Gene(self, grp,  gene):
        """ add group seqence statistics to gene statistics """
        grp_aseq, grp_mseq, am_stat  =  self.getCommSeq()
        self.hla_grp_num_prot[grp] = self.num_seq_in_stat
        self.hla_grp_amseq[grp] = grp_aseq, grp_mseq 
        _gene = self.hla_gene_am_stat[gene]
        for key in am_stat:
            if key in _gene:
                _gene_cnt , _gene_nt = _gene[key]
                _grp_cnt, _grp_nt  = am_stat[key]
                if _gene_nt == _grp_nt:
                    new_cnt = _gene_cnt + _grp_cnt 
                    _gene[key] = new_cnt, _gene_nt  
                else:
                    del[_gene[key]]
            else:
                _gene[key] = am_stat[key]
    
    
    def addGeneSeqStat2HLA(self, gene):
        _gene_stat = self.__getGeneSeqStat(gene)
        #print("# gene_sat:%s : %s" % (gene, _gene_stat) )
        for key in _gene_stat:
            if key in self.hla_am_stat:
                gene_cnt, gene_nt = _gene_stat[key]
                hla_cnt, hla_nt = self.hla_am_stat[key]
                if gene_nt == hla_nt:
                    self.hla_am_stat[key] = gene_cnt+hla_cnt, hla_nt
                else:
                    del self.hla_am_stat[key]
            else:
                self.hla_am_stat[key] = _gene_stat[key] 
            
            
    def __getGeneSeqStat(self, gene):
        return self.hla_gene_am_stat[gene]
 
 
    def getGeneNumProt(self, gene):
        return self.hla_gene_num_prot[gene] 
    
    
    def getAnchorMajorSeq(self, am_stat, num_prot, cutoff=0.9):
        aseq = {}
        mseq = {}
        num_cut = num_prot * cutoff
        for sid in list(am_stat.keys()):
            num_cnt, nt = am_stat[sid]
            if num_cnt == num_prot:
                aseq[sid] = nt 
            elif num_cnt > num_cut:
                mseq[sid] = nt 
        return aseq, mseq
    

    def setHLAAchorMajorSeq(self):
        self.hla_amseq = self.getAnchorMajorSeq(self.hla_am_stat, self.hla_num_prot)
        self.hla_gene_amseq = {}
        for gene in self.hla_genes:
            am_stat = self.__getGeneSeqStat(gene)
            num_prot = self.getGeneNumProt(gene)
            self.hla_gene_amseq[gene] = self.getAnchorMajorSeq(am_stat, num_prot)
             
              
    def getHLAAMSeq(self):
        return self.hla_amseq
    
    
    def getGeneAMseq(self, gene):
        return self.hla_gene_amseq[gene]
        
        
    def setGrpCommSeq(self, vb=0):
        self.hla_gene_am_stat = {}
        self.hla_am_stat = {}
        self.hla_gene_num_prot = {}
        self.hla_grp_num_prot = {}
        self.hla_num_prot = 0 
        self.hla_grp_amseq = {}
        self.hla_prot_seq = {}
        for gene in self.hla_genes:
            self.hla_gene_am_stat[gene] = {}
            gene_grps = self.getGeneGrps(gene) 
            gene_num_prot_seq = 0
            cnt =  0 
            for grp in gene_grps:
                self.__setProtSeqStat()
                prots = self.getGrpProts(grp)
                for _prot in prots:
                    self.addProt2SeqStat(_prot)
                gene_num_prot_seq +=  self.num_seq_in_stat 
                self.addGrpSeqStat2Gene(grp, gene)
                self.hla_gene_num_prot[gene] = gene_num_prot_seq 
            self.hla_num_prot += gene_num_prot_seq 
            self.addGeneSeqStat2HLA(gene)
            
            
    def getHLAGrpAmseq(self,grp):
        """ return HLA gene group's anchor sequence and major sequence """
        return self.hla_grp_amseq[grp]
    
    def getHLAProtAmseq(self, prot):
        """ protein's major sequence set empty for find the most possiblbe match,
            set anchor sequence = protein sequence, it only find identical match in sequence
        """
        return {},self.hla_prot_seq[prot]
    
    def loadHLAAMSeq(self):
        """ self.hla_gene_num_prot = {}
        self.hla_grp_num_prot = {}
        self.hla_num_prot = 0 
        
        self.hla_grp_amseq = {}
        self.hla_prot_seq = {}
        """
 
        fp = self.getAnchorMajorSeqFp()
        if self.isExistFile(fp) and not self.UPDATE:
            nums, amseqs = self.loadObj(fp, vb=1)
            self.hla_num_prot, self.hla_gene_num_prot, self.hla_grp_num_prot = nums
            self.hla_amseq, self.hla_gene_amseq, self.hla_grp_amseq, self.hla_prot_seq = amseqs
        else:
            self.setGrpCommSeq(vb=1)
            self.setHLAAchorMajorSeq()
            self.setHLAProtSeqs()
            amseqs = self.hla_amseq, self. hla_gene_amseq, self.hla_grp_amseq, self.hla_prot_seq 
            nums = self.hla_num_prot, self.hla_gene_num_prot, self.hla_grp_num_prot  
            self.dumpObj((nums,amseqs), fp, vb=1)
        self.prnHLAAMSeq() 
    
    
    def sortNamebyLastSetDigit(self, names, sep="*"):
        sort_names = []
        for name in names:
            tmp, sd = name.rsplit(sep, 1)
            sort_names.append((int(sd), name))
        sort_names.sort() #sorted(sort_names)
        return [ y for x,y in sort_names]
    
            
    def prnHLAAMSeq(self):
        lns = ["# HLA proteins: total_num= {:8d}".format(self.hla_num_prot)]
        for gene in self.hla_genes:
            gene_grps = self.getGeneGrps(gene)
            aseq, mseq = self.getGeneAMseq(gene)
            lns.append("# gene:{}  num_prot={:6d}  num_grp: {:3d}  num_anchor:{:4d} num_major={:4d}".format(gene, \
                    self.hla_gene_num_prot[gene], len(gene_grps),len(aseq), len(mseq)))
            cnt = 0 
            for grp in gene_grps: #= self.getGeneGrps(gene)
                cnt += 1
                aseq, mseq = self.getHLAGrpAmseq(grp)
                lns.append("# {:3d}  grp: {:6s}  num_prot:{:5d} num_res: anchor: {:3d} major: {:3d}".format(cnt, grp, self.hla_grp_num_prot[grp], \
                                                                                                       len(aseq), len(mseq)))
        #self.prnLines(lns)
        self.add2log(lns)
        
        
    def testGrpProtStat(self, grp):
        """
        # error in grp: A*36 , no major seqs
# error in grp: A*43 , no major seqs
# error in grp: A*69 , no major seqs
# error in grp: A*80 , no major seqs
# error in grp: B*81 , no major seqs
# error in grp: B*73 , no major seqs
# error in grp: B*59 , no major seqs
# error in grp: B*82 , no major seqs
# error in grp: B*83 , no major seqs
# error in grp: B*78 , no major seqs
# error in grp: B*67 , no major seqs
# error in grp: B*47 , no major seqs

        """
        self.__setProtSeqStat()
        prots = self.getGrpProts(grp)
        for _prot in prots:
            self.addProt2SeqStat(_prot)
        grp_aseq, grp_mseq, am_stat  =  self.getCommSeq()
        print(( "# grp prot stat:  aseq: %s" % grp_aseq ))
        print(( "# mseq: ", grp_mseq ))
        print(( "# am_stat: ", am_stat)) 
    
    #def setnonCommProtSeq(self):
    #    for prot in self.
    
    def protSeqStat(self, prots, reg=None):
        self.__setProtSeqStat()
        for _prot in prots:
            self.addProt2SeqStat(_prot)
        grp_aseq, grp_mseq, am_stat  =  self.getCommSeq()
       
        
    def findHLADupProtSeq(self):
        """ find duplicated in protein sequences in HLA """
        prot_aseq_str = {}
        #prots = self.getGrpProts("C*14")
        prots = self.getHLAProts()
        print(("# num_prot_in_HLA: {} ".format(len(prots))))
        for prot in prots:
        #prot = "C*14:57"
            #self.statProtAlleleLength(prot)
            aseq = self.statHLAUniqProtSeq(prot)
            aseq_str = self.aseq2Str(aseq)
            if aseq_str in prot_aseq_str:
                print((" prot: {}  ==  {} len: {} ".format(prot, prot_aseq_str[aseq_str],len(aseq_str))))
            else:
                prot_aseq_str[aseq_str] = prot 
    
    def setHLAProtSeqs(self):
        """ extracted protein sequences for  same protein in HLA """
        prots = self.getHLAProts()
        print(("# num_prot_in_HLA: {} ".format(len(prots))))
        self.hla_prot_seq = {}
        for prot in prots:
            aseq = self.statHLAUniqProtSeq(prot)
            self.hla_prot_seq[prot] = aseq 
            
        
    def statHLAAllele(self):
        self.loadAlignedSeq()
        self.setHLAAlignedSeqDt()
        self.__setAlleleGrps() 
        self.setRefSeq()
        #self.prnComonRefSeq()
        ##self.statProtLengh()
        #self.UPDATE = True
        self.loadHLAAMSeq()
        #self.setUniqProtSeq()
        #self.hla_prot_seq = {}
        #self.testGrpProtStat("A*36")
        #self.test()
        #self.findHLADupProtSeq()
        
    def alignPDBSeq2RefSeq(self, seq, MAX_NUM_RES=275):
        full_ref_seq = self.getRefSeq()
        start_id = 22 
        _rseq = full_ref_seq[start_id:305]
        #print(ref_seq[start_id:305])
        #print(seq) 
        alns = pairwise2.align.globalds(_rseq, seq, matrix, gap_open, gap_extend)
        aligned_seq = alns[0][1]
        self.add2log("# ref_seq: {}".format(alns[0][0]))
        self.add2log("# aln_seq: {}".format(aligned_seq))
        aseq = {}
        cnt = 0 
        for i in range(len(aligned_seq)):
            if aligned_seq[i] != "-":
                aseq[i+start_id] = aligned_seq[i]
                cnt += 1
                if cnt > MAX_NUM_RES:
                    break 
        #print(aseq)
        self.add2log("# num_res_in_aligned_seqdt: {}".format(len(aseq)))
        return aseq
    
    def prnCommRegion(self):
        creg = self.getCommSeqRegion()
        print("# MHCI common region: {} -{}".format(creg[0], creg[-1]))
    
    
    def getCommRegSeqdt(self, aseq_dt):
        com_aseq_dt = {}
        for i in self.getCommSeqRegion():
            if i in aseq_dt:
                com_aseq_dt[i] = aseq_dt[i]
        return com_aseq_dt
    