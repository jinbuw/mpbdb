from .libalnallele import HLAAllele

class Node(object):
    #nid = 0 
    def __init__(self, name, parent):
        self.name= name
        self.parent = parent 
        self.children={}
        self.nid  = None
        self.num_children = 0 
        self.DB =  False
        self.prob_cut = 0.9
        self.num_poss_match = 0 
        self.MAXPC = 20
        
    
    def isTermNode(self):
        if self.num_children == 0:
            return True
        else:
            return False
        
    def setNodeID(self, nid):
        self.nid = nid
        
    def setCommSeq(self, cseq):
        self.cseq = cseq
        self.cseq_len = len(cseq)
        
    def setPrevSib(self, node):
        self.prev_sib = node 
        
    def setAMSeq(self, amseq, vb=0):
        """ set anchor and major sequences for this node """
        self.anchor_seq, self.major_seq = amseq
        self.asids = list(self.anchor_seq.keys())
        self.msids = list(self.major_seq.keys())
        self.asids_len = float(len(self.asids))
        self.msids_len = float(len(self.msids))
        if vb: print(("# anchor_sids: {}  major_sids: {}".format(self.asids, self.msids)))
     
    
    def getPossMChildren(self, seq, vb=0):
        temp = []
        for _name  in self.children: #i in range(self.num_children):
            _node = self.children[_name]
            #print("# node_name: {}".format(_node.name))
            simi = self.children[_name].seq2Similarity(seq)
            if self.DB or vb: 
                print(("# child_name: %4s  simi: %6.2f" % (_name, simi)))
            if self.children[_name].isTermNode():
                if simi > self.prob_cut:
                    temp.append((simi, _name))
            else: 
                if simi > 0:
                    temp.append((simi, _name))
        temp.sort(reverse=True)
        return temp
        
    def getMatchChild(self, seq):
        pchildren = self.getPossMChildren(seq)
        results = []
        for _prob, name in pchildren:
            if not self.children[name].isTermNode():
                results.append((_prob,name))
                results += self.children[name].getMatchChild(seq) 
            else:
                if _prob > 0.95:
                    results.append((_prob,name))
        return results 
    
      
    def addChild(self, child):
        """ child node """
        self.children[child.name] = child
        self.num_children += 1 
        
    
    def __getAnchorSim(self, seq):
        cnt = 0
        num_res = 0
        for sid in self.asids:
            if sid in seq:
                num_res += 1 
                if seq[sid] == self.anchor_seq[sid]:
                    cnt += 1
        sim = cnt/num_res
        if self.DB: print(("#anchor: cnt: {}  len:{}  sim:{}".format(cnt, num_res,sim)))
        return sim 
    
    
    def __getMajorSim(self, seq):
        cnt = 0
        num_res = 0
        for sid in self.msids:
            if sid in seq:
                num_res += 1
                if seq[sid] == self.major_seq[sid]:
                    cnt += 1
        sim = cnt/num_res
        if self.DB: print(("#major: cnt: {}  len:{}  sim:{}".format(cnt, num_res,sim)))
        return sim
        
    
    def seq2Similarity(self, seq):
        if self.asids_len > 0: # has achor amino acids 
            asim = self.__getAnchorSim(seq)
            if asim < 0.99:
                return -1 #0.0
            else:
                if self.msids_len > 0:
                    return self.__getMajorSim(seq)
                else:
                    return asim
        else:
            return self.__getMajorSim(seq)
        
    def prnSeqSim(self, seq):
        sim = self.seq2Similarity(seq)
        print("# node: {} sim: {:3f}".format(self.name, sim))
        
        
    def getNodeID(self):
        return self.nid 
    
    def __str__(self):
        return "#nid= %4d name= %8s  children= %2d" % (self.nid, self.name, self.num_children)
    
    def getNodeName(self):
        return self.name

class Tree:
    def __init__(self):
        """ nid: node index"""
        self.initTree()
        
        
    def initTree(self):
        self.node_id = 0 
        self.nodes = {}
        self.name_nids = {}

    def get_index(self, position):
        for index, node in enumerate(self.nodes):
            if node.identifier == position:
                break
        return index
    
    
    def getNode(self, node_name=None, nid=None):
        return self.nodes[node_name] #pass 
        

    def createNode(self, name, identifier=None, parent=None):
        node = Node(name, identifier)
        self.node_id += 1
        node.setNodeID(self.node_id)
        #print(node)
        name = node.getNodeName()
        if name in self.nodes:
            print(("#*** Error duplicated node name in the tree: {}".format(name)))
        else:
            self.nodes[name] = node 
        return node

    def setRoot(self, node):
        self.root = node 
        
    #def __getitem__(self, key):
        #return self.nodes[self.get_index(key)]

    #def __setitem__(self, key, item):
        #self.nodes[self.get_index(key)] = item

    #def __len__(self):
        #return len(self.nodes)
    def prnTreeNodes(self):
        print("#num_nodes: %4d" % len(self.nodes))
        for nid in list(self.nodes.keys()):
            print(self.nodes[nid]) 

class HLATree(Tree, HLAAllele):
    dup_prots = ["C*07:722","A*02:768","C*03:460", "C*15:176","C*12:228"]
    def isDupProtSeq(self, prot):
        if prot in self.dup_prots:
            print(("# dupseq prot: {} excluded".format(prot)))
            return True
        else:
            return False
                        
    def createDTree(self):
        """ create HLA decision tree """
        self.initTree()
        root = self.createNode("HlA")
        root.setAMSeq(self.getHLAAMSeq())
        self.setRoot(root)
        for gene in self.hla_genes: 
            gene_node = self.createNode(gene, parent=root)
            gene_node.setAMSeq(self.getGeneAMseq(gene))
            root.addChild(gene_node)
            grps = self.getGeneGrps(gene)
            for grp in grps:
                grp_node = self.createNode(grp, parent=gene_node)
                aseq, mseq = self.getHLAGrpAmseq(grp)
                if len(mseq) == 0 and len(aseq)== 0 :
                    print("# error in grp: %s , no anchor and  major seqs" % grp) 
                grp_node.setAMSeq(self.getHLAGrpAmseq(grp))
                gene_node.addChild(grp_node) 
                for prot in self.getGrpProts(grp):
                    if self.isDupProtSeq(prot):
                        continue
                    prot_node =  self.createNode(prot, parent=grp_node)
                    grp_node.addChild(prot_node)
                    prot_node.setAMSeq(self.getHLAProtAmseq(prot))
        #self.prnTreeNodes()
        
        
    def testSeq2Allel(self):
        num_test = 0
        num_corr = 0 
        for gene in self.hla_genes:
            grps = self.getGeneGrps(gene) # for  in GroupAllele(grp)
            #grp = grps[0]
            for grp in grps:
                prot, grp_seq = self.getGrpSeqProt(grp)
                _name, prob = self.findMatchProt(grp_seq)
                num_test += 1
                if prot.startswith(_name):
                    num_corr += 1
                else:
                    print(("#Error in math: prot:{} -> {}".format(prot,_name)))
                    
                print(("# grp: {} prot: {} m: {} prob: {:.3f}".format(grp,prot,_name, prob)))
        print(("# num_test: {}, num_corr: {}  corr_ration: {:.3f} ".format(num_test, num_corr, num_corr/float(num_test))))
    
    def getProtSeq(self, prot):
        return self.hla_prot_seq[prot]
    
    def testMatchProt(self, prot):
        prot_seq = self.getProtSeq(prot)
        match_names, mns  = self.findMatchProt(prot_seq)
        #prob, _name = pn 
        print(("# match: len: {} prot:{} -> {}".format(len(prot_seq), prot,match_names)))
        if prot in match_names: #.startswith(_name):
            if len(match_names) > 1:
                print(("# multi match: {}".format(match_names)))
            return True
        else:
            print(("#Error in math: prot:{} -> {}".format(prot,_name)))
            print(("# mnodes:", mns))
            print() 
            return False
    
    def getNodeSeqSim(self, node_name, aseq):
        node = self.getNode(node_name=node_name)
        node.DB = 1
        #children = node.getPossMChildren(aseq)
        sim = node.seq2Similarity(aseq)
        print("# node_name: {}, prob: {:5.3f}".format(node_name,sim))
    
    def chkNodeChildren(self, node_name, aseq):
        node = self.getNode(node_name=node_name)
        node.DB = 1
        #print(node.anchor_seq)
        #print(aseq[124])
        children = node.getPossMChildren(aseq)
                    
    def chkNodeAnchor(self, node_name, aseq):
        node = self.getNode(node_name=node_name)
        node.DB = 1
        print("# node: {}  anchor res: ".format(node_name))
        print(node.anchor_seq)
        #print(aseq[124])
        node.prnSeqSim(aseq) #seq2Similarity(aseq)
        print(node.anchor_seq, len(node.anchor_seq))
        node.asids.sort() 
        print(node.asids)
        seq_dif = self.getSeq2SeqDiff(node.anchor_seq, aseq)
        print("# anchor_seqdif: ", seq_dif)
        #children = node.getPossMChildren(aseq)
        
    def matchProtAseq(self, aseq_dt):
        match_names, mns  = self.findMatchProt(aseq_dt)
        #print(("# match: len: {}  -> {}".format(len(aseq_dt), match_names)))
        #print(mns)
        if len(mns) == 0:
            return [],[],[]
        max_prot = mns[0][0]
        prob_cut =max_prot * 0.95
        max_num = 10
        cnt = 0 
        mdiffs = []
        for prob, prot_name in mns:
            if prob > prob_cut:
                cnt += 1 
                self.add2log("# prob: {:5.3f}  prot: {:8s}".format(prob, prot_name))
                seq_diff = self.getSeqDiff(prot_name, aseq_dt)
                if seq_diff:
                    ln = "#seq_diff:"
                    for sid, prot_nt, seq_nt in seq_diff:
                        ln += "  {}:{}/{} ".format(sid, prot_nt, seq_nt) 
                    #print(ln) 
                    self.add2log(ln)
                mdiffs.append((prob, prot_name, seq_diff))
                if cnt > max_num: break  
        return match_names, mns, mdiffs
    
    
    def getSeqDiff(self, prot, aseq):
        """ return sequence diff: [(sid, prot_seq, aseq)"""
        prot_seq = self.getProtSeq(prot)
        return self.getSeq2SeqDiff(prot_seq, aseq) 
    
    def getSeq2SeqDiff(self, prot_seq, aseq):
        """ return sequence diff: [(sid, prot_seq, aseq)"""
        seq_dif = []
        cnt = 0 
        #print(aseq)
        for sid in aseq:
            if sid in prot_seq:
                cnt += 1
                if aseq[sid] != prot_seq[sid]:
                    seq_dif.append((sid, prot_seq[sid], aseq[sid]))
            #else:
            #    print("# sid not in protseq: {}".format(sid))
        #print("# joint num_res: {}".format(cnt))
        return seq_dif
    
    
    def testGrpMatch(self, grp):
        num_test = 0
        num_corr = 0 
        for prot in self.getGrpProts(grp):
            prot_seq = self.getProtSeq(prot)
            #print("# protin: {}".format(prot))
            #print self.findMatchProt(prot_seq)
            pn, mns  = self.findMatchProt(prot_seq)
            prob, _name = pn 
            num_test += 1
            if prot.startswith(_name):
                num_corr += 1
            else:
                print(("#Error in math: prot:{} -> {}".format(prot,_name)))
                print(("# mnodes:", mns))
                print() 
            #break 
        print(("# grp: {} num_prot: {}   num_corr: {} Ratio:{:.3f}".format(grp, num_test, num_corr, num_corr/float(num_test))))
    
    def testHLAProtMatch(self):
        num_test = 0
        num_corr = 0 
        hla_prots = self.getHLAProts()
        num_hla_prots = len(hla_prots)
        print(("# test all proteins in HLA: num_prot: {}".format(num_hla_prots)))
        fn = "test_corr_prot.bin"
        if self.isExistFile(fn):
            test_corr_prots = self.loadObj(fn)
        else:
            test_corr_prots = []
        num_save = 0    
        for prot in hla_prots:
            if prot in test_corr_prots:
                 continue 
            if self.isDupProtSeq(prot):
                continue
            prot_seq = self.getProtSeq(prot)
            #print("# protin: {}".format(prot))
            match_names, mns  = self.findMatchProt(prot_seq)
            #prob, _name = pn 
            num_test += 1
            if prot in match_names: #.startswith(_name):
                num_corr += 1
                #print("prot: {}  -> {}".format(prot,_name))
                if len(match_names) > 1:
                    print(("# multi match: {}".format(match_names)))                
                test_corr_prots.append(prot)
                num_save += 1
                if num_save == 50:
                    self.dumpObj(test_corr_prots, fn, vb=0)
                    num_save = 0
            else:
                print(("#Error in math: prot:{} -> {}".format(prot,match_names)))
                #print("# mnodes:", mns)
                print() 
            #break 
        print(("#  num_prot: {}   num_corr: {} Ratio:{:.3f}".format( num_test, num_corr, num_corr/float(num_test))))
        if num_save > 0:
            self.dumpObj(test_corr_prots, fn, vb=0)
                 
        
    def findMatchProt(self, seq):
        #self.root.seq2Similarity(seq)
        match_nodes = self.root.getMatchChild(seq)
        match_prot = []
        for _prob, name in match_nodes:
            if ":" in name:
                match_prot.append((_prob, name))
        match_prot.sort(reverse=True)
        #min_num = min(10, len(match_prot))
        match_names = []
        for _prob, name in match_prot:
            if _prob > 0.999:
                match_names.append(name)
        return match_names, match_prot #match_prot[:min_num]
        