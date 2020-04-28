import os, sys, pickle, shutil
from string import digits
from re import split
from math import  log, sin, cos, sqrt, pi, fabs, radians,degrees, acos, asin, atan2, pi, exp, ceil 
from sys import exit, platform
from socket import gethostname
from time import strftime, sleep 
from copy import copy 

import numpy as Npy
_array = Npy.array
def err_exit(cmt):
    print("\n# **********   Error  Message  ************\n#")
    print(cmt) 
    print("#")
    print("# **********   Error  Exit  ************\n")
    exit(0)

def dumpObj(obj, fp, **rest):
    """ dump obj into file fn on disk"""
    fmode = rest.get("fmode", 'wb')
    if 'dirpath' in rest: 
        fp = os.path.join(rest['dirpath'],fp)
    elif 'dpath' in rest:
        fp = os.path.join(rest['dpath'],fp)
    cmt = rest.get('cmt', "")
    msg = '# dump < %s > to file: %s' % (cmt, fp )
    fh = open(fp,fmode)                
    pickle.dump(obj,fh,1)
    fh.close()
    vb = rest.get('vb',1)
    if vb:  print(msg)
    return

def loadObj(fp, **rest):
    if 'dpath' in rest: 
        fp = os.path.join(rest['dpath'],fp)
    elif 'dirpath' in rest:
        fp = os.path.join(rest['dirpath'],fp)
    vb=rest.get('vb', 1)
    if vb: print('# load file: %s' % fp) 
    fh = open(fp, 'rb')
    _obj =  pickle.load(fh)
    fh.close()
    return _obj

class ERRMSG:
      
    def err_exit(self, err_msg):
        print("\n# **********   Error  Message  ************\n#")
        print(err_msg) 
        print("#")
        print("# **********   Error  Exit  ************\n")
        exit(0)
    
    def war_msg(self, war_msg):
        print("# Warning: %s " % war_msg) 
    
    def warning(self, _msg):
        print("#-------------------------------------")
        print("# Warning: ")
        print("#          %s" %  _msg)
        print("# ")
        print("# ------------------------------------")
        
        
class DICTEXTRA:
    def compDict(self, dt1, dt2):
        comm_keys = []
        dif_keys1 = []
        for key1 in list(dt1.keys()):
            if key1 in dt2:
                comm_keys.append(key1)
            else:
                dif_keys1.append(key1)
        dif_keys2 = self.getNewItems(comm_keys, list(dt2.keys()))
        return comm_keys, dif_keys1, dif_keys2
                
    
    def getNewItems(self, sub_li, dst_li):
        new_items = []
        for _item in dst_li:
            if _item not in sub_li:
                new_items.append(_item)
        return new_items

    
class MIO(ERRMSG, DICTEXTRA):
    log = []
    vb = 0
    def __init__(self, args, options=None):
        self.args = args
        self.options = options 
        self.cwd = None
        self.setHostname()
        self.setCWD()
    
    def setArgs(self, args=None):
        if args is not None:
            self.args = args 
            if hasattr(args, "cmd_lns"):
                self.prnLines(args.cmd_lns)
        
        else:
            if not hasattr(self, "args"):
                self.err_exit("no args")
        argv_dt = args.__dict__
        for name in list(argv_dt.keys()):
            setattr(self, name, argv_dt[name])
        
        self.setCWD()
        self.setHostname()
        
        
    def getCommandLines(self, cmt="", vb=0):
        cmd = " ".join(sys.argv)
        cmd_lns = [ "\n#---------------------------"]
        cmd_lns += [ "#Runing: python %s" % cmd, "#host: %s     #working dir: %s" % (self.getHostname(), self.getCWD()),  "# %s" % cmt ]
        cmd_lns += ["#---------------------------\n"]        
        self.cmd_lns = cmd_lns 
        if vb: self.prnLines(cmd_lns)
        #print self.hostname
        return cmd_lns
    
        
        
    def setHostname(self):
        self.hostname = gethostname()
    
    def getHostname(self):
        if not hasattr(self, "hostname"):
            self.setHostname()
            return self.hostname
        
    def getFileSize(self, fp):
        return os.path.getsize(fp)
  
        
    def setArgsSrcDir(self, args):
        if not hasattr(self, "cwd"):
            self.setCWD()
        if args.src_dir is None:
            self.src_dir  = self.getCwd()
        else:
            self.src_dir = args.src_dir
        self.args = args  
        

    def isWindows(self):
        if platform != "linux2":
            return True
        else:
            return False
        
    def dumpObj(self, obj, fp, **rest):
        dumpObj(obj, fp, **rest)
        
    def loadObj(self, fp, **rest):
        return loadObj(fp, **rest)
        
        
    def addSeperator(self, lns):
        lns.append("# -------------------------------------")
        return lns 
    
    def getDateStr(self):
        """ return date in fomrat of Mon-Day-Year"""
        return strftime("%m_%d_%Y")
    
    
    def add2log(self, lns, vb=0, prn=0):
        if type(lns) != type([]):
            lns = [lns]
        if not hasattr(self, 'log'):
            self.log = []
        self.log += lns 
        if prn: self.prnLns(lns)
        if vb or self.vb: 
            self.prnLines(["\n## adding to log: " ] + lns + ["###  end of adding to log\n#********"])
    
    def prnLog(self):
        self.prnLines(self.log)
        
    def saveLog(self, fn_log, **kwargs):
        self.saveLines(self.log, fn_log, **kwargs)
        
    def prnLines(self, lines, index=None, **rest  ):
        #print "# prnlines: ", rest
        if not rest.get('prn', 1): 
            #print("# Warning: no prn=1 in the function prnLines(..,prn=1), no print")            
            return
        if not rest.get("vb", 1):
            return 
        cmt= rest.get("cmt", "")
        if cmt: print(cmt) 
        num = len(lines)
        chkrow = rest.get("chkrow", 1)
        line_thresh = rest.get("thresh", 400)
        _omit = False
        if num > line_thresh  and chkrow:
            num = 100
            _omit = True            
        if index:
            for i in range(num): print('# %3d %s ' % ( i + 1, rstrip(lines[i])))
        else:
            for i in range(num): print(lines[i].rstrip())
        if _omit: print("# lines no >%s, full txt, use < chkrow=0 >" % line_thresh)
    
    
    def prnDict(self, _dict, sort=False, cmt="", raw=0):
        if cmt: print(cmt) 
        dkeys = list(_dict.keys())
        if sort:
            dkeys.sort()
        if raw:
            print("{")
            for _key in dkeys[:-1]:
                print("    %s:%s," % (_key, _dict[_key]))
            last_key = dkeys[-1]
            print("    %s:%s" % (last_key, _dict[last_key]))
            print("}") 
        else:
            for _key in dkeys:
                _dict[_key].sort()
                print("# %s -> %s "  % (_key, _dict[_key]))

    
    def rmBrace(self, ln):
        ln = ln.replace("(", " ")
        ln = ln.replace(")", " ")
        return ln 
    

    def prnLns(self, lns, **rest):
        self.prnLines(lns, **rest)
    
    
    def readTxtFile(self, fp, vb=0):
        if vb: print("# reading fp: ", fp) 
        fh = open(fp,'r')
        txt = fh.readlines()
        #print "# test num_ln: %6d" %  len(txt)
        fh.close()
        return txt 
    
    def txt2col(self, ccs):    
        """ only works for block data """
        all_col_strs = [] 
        for line in ccs:
            line = strip(line)
            if line and line[0] != "#":
                col_strs = split(r'\s+',line)
                all_col_strs.append(col_strs)
        return all_col_strs 
    
                
    def txt2colDat(self, ccs, dat_type="arr"):    
        """ only works for block data """
        bdata = []  # block data (2 dimensions's list
        for line in ccs:
            line = strip(line)
            if line and line[0] != "#":
                col_strs = split(r'\s+',line)
                col_ln = [ float(_item) for _item in col_strs]
                bdata.append(col_ln)
        if dat_type == 'list':
            return bdata
        elif dat_type == 'arr':
            return npy.array(bdata)
        else:
            print("#  output data type: %s " % dat_type) 
            raise "# no matched data type(list, arr) "

    
    def saveLines(self, lines, fp,  fmode='w', vb=0, **rest):
        """ save lines list into a disk file """
        cmt = rest.get("cmt", "lines")
        if vb: print("##Save %s into file: %s " % (cmt, fp)) 
        
        fh = open(fp, fmode)
        for line in lines:
            line = line.strip() 
            fh.write(line + '\n')
        fh.close()
        prn = rest.get("prn", False)
        if prn: self.prnLines(lines)
        
 
  
    def prnVec(self, vec, cmt= " "):
        print("# vec %s : %s " % (cmt, self.getVecStr(vec)))
    
    def getWorkInf(self):
        lns = ["# work information: run time: %s " % self.getDateStr() ]
        lns +=["# cwd: %s " % os.getcwd()]
        lns +=["# python %s " % " ".join(self.argv)]
        lns +=["# -----------------------"]
        return lns 
    
      
    def intli(self, li):
        new_li = []
        for item in li:
            new_li.append(int(item))
        return new_li
    
    
    def prnList(self, var_list, index=0):
        if index:
            for i in range(len(var_list)):
                print("i= %2d %s " % (i, var_list[i]))
        else:
            for i in range(len(var_list)):
                print("%s " % (var_list[i]))    
    

    def getUnMatchElem(self, atoms1, atoms2):
        """  input: list atoms1, list atoms2
            return: unmatched elements in list atoms1 """
        um_atoms1 = []
        for _atom in atoms1:
            if _atom not in atoms2:
                um_atoms1.append(_atom)
        return um_atoms1        
    
    
    def getMatchElem(self, atoms1, atoms2):
        """  input: list atoms1, list atoms2
            return: unmatched elements in list atoms1 """
        um_atoms1 = []
        for _atom in atoms1:
            if _atom in atoms2:
                um_atoms1.append(_atom)
        return um_atoms1    

    def getBestMatchElem(self, atoms1, atoms2):
        if len(atoms1) < len(atoms2):
            return self.getMatchElem(atoms1, atoms2)
        else:
            return self.getMatchElem(atoms2, atoms1)
            
    
    def getUniqKey(self, atoms):
        atoms.sort()
        return self.getKey(atoms)
    
    def getKey(self, atoms):
        return "_".join(atoms) #atoms.join("_")
    
    def getRevKey(self, atoms):
        atoms.reverse()
        return self.getKey(atoms)
    
    def getRevDict(self, _dt):
        rev_dt = {}
        for key, val in list(_dt.items()):
            rev_dt[val] = key 
        return rev_dt
    
        

    def isExistFile(self, fn):
        if os.path.exists(fn):
            return True
        else:
            return False
        
    def isExist(self, fp):
        if os.path.exists(fp):
            return True
        else:
            return False
        
    
    def isExistDir(self,fp):
        if self.isExist(fp) and os.path.isdir(fp):
            return True
        else:
            return False
        
        
    def copyFile(self, src, dst, vb=0):
        
        shutil.copyfile(src, dst)
        if vb:
            print(("# copied file: %s -> %s" % (src, dst)))

    def getFnPre(self, fn, vb=0):
        """ return return filename without extension"""
        if '.' in fn:
            fn_pre, fn_ext = fn.rsplit(".", 1)
        else:
            fn_pre = fn 
            fn_ext = 'NA'
        if vb: print("# fn: %s  fn_pre: %s  fn_ext: %s" %  (fn, fn_pre, fn_ext))
        return fn_pre
    
    
    def getFnExt(self, fn, vb=0):
        """ return return filename without extension"""
        if '.' in fn:
            #tmp = fn.rsplit(".", 1)
            #print "# tmp: ", tmp 
            fn_pre, fn_ext = fn.rsplit(".", 1)
        else:
            fn_pre = fn 
            fn_ext = ''
        if vb: print("# fn: %s  fn_pre: %s  fn_ext: %s" %  (fn, fn_pre, fn_ext))
        return fn_pre, fn_ext    

   
    def getFnScreeOut(self, _fn):
        fn_pre = self.getFnPre(_fn)
        return "%s.sco"  %fn_pre        
   
       
    def backupFile(self, fn):
        if not self.isExistFile(fn):
            return 
        fn_pre, fn_ext = self.getFnExt(fn)
        for i in range(100):
            if fn_ext:
                new_fn = "%s_%s.%s" % (fn_pre, i, fn_ext)
            else:
                new_fn = "%s_%s" % (fn_pre, i)
            if not self.isExistFile(new_fn):
                break 
        #self.copyFile(fn, new_fn)
        print("# move file: %s -> %s " % (fn, new_fn))
        shutil.move(fn, new_fn)
        
        
                
    def getCWD(self):
        return os.getcwd()
    
    def getCwd(self):
        return os.getcwd()
    
    
    def joinPath(self, path, fn):
        return os.path.join(path, fn)
    
    
    def join3Path(self, path1, path2, fn):
        return os.path.join(path1, path2, fn)
    
    
    def getDirFiles(self, _dir, **rest):
        #print _dir
        files = os.listdir(_dir)
        fns = []
        #fn_pre = rest.get("fn_pre", None)
        #fn_suf = rest.get("fn_suf", None)
        #print "# rest: ", rest , "fn_pre: ", fn_pre, "fn_suf:", fn_suf
        for fn in files:
            #if fn_pre is not None and not fn.startswith(fn_pre): continue
            #if fn_suf is not None and  not fn.endswith(fn_suf): continue
            #print "# ffn: ", fn 
            if os.path.isfile(self.joinPath(_dir,fn)):
                fns.append(fn)
        vb = rest.get("vb", 0)
        #print fns 
        #if vb: 
        #    print "# return fns: ", fns 
        if rest:
            #print "# rest: ", rest 
            return self.filtFns(fns, **rest)
        else:
            return fns 
    

    def getDirs(self, _dir, **rest):
        files = os.listdir(_dir)
        fns = []
        for fn in files:
            if os.path.isdir(self.joinPath(_dir,fn)):
                fns.append(fn)
        vb = rest.get("vb", 0)
        if rest:
            return self.filtFns(fns, **rest)
        else:
            return fns 
    
        
    def getBaseFn(self, fn):
        if "/" in fn:
            fpath, fn = fn.rsplit("/", 1)
        return fn 
    
    
    def filtFns(self, fns, fn_pre=None, fn_ext=None,fn_mid=None,  vb=0):
        #print "# filter: fn_pre: %s  fn_ext: %s" % (fn_pre, fn_ext)
        new_fns = []
        for fn in fns:
            fn = self.getBaseFn(fn)
            #print fn 
            if fn_pre is not None and  not fn.startswith(fn_pre): continue
            #print fn 
            if fn_ext is not None and  not fn.endswith(fn_ext): continue   
            #print "# found fn: %s" % fn 
            if fn_mid is not None and fn_mid not in fn:
                continue
            new_fns.append(fn)
        return new_fns
    
    def setCWD(self):
        if not hasattr(self, "cwd") or self.cwd is None:
            self.cwd = self.getCWD()
            self.cwdn = os.path.basename(self.cwd)        
        
        
    def setDestDir(self, dst_dir=None, vb=0):
        if not hasattr(self, "cwd"): self.setCWD()
        if dst_dir is not None:
            pass
        elif hasattr(self, 'options'): 
            dst_dir = self.options.dstdir
        if dst_dir is None:
            if hasattr(self, "args") and len(self.args) == 1:
                dst_dir = self.args[0]
                dst_path = self.joinPath(self.cwd, dst_dir)
                if self.isExistDir(dst_path):
                    self.dst_path = dst_path
                else:
                    self.dst_path = self.cwd 
            else:
                self.dst_path = self.cwd
                self.dst = os.path.basename(self.cwd)
        else:
            self.dst_path = self.joinPath(self.cwd, dst_dir)  
        self.dst_path = self.rmlastSlash(self.dst_path)
        self.dstn = os.path.basename(self.dst_path)    
        self.chkDir(self.dst_path)
        if vb:
            print("# cwd: %s" % self.cwd) 
            print("# dst_path: %s" % self.dst_path)
            print("# dst_dirn: %s" % self.dstn)
            
        
    def setDestDirFns(self, **rest):
        self.fns_dstdir = self.getDirFiles(self.dst_path)
        vb = rest.get("vb", 0)
        if vb:
            print("# dst_path: %s" % self.dst_path)
            print("# self.fns_dstdir: ", self.fns_dstdir)
    
    def getDestDirFns(self, **rest):
        if not hasattr(self, 'fns_dstdir'):
            self.setDestDirFns()
        return self.filtFns(self.fns_dstdir, **rest)
    
    
    def rmlastSlash(self, fn):
        if fn is None: return 
        if fn[-1] == "/":
            return fn[:-1]
        else:
            return fn 
        
        
    def setSrcDir(self, src_dir=None, vb=0):
        if not hasattr(self, "cwd"):  self.setCWD()
        if hasattr(self, 'options'):  src_dir = self.options.srcdir
        
        if src_dir is None:
            self.src_path = self.cwd
            self.srcn = self.cwdn
        else:
            src_dir = self.rmlastSlash(src_dir)
            self.src_path = self.joinPath(self.cwd, src_dir)  
            if not self.isExistDir(self.src_path):
                self.err_exit("Src path not existed: %s" % self.src_path)
        self.src_dirn = os.path.basename(self.src_path)    
        if vb:
            print("# cwd: %s" % self.cwd) 
            print("# src_path: %s" % self.src_path)
            print("# src_dirn: %s" % self.src_dirn)        

        
    def setSrcDirFns(self, **rest):
        self.fns_src = self.getDirFiles(self.src_path)
        vb = rest.get("vb", 0)
        if vb:
            print("# src_path: %s" % self.src_path)
            print("# self.fns_src: ", self.fns_src)     
            

    def getSrcDirFns(self, **rest):
        if not hasattr(self, 'fns_src'):
            self.setSrcDirFns()
        return self.filtFns(self.fns_src, **rest)    

    
    def getDestFp(self, fn):
        if fn is None: return 
        return self.joinPath(self.dst_path, fn)
    
    
    def getSrcFp(self, fn):
        if fn is None: return
        if hasattr(self, "src_path"):
            return self.joinPath(self.src_path, fn)
        else:
            return self.joinPath(self.dst_path, fn)    
    

    def getBaseName(self, fn):
        return os.path.basename(fn)
    
    
    def getEqSplit(self, ln):
        return split("\s*=\s*", ln)
    
    def getFileTime(self, fp):
        return os.path.getmtime(fp)
        
    
    def cpFiles(self, fns):
        cnt = 1
        for src_fn, fn in fns:
            print("# copy %4d: %s -> %s" % (cnt, src_fn, fn ))
            cnt += 1
        ans = input("Are you sure to add these files? [Y/N]: " )
        if ans in ["Y","y"]:
            for src_fn, fn in fns:
                shutil.copy(src_fn, fn)
                #pass         
                
    
    def getPDBFns(self, fns, fn_pre=None, vb=0):
        pdb_fns = []
        for fn in fns:
            #print fn 
            if fn.endswith(".pdb"):
                if fn_pre is not None:
                    if fn.startswith(fn_pre):
                        pdb_fns.append(fn)
                else:
                    pdb_fns.append(fn)
        self.num_pdb = len(pdb_fns)
        print("# num_pdb= %4d " % self.num_pdb)
        if vb:
            print("# pdb_fns: ", pdb_fns)
        return pdb_fns
    
    
    def getOptionAttr(self, attr):
        if hasattr(self, "options") and hasattr(self.options, attr):
            return getattr(self.options, attr)
            
    
    def setFnPre(self, fn_pre=None, vb=0):
        if fn_pre is None:
            self.fn_pre = self.getOptionAttr("fn_pre")
        else:
            self.fn_pre = fn_pre 
        if vb:
            print("# fn_pre: ", self.fn_pre) 
        
    def mvFiles(self, fps, dst, test=1, vb=0):    
        for fp in fps:
            fn = os.path.basename(fp)
            dst_fp = os.path.join(dst, fn)
            if not test:
                if vb: print("# move: %s --> %s " % (fp, dst_fp))
                shutil.move(fp, dst_fp)
            else:
                print("# Testing, not move: %s --> %s " % (fp, dst_fp))

    def mvFile2Dst(self, old_new_fps):
        for fp_old, fp_new in old_new_fps:
            shutil.move(fp_old, fp_new)
            
    def mvFile(self, src_fp, dst_fp, vb=0):
        if vb: print("# move: %s --> %s " % (src_fp, dst_fp))
        shutil.move(src_fp, dst_fp)        
         
            
    def chkDir(self, dst, vb=1):
        if self.isExist(dst):
            if vb: print("#*** Warning: already existed: %s " % dst)
        else:
            os.makedirs(dst)
        
      
    def wait(self, tm):
        tval, tunit = tm
        if tunit =="h":
            tval *= 3600
        elif tunit == "m":
            tval *= 60
        print("# wait %s (secs) " % tval)
        sleep(tval)
                    
                    
    def getStdTxtOut(self, _std_out, name, num_sample=0):
        """ return standard deviation analysis results """
        _av, _std, _min, _max = _std_out
        ln = "# %s(av)= %7.2f %s(std)= %7.3f  %s_min= %7.3f %s_max= %7.3f" % (name, _av, name, _std, name, _min,  name, _max)
        if num_sample > 0:
            ln += " num_sample= %4d" % num_sample
        return ln
        
    
    def fn2pdbid(self, fn):
        fn = os.path.basename(fn)
        fn_pre = self.getFnPre(fn)        
        if fn.endswith(".pdb"):
            #print "# **********"
            if len(fn_pre) == 4:
                pdbid = fn_pre
                return pdbid
            elif "_" in fn_pre:
                print(fn_pre) 
                if  fn_pre.endswith("_m2"):
                    pdbid, tt = fn_pre.split("_")
                elif fn_pre.endswith("_wH"):
                    t1, pdbid, tt = fn_pre.split("_")
                elif fn_pre.startswith("mhcII") or fn_pre.startswith("aligned_") or fn_pre.startswith("am_"):
                    tt, pdbid = fn_pre.split("_")
                else:
                    print("# Error: could not parse pdbid from fn: %s" % fn) 
                    return 
                assert len(pdbid) == 4
                return pdbid
            else:
                print("#Error: non parsed pdbid: %s" % fn_pre)
        elif fn.endswith(".ent") and fn.startswith("pdb") and len(fn_pre) == 7:
            pdbid = fn_pre[3:]
            #print fn_pre , pdbid
            assert len(pdbid) == 4
            return pdbid
        else:
            print("#Error: could not parse PDBID from filename: %s  " % fn) 

    def chkCmd(self, cmd):
        if self.hostname.startswith("KU"): 
            print("# in windows, no checking for cmd: %s" % cmd) 
            return 
        cmd = cmd.strip()
        elem = cmd.split()
        if len(elem) > 1:
            command = elem[0]
        else:
            command = cmd 
        if not self.isExist(command):
            self.err_exit("# command not exsited: %s" %  command)
         
    def notExistPath(self, _path, vb=0, cmt=""):
        if not self.isExistDir(_path):
            if vb:
                print("#* Warning: %s Path NOT  existed: %s" % (cmt, _path))                 
            return True
        else:
            return False 
    
    def notExistFp(self, _fp, vb=0, cmt=""):
        if not self.isExistFile(_fp):
            if vb:
                print("#* Warning: %s File NOT  existed: %s" % (cmt, _fp))                 
            return True
        else:
            return False 

 
    def fp2FpTxt(self, fp):
        return "{}.txt".format(fp[:-4])
    
    def fp2FpBin(self, fp):
   
        return "{}.bin".format(fp[:-4])
   
    def fp2FpLog(self, fp):
        return "{}.log".format(fp[:-4])
    
    
    def isNew(self, fp):
        if self.isExistFile(fp):
            if hasattr(self, "UPDATE") and self.UPDATE:
                return True
            else:
                return False
        else:
            return True
        
               
                                    
if (__name__ == "__main__"):
    pass 