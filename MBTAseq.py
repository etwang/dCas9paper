import os, sys, subprocess, gzip
from pylab import *

def splitBarcodesSeparateFiles(filelist_f, plasmidcodes_f, fastqdir, out_f):

    files = []
    for line in open(filelist_f):
        files.append(line.strip())
    
    i = 0 
    codeToPlasmid = {} 
    plasmids = []
    for line in open(plasmidcodes_f):
        plasmid, code = line.strip().split()    
        codeToPlasmid[code] = i
        plasmids.append(plasmid) 
        i += 1

    mat = zeros((len(files), len(plasmids)), dtype='int')

    fcode = 0
    for f in files:
        print f
        fh = gzip.open(os.path.join(fastqdir, f))
        line = fh.readline()
        i = 0
        while len(line) > 0:
            if line.startswith("@"):
                seq = fh.readline() 
                fh.readline()
                fh.readline()
                line = fh.readline()
                if i%100000 == 0:
                    print mat
                pos = seq.find("AACGTGGATTGGGGTTGTTG")
                scode = seq[pos - 8 : pos]
                if scode in codeToPlasmid:
                    plasmid = codeToPlasmid[scode]
                    mat[fcode, plasmid] += 1
                i += 1
        fcode += 1

    out = open(out_f, 'w')
    out.write("#Samples\t" + "\t".join(plasmids) + "\n")
    for i in range(mat.shape[0]):
        out.write(files[i] + "\t" + "\t".join(map(str, mat[i, :])) + "\n")
    out.close()








def splitBarcodes(filecodes_f, plasmidcodes_f, fastq_f, out_f):
    
    codeToF = {}
    i = 0
    files = []
    for line in open(filecodes_f):
        f, code = line.strip().split()
        files.append(f)
        codeToF[code] = i
        i += 1 
  
    i = 0 
    codeToPlasmid = {} 
    plasmids = []
    for line in open(plasmidcodes_f):
        plasmid, code = line.strip().split()    
        codeToPlasmid[code] = i
        plasmids.append(plasmid) 
        i += 1

    mat = zeros((len(files), len(plasmids)), dtype='int')
    fh = gzip.open(fastq_f)
    line = fh.readline()
    i = 0
    while len(line) > 0:
        if line.startswith("@"):
            fcode = line.strip().split(":")[-1]
            seq = fh.readline() 
            fh.readline()
            fh.readline()
            line = fh.readline()
            if i%100000 == 0:
                print mat
            if fcode in codeToF:
                f = codeToF[fcode]
                pos = seq.find("AACGTGGATTGGGGTTGTTG")
                scode = seq[pos - 8 : pos]
                if scode in codeToPlasmid:
                    plasmid = codeToPlasmid[scode]
                    mat[f, plasmid] += 1
            i += 1

    out = open(out_f, 'w')
    out.write("#Samples\t" + "\t".join(plasmids) + "\n")
    for i in range(mat.shape[0]):
        out.write(files[i] + "\t" + "\t".join(map(str, mat[i, :])) + "\n")
    out.close()





