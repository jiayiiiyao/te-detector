from __future__ import print_function

import bisect
import signal
import sys
import optparse
import operator
import signal

class Alignment:
    def __init__(self, refChr, refStart, refEnd, refLen, refSeq,
                 qryName, qryStart, qryEnd, qryLen, qrySeq):
        self.refChr, self.refStart, self.refEnd, self.refLen, self.refSeq = refChr, refStart, refEnd, refLen, refSeq
        self.qryName, self.qryStart, self.qryEnd, self.qryLen, self.qrySeq = qryName, qryStart, qryEnd, qryLen, qrySeq

class Transposon:
    def __init__(self, genoName, genoStart, genoEnd, strand,
                  repName, repClass, repFamily):
        self.genoName, self.genoStart, self.genoEnd, self.strand = genoName, genoStart, genoEnd, strand
        self.repName, self.repClass, self.repFamily = repName, repClass, repFamily

class Insertion:
	def __init__(self, targetChr, targetSite, readName, insertionStart, insertionEnd, TSD):
		self.targetChr, self.targetSite, self.readName, self.insertionStart, self.insertionEnd, self.TSD = targetChr, targetSite, readName, insertionStart, insertionEnd, TSD
	def getInfo(self):
		data = self.targetChr, self.targetSite, self.readName,self.insertionStart, self.insertionEnd, self.TSD.length, self.TSD.seq, self.TSD.left_flanking, self.TSD.right_flanking
		return data

class TSD:
	def __init__(self, length, seq, left_flanking, right_flanking):
		self.length, self.seq, self.left_flanking, self.right_flanking = length, seq, left_flanking, right_flanking

def openAndLog(fileName):
    f = open(fileName)
    sys.stderr.write("reading " + fileName + "..." + "\n")
    return f

def writeAndLog(fileName):
    f = open(fileName, 'w')
    sys.stderr.write("writing " + fileName + "..." + "\n")
    return f

# For reading Alignments
def readBlockFromMaf(lines):
    block = []
    for line in lines:
        if line.isspace():
            if block:
                yield block
                block = []
        elif line[0] != '#':
            block.append(line)
    if block:
        yield block

def parseAlignments(fields):
    name = fields[1]
    start, lens = int(fields[2]), int(fields[3])
    seq = fields[6]
    end = start + len(seq) - seq.count('-')
    return name, start, end, lens, seq

def readAlignments(filename):
    alignments = {}
    blocks = readBlockFromMaf(openAndLog(filename))
    for block in blocks:
        refname, refstart, refend, reflen, refseq = parseAlignments(block[1].split())
        qryname, qrystart, qryend, qrylen, qryseq = parseAlignments(block[2].split())
        align = Alignment(refname, refstart, refend, reflen, refseq, qryname, qrystart, qryend, qrylen, qryseq)
        alist = alignments.get(refname)
        if alist == None:
            alist = []
        alist.append(align)
        alignments[refname] = alist
    return alignments

# read repeats
def readRepeats(filename):
    repeats = {}
    file = openAndLog(filename)
    for line in file:
        attributes = line.split()
        genoName = attributes[5]
        genoStart, genoEnd = int(attributes[6]), int(attributes[7])
        strand, repName, repClass, repFamily= attributes[9:13]
        if not repFamily == ' ':
            transposon = Transposon(genoName, genoStart, genoEnd, strand, repName,
                            repClass, repFamily)
            tlist = repeats.get(genoName)
            if tlist == None:
                tlist = []
            tlist.append(transposon)
            repeats[genoName] = tlist
    file.close()
    for chr in repeats.keys():
       repeats[chr].sort(key = operator.attrgetter('genoStart'))
    return repeats


# For bisecting 
def getRefStartPoints(alignlist):
    alist = []
    for align in alignlist:
        alist.append(align.refStart)
    return alist