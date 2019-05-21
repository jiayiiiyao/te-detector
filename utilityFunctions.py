from __future__ import print_function
from itertools import groupby
import bisect
import logging
import collections
import signal
import sys
import optparse
import operator


class Alignment:
    aligncount = 0
    def __init__(self, refChr, refStart, refEnd, refLen, refSeq,
                 qryName, qryStart, qryEnd, qryLen, qrySeq):
        self.refChr, self.refStart, self.refEnd, self.refLen, self.refSeq = refChr, refStart, refEnd, refLen, refSeq
        self.qryName, self.qryStart, self.qryEnd, self.qryLen, self.qrySeq = qryName, qryStart, qryEnd, qryLen, qrySeq
        Alignment.aligncount += 1

class Transposon:
    def __init__(self, genoName, genoStart, genoEnd, strand,
                  repName, repClass, repFamily):
        self.genoName, self.genoStart, self.genoEnd, self.strand = genoName, genoStart, genoEnd, strand
        self.repName, self.repClass, self.repFamily = repName, repClass, repFamily

class Insertion:
	def __init__(self, targetChr, targetSite, readName, insertionStart, insertionLength, TSD):
		self.targetChr, self.targetSite, self.readName, self.insertionStart, self.insertionLength, self.TSD = targetChr, targetSite, readName, insertionStart, insertionLength, TSD
	def getInfo(self):
        	data = self.targetChr, self.targetSite, self.readName, self.insertionStart, self.insertionLength
        	return data

class TSD:
	def __init__(self, length, seq, left_flanking, right_flanking):
		self.length, self.seq, self.left_flanking, self.right_flanking = length, seq, left_flanking, right_flanking

def openAndLog(fileName):
    logging.info("opening " + fileName + "..." + "\n")
    return open(fileName)

def writeAndLog(fileName):
    logging.info("writing " + fileName + "..." + "\n")
    return open(fileName, 'w')

def logReads(f, results):
    for r in results:
        print(r, end='\n', file=f)

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

def readAlignments(f):
    alignments = collections.defaultdict(list)
    blocks = readBlockFromMaf(f)
    for block in blocks:
        refname, refstart, refend, reflen, refseq = parseAlignments(block[1].split())
        qryname, qrystart, qryend, qrylen, qryseq = parseAlignments(block[2].split())
        align = Alignment(refname, refstart, refend, reflen, refseq, qryname, qrystart, qryend, qrylen, qryseq)
        alignments[qryname].append(align)
    return alignments

def readRepeats(f):
    repeats = collections.defaultdict(list)
    for line in f:
        attributes = line.split()
        genoName = attributes[5]
        genoStart, genoEnd = int(attributes[6]), int(attributes[7])
        strand, repName, repClass, repFamily = attributes[9:13]
        if not repFamily == 'Simple_repeat':
            transposon = Transposon(genoName, genoStart, genoEnd, strand, repName, repClass, repFamily)
            repeats[genoName].append(transposon)
    for genoname in repeats.keys():
       repeats[genoname].sort(key = operator.attrgetter('genoStart'))
    return repeats


def logNegative(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align, insert = r[0:2]
        target_chr, target_site, donor_chr, donor_start, donor_end = insert.targetChr, insert.targetSite, align.refChr, align.refStart, align.refEnd
        data = donor_chr, donor_start, donor_end, align.qryName, align.qryStart, align.qryEnd
        print(*data, sep='\t', end='\n', file=outfile)
    outfile.close()

def logMapped(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align, insert = r[0:2]
        data = align.refChr, align.refStart, align.refEnd, align.qryName, align.qryStart, align.qryEnd, insert.targetChr, insert.insertionStart, insert.insertionLength
        print(*data, sep='\t', end='\n', file=outfile)
    outfile.close()

def logPositive(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align, te, insert, transduction = r[0:4]
        tsd = insert.TSD
        target_chr, target_site, donor_chr, donor_start, donor_end = insert.targetChr, insert.targetSite, align.refChr, align.refStart, align.refEnd
        data1 = target_chr, target_site, donor_chr, donor_start, donor_end, align.qryName, align.qryStart, align.qryEnd, insert.insertionStart, insert.insertionLength
        data2 = te.repName, te.repClass, te.repFamily, te.genoStart, te.genoEnd
        data3 = tsd.length, tsd.seq, tsd.left_flanking, tsd.right_flanking, transduction[0], transduction[1]
        print(*data1, sep='\t', end='\n', file=outfile)
        print(*data2, sep='\t', end='\n', file=outfile)
        print(*data3, sep='\t', end='\n', file=outfile)
        outfile.write('\n')
    outfile.close()

def logUnmapped(filename, results):
    outfile = writeAndLog(filename)
    for r in results:
        data = r.getInfo()
        print(*data, sep='\t', end='\n', file=outfile)   

def logInsertions(results, f):
    for name in results.keys():
        rlist = results[name]
        for r in rlist:
            data = r.getInfo()
            print(*data, sep='\t', end='\n', file=f)
