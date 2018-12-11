from __future__ import print_function

import bisect
import signal
import sys
import optparse
import operator
import signal

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
    for name in alignments.keys():
        alignments[name].sort(key=operator.attrgetter('refStart'))
    return alignments

# For bisecting 
def getRefStartPoints(alignlist):
    startpoints = []
    for align in alignlist:
        startpoints.append(align.refStart)
    return startpoints

def getAdjacent(alist, end):
    candidates = []
    tsd = 15 # target site duplication here 
    startpoints = getRefStartPoints(alist)
    pos = bisect.bisect_right(startpoints, end)
    if pos != len(startpoints) and pos >= 2:
        for next in alist[pos - 1], alist[pos]:
            if abs(next.refStart - end) <= tsd:
                candidates.append(next)
    return candidates

def getAdjacent_2(alist, end):
    tsd = 30
    for a in alist:
        if abs(a.refStart - end) <= tsd:
            return a
    return None

def checkGap(align1, align2):
    bias = 100
    if align1.qryName == align2.qryName:
        if (align2.qryStart - align1.qryEnd) > bias:
            return True
    return False

# return list of Insertion entities
def findInsertion(alignments, outfile):
    sys.stderr.write("Start finding insertions..." + "\n")
    insertionList = []
    for chr in alignments.keys():
        alist = alignments[chr]
        for align1 in alist:
            end = align1.refEnd
            align2 = getAdjacent(alist, end)
            if align2:
                isgap = checkGap(align1, align2)
                if isgap:
                    length, seq, left, right = checkTSD(align1, align2)
                    tsd = TSD(length, seq, left, right)
                    insertion = Insertion(align1.refChr, align1.refEnd, align1.qryName, align1.qryEnd, align2.qryStart, tsd)
                    #insertionList.append(insertion)
                    logSingleInsertion(outfile, insertion)
    #return insertionList
    return

def logInsertionReadName(alignments):
    readlist = []
    for chr in alignments.keys():
        alist = alignments[chr]
        for align1 in alist:
            end = align1.refEnd
            candidates = getAdjacent(alist, end)
            if len(candidates) > 0:
                for align2 in candidates:
                   # sys.stderr.write("Start checking candidates..." + "\n")
                    isgap = checkGap(align1, align2)
                    if isgap:
                        readlist.append(align1.qryName)
    return readlist

def findInsertion_2(alignments):
    sys.stderr.write("Start finding insertions..." + "\n")
    insertionList = []
    for chr in alignments.keys():
        alist = alignments[chr]
        for align1 in alist:
            end = align1.refEnd
            align2 = getAdjacent_2(alist, end)
            if align2:
                isgap = checkGap(align1, align2)
                if isgap:
                    length, seq, left, right = checkTSD(align1, align2)
                    tsd = TSD(length, seq, left, right)
                    insertion = Insertion(align1.refChr, align1.refEnd, align1.qryName, align1.qryEnd, align2.qryStart, tsd)
                    insertionList.append(insertion)
    return insertionList

def checkTSD(align1, align2):
    tsd_length = align1.refEnd - align2.refStart
    if tsd_length > 0:
        tsd_seq = align2.refSeq[0 : tsd_length]
        tsd_left = align1.refSeq[-1 * (tsd_length + 10) : -1 * tsd_length]
        tsd_right = align2.refSeq[tsd_length : tsd_length + 10]
    else:
        tsd_seq, tsd_left, tsd_right = "N", "N", "N"
    return tsd_length, tsd_seq, tsd_left, tsd_right

def logInsertion(filename, results):
    outfile = writeAndLog(filename)
    for r in results:
        data = r.getInfo()
        print(*data, sep='\t', end='\n', file=outfile)


def logReads(outfile, reads):
    f = open(outfile, 'a+')
    for r in reads:
        print(r, end='\n', file=f)
    f.close()


def tedet_insertions(args):
    alignments = readAlignments(args[0])
    sys.stderr.write("alignmentCount " + str(Alignment.aligncount) + "\n")
    #insertions = findInsertion(alignments, args[1])
    #findInsertion(alignments, args[1])
    #insertions = findInsertion_2(alignments)
    #logInsertion(args[1], insertions)
    outfile = args[1]
    readnames = logInsertionReadName(alignments)
    logReads(outfile, readnames)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] maf-file"
    description = "Find insertions in reads."
    op = optparse.OptionParser(description=description)
    #op.add_option("-v", "--verbose", action="count", help="show progress messages")
    args = op.parse_args()[1]
    #sys.stderr.write("args " + args[1] + "..." + "\n")
    # if len(args) < 2:
    #     op.error("please give me alignments in MAF format and Output Directory")
    tedet_insertions(args)




