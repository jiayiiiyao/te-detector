from __future__ import print_function

import optparse
import operator
import bisect
import signal
import sys

class Alignment:
    aligncount = 0

    def __init__(self, refChr, refStart, refEnd, refLen, refSeq,
                 myName, myStart, myEnd, myLen, mySeq):

        self.refChr, self.refStart, self.refEnd, self.refLen, self.refSeq = refChr, refStart, refEnd, refLen, refSeq
        self.myName, self.myStart, self.myEnd, self.myLen, self.mySeq= myName, myStart, myEnd, myLen, mySeq
        Alignment.aligncount += 1


class Transposon:
    def __init__(self, genoName, genoStart, genoEnd, strand,
                  repName, repClass, repFamily, repStart, repEnd):

        self.genoName, self.genoStart, self.genoEnd, self.strand = genoName, genoStart, genoEnd, strand
        self.repName, self.repClass, self.repFamily, self.repStart, self.repEnd = repName, repClass, repFamily, repStart, repEnd


def openAndLog(fileName):
    f = open(fileName)
    sys.stderr.write("reading " + fileName + "..." + "\n")
    return f

def writeAndLog(fileName):
    f = open(fileName, 'w')
    sys.stderr.write("writing " + fileName + "..." + "\n")
    return f

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

def getRefStartPoints(alignlist):
    alist = []
    for align in alignlist:
        alist.append(align.refStart)
    return alist


def getAdjacent(alist, end):
    tsd = 15 # target site duplication here 
    startpoints = getRefStartPoints(alist)
    pos = bisect.bisect_right(startpoints, end)
    if pos != len(startpoints) and pos >= 2:
        for next in alist[pos - 1], alist[pos]:
            if abs(next.refStart - end) <= tsd:
                return next
    return None

def checkGap(align1, align2):
    bias = 200
    myname1, myname2 = align1.myName, align2.myName
    if myname1 == myname2:
        end = align1.myEnd
        start = align2.myStart
        #if abs(start - end) > bias:
        if (start - end) > bias:
            return True
    return False

def findInsertion(alignments):
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
                    insertion = align1, align2
                    insertionList.append(insertion)
    return insertionList

def logInsertion(filename, results):
    outfile = writeAndLog(filename)
    for r in results:
        align1, align2 = r[0], r[1]
        acceptor = align1.refChr, align1.refEnd, align1.refStart
        read = align1.myName, align1.myEnd, align2.myStart
        if align1.refEnd > align2.refStart:
            tsd_length = align1.refEnd - align2.refStart
            reverse = -1*tsd_length
            tsd_1 = align1.mySeq[reverse:0]
            tsd_2 = align2.mySeq[0:tsd_length]
        else:
            tsd_length = 0
            tsd_1, tsd_2 = " ", " "
        tsd = tsd_length, tsd_1, tsd_2
        print(*acceptor, sep='\t', end='\n', file=outfile)
        print(*read, sep='\t', end='\n', file=outfile)
        print(*tsd, sep='\t', end='\n', file=outfile)
        outfile.write('\n')


def tedet_insertions(args):
    alignments = readAlignments(args[0])
    insertions = findInsertion(alignments)
    outputfile = args[1]
    logInsertion(outputfile, insertions)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "alignments.maf outputfile"
    description = "Try to find insertions."
    op = optparse.OptionParser (description=description)
    args = op.parse_args()[1]
    if len(args) < 2:
        op.error("please give me MAF alignments and outputfile")
    tedet_insertions(args)




