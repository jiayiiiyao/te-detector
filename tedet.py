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

    def getTransInfo(self):
        data = self.repName, self.repClass, self.repFamily, self.genoName, \
               self.genoStart, self.genoEnd
        return data

def openAndLog(fileName):
    f = open(fileName)
    sys.stderr.write("reading" + fileName + "..." + "\n")
    return f

def writeAndLog(fileName):
    f = open(filename, 'w')
    sys.stderr.write("writing" + fileName + "..." + "\n")
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

def readAlignments(filename):
    alignments = {}
    blocks = readBlockFromMaf(openAndLog(filename))
    for block in blocks:
        refname, refstart, refend, reflen, refseq = parseAlignments(block[1].split())
        #if refname.startswith("chrUn") == False:
        qryname, qrystart, qryend, qrylen, qryseq = parseAlignments(block[2].split())
        align = Alignment(refname, refstart, refend, reflen, refseq, qryname, qrystart, qryend, qrylen, qryseq)
        alist = alignments.get(refname)
        if alist == None:
            alist = []
        alist.append(align)
        alignments[refname] = alist
    return alignments

def parseAlignments(fields):
    name = fields[1]
    start, lens = int(fields[2]), int(fields[3])
    seq = fields[6]
    end = start + len(seq) - seq.count('-')
    return name, start, end, lens, seq

def readRepeats(filename):
    repeats = {}
    file = openAndLog(filename)
    for line in file:
        attributes = line.split()
        genoName = attributes[5]
        genoStart, genoEnd = int(attributes[6]), int(attributes[7])
        strand, repName, repClass, repFamily, repStart, repEnd = attributes[9:15]
        transposon = Transposon(genoName, genoStart, genoEnd, strand, repName,
                        repClass, repFamily, repStart, repEnd)
        tlist = repeats.get(genoName)
        if tlist == None:
            tlist = []
        tlist.append(transposon)
        repeats[genoName] = tlist
    file.close()
    for chr in repeats.keys():
       repeats[chr].sort(key = operator.attrgetter('genoStart'))
    return repeats

def checkInRmsk(repeats, chr, start, end):
    bias = 20
    startpoints = []
    telist = repeats[chr]
    for te in telist:
        startpoints.append(te.  genoStart)
    size = len(startpoints)
    pos = bisect.bisect_right(startpoints, start)
    if pos != size:
        te = telist[pos]
        if te.repFamily != 'Simple_repeat':
            tend, tstart = te.genoEnd, te.genoStart
            if end - tend <= bias and start - tstart <= bias:
                return te
    else:
        return None

def getRefStartPoints(alignlist):
    alist = []
    for align in alignlist:
        alist.append(align.refStart)
    return alist


def findAnnotated(alignments, repeats):
    result = []
    for alignlist in alignments.values():
        for a in alignlist:
            refchr, refstart, refend = a.refChr, a.refStart, a.refEnd
            te = checkInRmsk(repeats, refchr, refstart, refend)
            if te:
                reference = refchr, refstart, a.refSeq
                query = a.myName, a.myStart, a.mySeq
                repeat = te.repName, te.repClass, te.repFamily
                found = reference, query, repeat
                print("foundAnnotated")
                result.append(found)
    return result

def getAdjacent(alist, end):
    bias = 15
    startpoints = getRefStartPoints(alist)
    pos = bisect.bisect_right(startpoints, end)
    if pos != len(startpoints) and pos >= 2:
        for next in alist[pos - 2:pos]:
            if abs(next.refStart - end) < bias:
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
        alist.sort(key=operator.attrgetter('refStart'))
        for align1 in alist:
            end = align1.refEnd
            align2 = getAdjacent(alist, end)
            if align2:
                isgap = checkGap(align1, align2)
                if isgap:
                    insertion = align1, align2
                    insertionList.append(insertion)
    return insertionList

def checkIfMapped(insertions, alignments):
    unmapped = []
    mapped = []
    for insertion in insertions:
        candidates = []
        bias = 50
        #align1 = insertion[0], align2 = insertion[1]
        name, start, end = insertion[0].myName, insertion[0].myEnd, insertion[1].myStart
        for alignmentlist in alignments.values():
            for aligment in alignmentlist:
                if aligment.myName == name and abs(start - aligment.myStart) < bias and abs(end - aligment.myEnd) < bias:
                    candidates.append(aligment)
        if len(candidates) == 0:
            gap = name, start, end
            unmapped.append(gap)
        else:
            mappedData = insertion, candidates
            mapped.append(mappedData)
    return unmapped, mapped

def checkMappedInRmsk(mapped, repeats):
    TEpositive = []
    TEnegative = []
    for map in mapped:
        candidates = map[1]
        for c in candidates:
            chr, start, end = c.refChr, c.refStart, c.refEnd
            te = checkInRmsk(repeats, chr, start, end)
            if te:
                positive = map[0], te
                TEpositive.append(positive)
            else:
                TEnegative.append(map)
    return TEpositive, TEnegative

def logGap(outfile, gap):
    print(*gap, sep='\t', end='\n', file=outfile)

def logAlignment(outfile, align):
    refdata = align.refChr, align.refStart, align.refEnd
    mydata = align.myName, align.myStart, align.myEnd
    print(*refdata, sep='\t', end='\n', file=outfile)
    print(*mydata, sep='\t', end='\n', file=outfile)

def logTE(outfile, te):
    tedata = te.genoName, te.genoStart, te.genoEnd, te.repName, te.repClass, te.repFamily
    print(*tedata, sep='\t', end='\n', file=outfile)


def logResult(filename, result, option):
    outfile = writeAndLog(filename)
    if option == 'unmapped':
        for r in result:
            logGap(outfile, gap)
    elif option == 'TEmapped':
        for r in result:
            align1, align2 = r[0][0:2]
            te = r[1]
            logAlignment(outfile, align1)
            logAlignment(outfile, align2)
            logTE(outfile, te)
            outfile.write('\n')
    elif option == 'TEunmapped':
        for r in result:
            align1, align2 = r[0][0:2]
            logAlignment(outfile, align1)
            logAlignment(outfile, align2)
            outfile.write('\n')
    elif option == 'statistics':
        print(*result, sep='\n', end='\n', file=outfile)
    outfile.close()


def findTE(args):
    repeats = readRepeats(args[0])
    alignments = readAlignments(args[1])
    alignment_count = Alignment.aligncount

    insertions = findInsertion(alignments)

    unmapped, mapped = checkIfMapped(insertions, alignments)

    outputfile = args[2]
    logResult(outputfile + "_unmapped", unmapped, 'unmapped')
    TEpositve, TEnegative= checkMappedInRmsk(mapped, repeats)
    logResult(outputfile_dir + "_TEmapped", TEpositve, 'TEmapped')
    logResult(outputfile_dir + "_TEunmapped", TEnegative, 'TEunmapped')

    statistics = Alignment.aligncount, len(insertions), len(unmapped), len(mapped), len(TEpositve), len(TEnegative)
    logResult(outputfile_dir + "_statistics", statistics, 'statistics')


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "repeats.txt alignments.maf outputFile"
    description = "Try to find transposable elements insertion."
    op = optparse.OptionParser (description=description)
    #op.add_option("-d", "--detail", help="Show detailed TE sequence")
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("please give me repeats, MAF alignments and outputfilename")
    findTE(args)




