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

# read/write file

def openAndLog(fileName):
    f = open(fileName)
    sys.stderr.write("reading " + fileName + "..." + "\n")
    return f

def writeAndLog(fileName):
    f = open(fileName, 'w')
    sys.stderr.write("writing " + fileName + "..." + "\n")
    return f

# read Alignments

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

# read repeats
def readRepeats(filename):
    repeats = {}
    file = openAndLog(filename)
    for line in file:
        attributes = line.split()
        genoName = attributes[5]
        genoStart, genoEnd = int(attributes[6]), int(attributes[7])
        strand, repName, repClass, repFamily, repStart, repEnd = attributes[9:15]
 #       if repFamily != 'Simple_repeat' and repFamily != 'Low_complexity':
        if not (repFamily == 'Simple_repeat' or 'Low_complexity'):
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

# check if annotated

def getRefStartPoints(alignlist):
    alist = []
    for align in alignlist:
        alist.append(align.refStart)
    return alist

def getCandidateTE(telist, refstart):
    startpoints = []
    tes = []
    for te in telist:
        startpoints.append(te.genoStart)
    size = len(startpoints)
    pos = bisect.bisect_right(startpoints, refstart)
    for i in range(pos - 5, pos + 5):
        if i >= 0 and i < size:
            tes.append(telist[i])
    return tes

def checkCandidates_1(pairs, repeats):
    bias = 100
    positive_full = []
    positive_partial = []
    negative = []
    for pair in pairs:
        a, insertion = pair[0], pair[1]
        refchr = a.refChr.split('_')[0]
        refstart, refend, readname, readstart, readend = a.refStart, a.refEnd, a.myName, int(a.myStart), int(a.myEnd)
        if refchr in repeats.keys(): 
            cans = getCandidateTE(repeats[refchr], refstart)
            for te in cans:
                testart, teend = te.genoStart, te.genoEnd
                if refstart <= testart and refend >= teend:
                    transduction = checkTransduction(a, te)
                    res = a, te, insertion, transduction
                    positive_full.append(res)
                    continue
                elif refstart >= testart and refend <= teend:
                    transduction = "none"
                    res = a, te, transduction, insertion
                    positive_partial.append(res)
                    continue
                elif (testart - refstart) <= bias and teend <= refend:
                    transduction = checkTransduction(a, te)
                    res = a, te, insertion, transduction
                    positive_partial.append(res)
                    continue       
                elif testart > refstart and (teend - refend) <= bias:
                    transduction = "none"
                    res = a, te, transduction, insertion
                    positive_partial.append(res)
                    continue
                else:
                    negative.append(pair)
                    continue
        else:
            negative.append(pair)
    return positive_full, positive_partial, negative

def checkCandidates(pairs, repeats):
    bias = 100
    positive_full = []
    positive_partial = []
    negative = []
    for pair in pairs:
        a, insertion = pair[0], pair[1]
        #insertionstart -> donarstart
        insertionstart, insertionend = int(insertion[2]), int(insertion[3])
        donarstart = a.refStart + (insertionstart - a.myStart)
        donarend = a.refEnd - (a.myEnd - insertionend)
        refchr = a.refChr.split('_')[0]
        #refstart, refend, readname, readstart, readend = a.refStar, a.refEnd, a.myName, int(a.myStart), int(a.myEnd)
        if refchr in repeats.keys(): 
            cans = getCandidateTE(repeats[refchr], donarstart)
            #cans = repeats[refchr]
            for te in cans:
                testart, teend = te.genoStart, te.genoEnd
                if donarstart <= testart and donarend >= teend:
                    transduction = checkTransduction(a, donarend, te)
                    res = a, te, insertion, transduction
                    positive_full.append(res)
                    break
                elif (testart <= donarstart and teend >= donarstart and teend <= donarend) or (testart >= donarstart and testart <= donarend and teend >= donarend) or (testart <= donarstart and teend >= donarend):
                    transduction = checkTransduction(a, donarend, te)
                   # transduction = "none"
                    res = a, te, insertion, transduction
                    positive_partial.append(res)
                    break       
                else:
                    res = a, insertion
                    negative.append(res)
                    break
        else:
            res = a, insertion
            negative.append(res)
    return positive_full, positive_partial, negative


# insertions: dictionary
# accpetor_chr, acceptor_site, readstart, readend, tsdlength, tsdseq
def readInsertions(insertionfile):
    insertions = {}
    lines = openAndLog(insertionfile)
    for line in lines:
        fields = line.split()
        readname = fields[3]
        ilist = insertions.get(readname)
        if ilist == None:
            ilist = []
        if len(fields) == 9:
            data = fields[0], fields[1], fields[4], fields[5], fields[6], fields[7], fields[8]
        else:
            data = fields[0], fields[1], fields[4], fields[5], fields[6]
        ilist.append(data)          
        insertions[readname] = ilist
    return insertions

def selectedFromAlignments(alignments, insertions):
    selected = []
    for alignlist in alignments.values():
        for align in alignlist:
            readname, readstart, readend = align.myName, align.myStart, align.myEnd
            if readname in insertions.keys():
                ilist = insertions[readname]
                for i in ilist:
                    insertstart, insertend = int(i[2]), int(i[3])
                    if not (readstart >= insertend or readend <= insertstart) :
                        pair = align, i
                        selected.append(pair)
    return selected


def checkTransduction(align, donarend, te):
    threshold = 20
    length = donarend - te.genoEnd
    if length >= threshold:
        startseq, endseq = te.genoEnd - align.refStart, donarend - align.refStart
        transduction = te.genoEnd, donarend, length, align.mySeq[startseq:endseq]
    else:
        transduction = "None" 
    return transduction


def logNegative(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
       # sys.stderr.write(str(len(r)) + '\n')
        align, insert = r[0:2]
        acc_chr, acc_site, donar_chr, donar_start, donar_end = insert[0], insert[1], align.refChr, align.refStart, align.refEnd
        #data = align.refChr, align.refStart, align.refEnd
        data = donar_chr, donar_start, donar_end, align.myName, align.myStart, align.myEnd
        print(*data, sep='\t', end='\n', file=outfile)
    outfile.close()

def logPositive(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align, te, insert, transduction = r[0:4]
        acc_chr, acc_site, donar_chr, donar_start, donar_end = insert[0], insert[1], align.refChr, align.refStart, align.refEnd
        data1 = acc_chr, acc_site, donar_chr, donar_start, donar_end, align.myName, align.myStart, align.myEnd
        data2 = te.repName, te.repClass, te.repFamily, te.genoStart, te.genoEnd
        data3 = insert[4:],transduction[0:]
        print(*data1, sep='\t', end='\n', file=outfile)
        print(*data2, sep='\t', end='\n', file=outfile)
        print(*data3, sep='\t', end='\n', file=outfile)
        outfile.write('\n')
    outfile.close()


def findTE(args):
    alignments = readAlignments(args[0])
    repeats = readRepeats(args[1])
    insertions = readInsertions(args[2])
    sys.stderr.write("insertions: " + str(len(insertions)) + '\n')
    outputfile = args[3]
    selectedPairs = selectedFromAlignments(alignments, insertions)
    sys.stderr.write("selectedPairs: " + str(len(selectedPairs)) + '\n')
    positive_full, positive_partial, negative = checkCandidates(selectedPairs, repeats)
    sys.stderr.write("positive_full: " + str(len(positive_full)) + '\n')
    sys.stderr.write("positive_partial: " + str(len(positive_partial)) + '\n')
    sys.stderr.write("negative: " + str(len(negative)) + '\n')
    logNegative(outputfile+"_negative", negative)
    logPositive(outputfile+"_positive_full", positive_full)
    logPositive(outputfile+"_positive_partial", positive_partial)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "alignments.maf repeats.txt insertionList outputfile"
    description = "Try to find transposable elements insertion."
    op = optparse.OptionParser (description=description)
    #op.add_option("-d", "--detail", help="Show detailed TE sequence")
    args = op.parse_args()[1]
    if len(args) < 4:
        op.error("please give me MAF alignments, repeat annnotation, insertionlist, and outputfilename")
    findTE(args)




