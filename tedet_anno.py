from __future__ import print_function

import optparse
import operator
import bisect
import signal
import sys

class Alignment:
    aligncount = 0
    def __init__(self, refChr, refStart, refEnd, refLen, refSeq,
                 qryName, qryStart, qryEnd, qryLen, qrySeq):
        self.refChr, self.refStart, self.refEnd, self.refLen, self.refSeq = refChr, refStart, refEnd, refLen, refSeq
        self.qryName, self.qryStart, self.qryEnd, self.qryLen, self.qrySeq= qryName, qryStart, qryEnd, qryLen, qrySeq
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
        # ignore Simple_repeat and Low_complexity
        #if repFamily != 'Simple_repeat' and repFamily != 'Low_complexity':
        if not repFamily == ' ':
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
    for i in range(pos - 2, pos + 2):
        if i >= 0 and i < size:
            tes.append(telist[i])
    return tes

def checkCandidates(pairs, repeats):
    bias = 100
    positive_full = []
    positive_partial = []
    negative = []
    for pair in pairs:
        a, insertion = pair[0], pair[1] 
        insertionstart, insertionend = int(insertion[3]), int(insertion[4])
        donarstart = a.refStart
        donarend = a.refEnd
        refchr = a.refChr
        if refchr in repeats.keys(): 
            cans = getCandidateTE(repeats[refchr], donarstart)
            #cans = repeats[refchr]
            for te in cans:
                testart, teend = te.genoStart, te.genoEnd
                if donarstart <= te.genoStart and donarend >= te.genoEnd:
                    transduction = checkTransduction(a, insertionend, te) 
                    res = a, te, insertion, transduction                   
                    positive_full.append(res)
                    
                elif (te.genoStart <= donarstart and te.genoEnd <= donarend and te.genoEnd > donarstart):
                    transduction = checkTransduction(a, insertionend, te)
                    res = a, te, insertion, transduction
                    positive_partial.append(res)
                      
                else:
                    res = a, insertion
                    negative.append(res)
                    
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
        if len(line):
            fields = line.split()
            readname = fields[3]
            ilist = insertions.get(readname)
            if ilist == None:
                ilist = []
            # if len(fields) == 10:
            #     data = fields[0], fields[1], fields[3], fields[4], fields[5], fields[6], fields[7], fields[8], fields[9]
            # else:
            data = fields[0], fields[1], fields[3], fields[4], fields[5], fields[6], fields[7], fields[8], fields[9]
            ilist.append(data)          
            insertions[readname] = ilist
    return insertions

def selectedFromAlignments(alignments, insertions):
    selected = []
    for alignlist in alignments.values():
        for align in alignlist:
            readname, readstart, readend = align.qryName, align.qryStart, align.qryEnd
            if readname in insertions.keys():
                ilist = insertions[readname]
                for i in ilist:
                    insertstart, insertend = int(i[3]), int(i[4])
                    if (abs(insertstart - readstart) <= 50 and abs(insertend - readend) <= 50) :
                        pair = align, i, readname
                        selected.append(pair)
    return selected


def checkTransduction(align, insertionend, te):
    length = align.refEnd - te.genoEnd -(align.qryEnd - insertionend)
    if length:
       start = te.genoEnd - align.refStart
       end = start + length
       transductionSeq = align.refSeq[start:end]
       transduction = length, transductionSeq
    else:
        transduction = 0, "None" 
    return transduction


def logNegative(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align, insert = r[0:2]
        acc_chr, acc_site, donar_chr, donar_start, donar_end = insert[0], insert[1], align.refChr, align.refStart, align.refEnd
        data = donar_chr, donar_start, donar_end, align.qryName, align.qryStart, align.qryEnd
        print(*data, sep='\t', end='\n', file=outfile)
    outfile.close()

def logPositive(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align, te, insert, transduction = r[0:4]
        acc_chr, acc_site, donar_chr, donar_start, donar_end = insert[0], insert[1], align.refChr, align.refStart, align.refEnd
        data1 = acc_chr, acc_site, donar_chr, donar_start, donar_end, align.qryName, align.qryStart, align.qryEnd, insert[3], insert[4]
        data2 = te.repName, te.repClass, te.repFamily, te.genoStart, te.genoEnd
        data3 = insert[5], insert[6], insert[7], insert[8], transduction[0], transduction[1]
        print(*data1, sep='\t', end='\n', file=outfile)
        print(*data2, sep='\t', end='\n', file=outfile)
        print(*data3, sep='\t', end='\n', file=outfile)
        outfile.write('\n')
    outfile.close()

def logInsertion(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align = r[0]
        aligndata = align.qryName, align.qryStart, align.qryEnd, align.refChr, align.refStart, align.refEnd
        insert = r[1][0:]
        print(*aligndata, sep='\t', end='\n', file=outfile)
        print(*insert, sep='\t', end='\n', file=outfile)
        outfile.write('\n')
    outfile.close()


def findTE(args):
    alignments = readAlignments(args[0])
    repeats = readRepeats(args[1])
    insertions = readInsertions(args[2])
    sys.stderr.write("insertions: " + str(len(insertions)) + '\n')
    outputfile = args[3]
    selectedPairs = selectedFromAlignments(alignments, insertions)
    logInsertion(outputfile+"_insertion", selectedPairs)
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
    usage = "alignments.maf repeatsAnnotation insertionList outputdir"
    description = "Try to find transposable elements insertion."
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 4:
        op.error("please give me alignments in MAF format, RepeatAnnnotation, Insertionlist, and OutputDirectory")
    findTE(args)




