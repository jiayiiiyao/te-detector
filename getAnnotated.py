from __future__ import print_function
import signal
import sys
import bisect
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
    alignments = {}
    blocks = readBlockFromMaf(f)
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
def readRepeats(f):
    repeats = {}
    for line in f:
        attributes = line.split()
        genoName = attributes[5]
        genoStart, genoEnd = int(attributes[6]), int(attributes[7])
        strand, repName, repClass, repFamily = attributes[9:13]
        if repClass == 'SINE' or repClass == 'LINE':
            transposon = Transposon(genoName, genoStart, genoEnd, strand, repName, repClass, repFamily)
            tlist = repeats.get(genoName)
            if tlist == None:
                tlist = []
            tlist.append(transposon)
            repeats[genoName] = tlist
    for chr in repeats.keys():
       repeats[chr].sort(key = operator.attrgetter('genoStart'))
    return repeats

def getRefStartPoints(alignlist):
	startpoints = []
	for align in alignlist:
		startpoints.append(align.refStart)
	return startpoints

def openAndLog(fileName):
    f = open(fileName)
    sys.stderr.write("Reading " + fileName + "..." + "\n")
    return f

def writeAndLog(fileName):
    f = open(fileName, 'w')
    sys.stderr.write("writing " + fileName + "..." + "\n")
    return f

def intersection(lst1, lst2):
	lst3 = [value for value in lst1 if value in lst2]
	return lst3

def getAnnotatedTEs(args):
	alignments = readAlignments(openAndLog(args[0]))
	repeats = readRepeats(openAndLog(args[1]))
	results = []
	#print(repeats.keys())
	#reflist = list(set(alignments.keys()) & set(repeats.keys()))

	for refchr in alignments.keys():
		if refchr in repeats.keys():
			alist = alignments[refchr]
			alist.sort(key = operator.attrgetter('refStart'))
			rlist = repeats[refchr]
			startpoints = getRefStartPoints(alist)
			asize = len(alist)
			for r in rlist:
				pos = bisect.bisect_right(startpoints, r.genoStart)
				if pos < asize:
					candidate = alist[pos-1]
					if candidate.refStart <= r.genoStart and candidate.refEnd >= r.genoEnd:
						res = candidate.refChr, candidate.refStart, candidate.refEnd, candidate.qryName, r.repName, r.repFamily, r.genoName, r.genoStart, r.genoEnd
						results.append(res)
						break
	outfile = writeAndLog(args[2])
	for r in results:
		print(*r, sep='\t', end='\n', file=outfile)



if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] maf-file"
    description = "Get annotated transposons from the reads"
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("please give me alignments in MAF format, TE annotation and Output name")
    getAnnotatedTEs(args)



