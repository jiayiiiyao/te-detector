from __future__ import print_function
from utilityFunctions import*

# For bisecting 
def getRefStartPoints(alignlist):
    startpoints = []
    for align in alignlist:
        startpoints.append(align.refStart)
    return startpoints

def getAdjacent(alist, end):
    candidates = []
    tsd = 30 # target site duplication here 
    startpoints = getRefStartPoints(alist)
    pos = bisect.bisect_right(startpoints, end)
    if pos != len(startpoints) and pos >= 2:
        for next in alist[pos - 1], alist[pos]:
            if abs(next.refStart - end) <= tsd:
                candidates.append(next)
    return candidates

def logInsertionReads(f, alignments):
    readlist = []
    for chr in alignments.keys():
        alist = alignments[chr]
        for align1 in alist:
            end = align1.refEnd
            candidates = getAdjacent(alist, end)
            if len(candidates) > 0:
                for align2 in candidates:
                    isgap = checkGap(align1, align2)
                    if isgap:
                        readlist.append(align1.qryName)
    logReads(f, list(set(readlist)))


def tedet_insertions(args):
    f = openAndLog(args[0])
    alignments = readAlignments(f)
    f.close()
    f2 = writeAndLog(args[1])
    logInsertionReads(f2, alignments)
    f2.close()

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] maf-file"
    description = "Find insertions in reads."
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 2:
        op.error("please give me alignments in MAF format and Output Directory")
    tedet_insertions(args)




