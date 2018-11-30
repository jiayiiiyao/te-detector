from utilityFuntions import*

def getAdjacent(alist, end):
    tsd = 30 # target site duplication here 
    startpoints = getRefStartPoints(alist)
    pos = bisect.bisect_right(startpoints, end)
    if pos != len(startpoints) and pos >= 2:
        for next in alist[pos - 1], alist[pos]:
            if abs(next.refStart - end) <= tsd:
                return next
    return None

def checkGap(align1, align2):
    bias = 200
    qryname1, qryname2 = align1.qryName, align2.qryName
    if qryname1 == qryname2:
        end = align1.qryEnd
        start = align2.qryStart
        if (start - end) > bias:
            return True
    return False

# return list of Insertion entity
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

def tedet_insertions(args):
    alignments = readAlignments(args[0])
    insertions = findInsertion(alignments)
    outputfile = args[1]
    logInsertion(outputfile, insertions)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "Alignments.maf OutputDirectory"
    description = "Try to find insertions in reads."
    op = optparse.OptionParser (description=description)
    args = op.parse_args()[1]
    if len(args) < 2:
        op.error("please give me alignments in MAF format and Output Directory")
    tedet_insertions(args)




