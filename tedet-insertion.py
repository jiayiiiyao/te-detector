from utilityFunctions import*

def getReferenceChr(alignment):
    return alignment.refChr

def tinyOverlap(a1, a2):
    if abs(a2.refStart - a1.refEnd) <= 20 and (a2.qryStart - a1.qryEnd) >= 100:
        return True
    else:
        return False

def largeOverlap(a1, a2):
    unaligned = a2.qryStart - a1.qryEnd
    overlap = a1.refEnd - a2.refStart
    if(a1.refStart < a2.refStart and a2.refEnd > a1.refEnd and a2.refStart < a1.refEnd) and (unaligned >= 100):
        return True
    else:
        return False

def largeGap(a1, a2):
    unaligned = a2.qryStart - a1.qryEnd
    gap = a2.refStart - a1.refEnd
    if gap > 20 and unaligned >= (gap+100):
        return True
    else:
        return False

def findInsertions(f, alignments):
    logging.info("Checking insertions...\n")
    #readlist = []
    insertions = collections.defaultdict(list)
    #gap = 150
    #gap = 50
    for readname in alignments.keys():
        alist = alignments[readname]
        if len(alist) >= 2:
            alist.sort(key=operator.attrgetter('qryStart'))
            for i in range(len(alist)-1):
                a1, a2 = alist[i], alist[i+1]
                if(a1.refChr == a2.refChr):
                    ref = a1.refChr
                    if tinyOverlap(a1, a2):
                        inslength = a2.qryStart - a1.qryEnd
                        tsdd = TSD(a2.refStart - a1.refEnd, None, None, None)
                        insertion = Insertion(ref, a2.refStart, readname, a1.qryEnd, inslength, tsdd)
                        insertions[ref].append(insertion)
                    elif largeOverlap(a1, a2):
                        unaligned = a2.qryStart - a1.qryEnd
                        overlap = a1.refEnd - a2.refStart
                        inslength = unaligned + overlap
                        tsdd =  TSD(-1*overlap, None, None, None)
                        insertion = Insertion(ref, a2.refStart, readname, a1.qryEnd, inslength, tsdd)
                        insertions[ref].append(insertion)
                    elif largeGap(a1, a2):
                        unaligned = a2.qryStart - a1.qryEnd
                        gap = a2.refStart - a1.refEnd
                        inslength = unaligned - gap
                        tsdd =  TSD(gap, None, None, None)
                        insertion = Insertion(ref, a2.refStart, readname, a1.qryEnd, inslength, tsdd)
                        insertions[ref].append(insertion)
    for ref in insertions.keys():
        insertions[ref].sort(key=operator.attrgetter('targetSite'))
    logInsertions(insertions, f)


def tedet_insertions(args): 
    alignments = readAlignments(openAndLog(args[0]))
    findInsertions(writeAndLog(args[1]), alignments)

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] maf-file"
    description = "Log insertions. "
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 2:
        op.error("Please give me alignments in MAF format and OutputFile")
    tedet_insertions(args)

