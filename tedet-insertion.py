from __future__ import print_function
import itertools as it
from utilityFunctions import*


def getReferenceChr(alignment):
    return alignment.refChr

def getStartCoordinates(list):
    starts = [a.refStart for a in list]
    for a in list:
        starts.append(a.refStart)
    return starts

def getAjacents(list, a1):
    res = []
    starts = [a.refStart for a in list]
    a1end = a1.refEnd
    i = bisect.bisect(starts, a1end)
    if i <= len(list)-1:
        res.append(i-1)
        res.append(i)
    else:
        res.append(i-1)
    return res

def collectInsertions(f, alignments):
    logging.info("Checking insertions...\n")
    insertions = collections.defaultdict(list) 
    for readname in alignments.keys():
        alist = alignments[readname]
        groupedlistbychr = collections.defaultdict(list)
        for a in alist:
            groupedlistbychr[a.refChr].append(a)
        for refChr in groupedlistbychr.keys():
            alignsOneRef = groupedlistbychr[refChr]
            if len(alignsOneRef) >= 2:
                alignsOneRef.sort(key=operator.attrgetter('refStart'))
                for a1 in alignsOneRef:
                    for i in getAjacents(alignsOneRef, a1):
                        a2 = alignsOneRef[i]
                        og = a1.refEnd - a2.refStart
                        unaligned = a2.qryStart - a1.qryEnd
                        if abs(og) <= 500 and unaligned >= 50:
                            tsdd = TSD(og, None)
                            insertion = Insertion(refChr, a1.refEnd, readname, a1.qryEnd, unaligned, tsdd)
                            insertions[refChr].append(insertion)
    for ref in insertions.keys():
        insertions[ref].sort(key=operator.attrgetter('targetSite'))
    logInsertions(insertions, f)    


def tedet_insertions(args): 
    alignments = readAlignments(openAndLog(args[0]))
    collectInsertions(writeAndLog(args[1]), alignments)

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] maf-file"
    description = "Log insertions. "
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 2:
        op.error("Please give me alignments in MAF format and OutputFile")
    tedet_insertions(args)

