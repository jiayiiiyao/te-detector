from __future__ import print_function
from utilityFunctions import*
from operator import itemgetter


def checkTSD(align1, align2):
    tsd_length = align1.refEnd - align2.refStart
    if tsd_length > 0:
        tsd_seq = align2.refSeq[0:tsd_length]
        #tsd_left = align1.refSeq[-1 * (tsd_length + 10) : -1 * tsd_length]
        #tsd_right = align2.refSeq[tsd_length : tsd_length + 10]
    else:
        tsd_seq = 'N'
    return tsd_length, tsd_seq

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


def isIncluded(start, end, qrystart, qryend):
    if (start <= qrystart and qrystart < end) or (start < qryend and qryend <= end):
        return True
    else:
        return False

def ismapped(alist, start, end, strand):
    mapped = []
    for a in alist:
        if a.qryStrand == strand:
            qrystart = a.qryStart
            qryend = a.qryEnd
        else:
            qryend = a.qryFullLen - a.qryStart
            qrystart = a.qryFullLen - a.qryStart - a.qryLen            
        if isIncluded(start, end, qrystart, qryend):
            mapped.append(a)
    return mapped

def threeprimeintact(te, m):
    if te.genoStart >= m.refStart and te.genoStart <= m.refEnd:
        return True
    else:
        return False

def fiveprimeintact(te, m):
    if te.genoEnd >= m.refStart and te.genoEnd <= m.refEnd:
        return True
    else:
        return False

def twoendtruncated(te, m):
    if m.refStart >= te.genoStart and m.refEnd <= te.genoEnd:
        return True
    else:
        return False

def getAnnotatedTE(telist, m):
    annotated = []
    candidates = getCandidateTE(telist, m.refStart)
    for te in candidates:
        if (te.genoStart >= m.refStart and te.genoStart <= m.refEnd) or (te.genoEnd >= m.refStart and te.genoEnd <= m.refEnd) or (m.refStart >= te.genoStart and m.refEnd <= te.genoEnd):
        #if threeprimeintact(candidate, m) or fiveprimeintact(candidate, m) or twoendtruncated(candidate, m):
            annotated.append(te)
    return annotated

def getCandidateTE(telist, refstart):
    startpoints = []
    tes = []
    for te in telist:
        startpoints.append(te.genoStart)
    size = len(startpoints) 
    pos = bisect.bisect_right(startpoints, refstart)
    for i in range(pos - 1, pos + 1):
        if i > 0 and i < size:
            tes.append(telist[i])
    return tes

def checkTransduction(m, insend, te):
    if insend > m.qryEnd:
        transLength = m.refEnd - te.genoEnd
        seqstart = m.qryEnd - m.qryStart - transLength
        transSeq = m.qrySeq[seqstart:]
    else:
        transLength = 0
        transSeq = 'None'
    return transLength, transSeq

# check if mapped here 
def catInsertionsFromAlignments(alignments, telists):
    insertions = collections.defaultdict(list)
    temapped, tenotmapped, unmapped = [], [], []
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
                            insstart, insend = a1.qryEnd, a2.qryStart
                            strand = a1.qryStrand
                            insmapped = ismapped(alist, insstart, insend, strand)
                            if len(insmapped) == 0:
                                #alien insertion
                                unmappedIns = refChr, a1.refEnd, readname, insstart, unaligned
                                unmapped.append(unmappedIns)
                            else:
                                for m in insmapped:
                                    tsdLen, tsdSeq = checkTSD(a1, a2)
                                    mref = m.refChr
                                    tlist = telists[mref]
                                    tes = getAnnotatedTE(tlist, m)
                                    if len(tes) == 0:
                                        #te negative insertion
                                        res = refChr, a1.refEnd, a1.qryName, insstart, unaligned, tsdLen, m.qryStart, m.qryEnd-m.qryStart, mref, m.refStart, m.refEnd-m.refStart
                                        tenotmapped.append(res)
                                    else:
                                        for te in tes:
                                            transductionLength, transductionSeq = checkTransduction(m, insend, te)
                                            #te positive insertion
                                            res = refChr, a1.refEnd, a1.qryName, insstart, unaligned, tsdLen, tsdSeq, m.qryStart, m.qryEnd-m.qryStart, mref, m.refStart, m.refEnd-m.refStart, te.repName, te.genoStart, te.genoEnd, transductionLength, transductionSeq
                                            temapped.append(res)
    return temapped, tenotmapped, unmapped



def collectTransposon(args):
    print('start reading alignments...')
    readlist = readReadnames(openAndLog(args[1]))
    alignments = readtargetAlignments(openAndLog(args[0]), readlist)
    print('start reading repeats...')
    repeats = readRepeats(openAndLog(args[2]))
    outdir = args[3]
    print('start checking insertions...')
    temapped, tenotmapped, unmapped = catInsertionsFromAlignments(alignments,repeats)
    outfilea = writeAndLog(outdir+'-alien')
    for unmap in unmapped:
        print(*unmap, sep='\t', end='\n', file=outfilea)
    outfilea.close()
    outfiletn = writeAndLog(outdir+'-te-negative')
    for tenotmap in tenotmapped: 
        #targetchr, targetsite, readname, insertionstart, insertionlength, tsdlength
        #mapqrystart, mapqrylen, maprefchr, maprefstart, mapreflen
        res1 = tenotmap[0:6]
        res2 = tenotmap[6:] 
        print(*res1, sep='\t', end='\n', file=outfiletn)
        print(*res2, sep='\t', end='\n', file=outfiletn)
    outfiletn.write('\n')
    outfiletn.close()
    outfiletp = writeAndLog(outdir+'-te-positive')
    for temap in temapped: 
        #targetchr, targetsite, readname, insertionstart, insertionlength, tsdlength, tsdSeq
        #mapqrystart, mapqrylen, maprefchr, maprefstart, mapreflen
        #tename, testart, teend
        #transductionLength, transductionSeq 
        res1 = temap[0:7]
        res2 = temap[7:12] 
        res3 = temap[12:15] 
        res4 = temap[15:]
        print(*res1, sep='\t', end='\n', file=outfiletp)
        print(*res2, sep='\t', end='\n', file=outfiletp)
        print(*res3, sep='\t', end='\n', file=outfiletp)
        print(*res4, sep='\t', end='\n', file=outfiletp)
    outfiletp.write('\n')
    outfiletp.close()


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "alignments.maf repeatsAnnotation outputdir"
    description = "Try to find transposable elements insertion."
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("please give me alignments in MAF format, clusters, RepeatAnnotation, OutputDirectory")
    collectTransposon(args)    




