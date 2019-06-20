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

def ismapped(alist, start, end):
    mapped = []
    for a in alist:
        qrystart = a.qryStart
        qryend = a.qryEnd           
        if isIncluded(start, end, qrystart, qryend):
            mapped.append(a)
    return mapped

def ismapped_strand(alist, start, end, strand):
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

def teStartpoints(telist):
    startpoints = []
    for te in telist:
        startpoints.append(te.genoStart)   
    return startpoints

def getAnnotatedTE(telist, m, tslist):
    annotated = []
    candidates = getCandidateTE(telist, m.refStart, tslist)
    for te in candidates:
        if (te.genoStart >= m.refStart and te.genoStart <= m.refEnd) or (te.genoEnd >= m.refStart and te.genoEnd <= m.refEnd) or (m.refStart >= te.genoStart and m.refEnd <= te.genoEnd):
        #if threeprimeintact(candidate, m) or fiveprimeintact(candidate, m) or twoendtruncated(candidate, m):
            annotated.append(te)
    return annotated

def getCandidateTE(telist, refstart, startpoints):
    tes = []
    size = len(startpoints) 
    pos = bisect.bisect_right(startpoints, refstart)
    for i in range(pos - 1, pos + 1):
        if i > 0 and i < size:
            tes.append(telist[i])
    return tes

def checkTransduction(m, insend, te):
    if insend > m.qryEnd:
        #transLength = m.refEnd - te.genoEnd
        transLength = insend - m.qryEnd
        #seqstart = m.qryEnd - m.qryStart - transLength
       # transSeq = m.qrySeq[seqstart:]
    else:
        transLength = 0
        #transSeq = 'None'
    #return transLength, transSeq
    return transLength

# check if mapped here 
def catInsertionsFromAlignments(alignments, repeats):
    insertions = collections.defaultdict(list)
    #get sorted TE startpoints array
    teStartPoints = collections.defaultdict(list)
    for genoname in repeats.keys():
        tlist = repeats[genoname]
        for te in tlist:
            teStartPoints[genoname].append(te.genoStart)
    #start
    temapped, tenotmapped, unmapped = [], [], []
    for readname in alignments.keys():
        alignmentsOnOneRead = alignments[readname]
        groupedAlignbychr = collections.defaultdict(list)
        for align in alignmentsOnOneRead:
            groupedAlignbychr[align.refChr].append(align)
        #groupedAlignbychr
        for refChr in groupedAlignbychr.keys():
            alignToeachChr = groupedAlignbychr[refChr]
            if len(alignToeachChr) >= 2:
                alignToeachChr.sort(key=operator.attrgetter('refStart'))
                for a1 in alignToeachChr:
                    for i in getAjacents(alignToeachChr, a1):
                        a2 = alignToeachChr[i]
                        og = a1.refEnd - a2.refStart
                        if not a2.qryStrand == a1.qryStrand:
                            a2qryStart = a2.qryFullLen - a2.qryStart - a2.qryLen 
                        else:
                            a2qryStart = a2.qryStart 
                        unaligned = a2qryStart - a1.qryEnd
                        if abs(og) <= 200 and unaligned >= 100:
                            insertStart, insertEnd = a1.qryEnd, a2.qryStart
                            strand = a1.qryStrand
                            insMapped = ismapped_strand(alignmentsOnOneRead, insertStart, insertEnd, strand)
                            if len(insMapped) == 0:
                                #alien insertion
                                unmappedIns = refChr, a1.refEnd, readname, insertStart, unaligned
                                unmapped.append(unmappedIns)
                            else:
                                for mapped in insMapped:
                                    tsdLen, tsdSeq = checkTSD(a1, a2)
                                    mappedRef = mapped.refChr 
                                    #print(mappedRef)
                                    tlist = repeats[mappedRef]
                                    tslist = teStartPoints[mappedRef]
                                    tes = getAnnotatedTE(tlist, mapped, tslist)
                                    if len(tes) == 0:
                                        #te negative insertion
                                        res = refChr, a1.refEnd, a1.qryName, insertStart, unaligned, tsdLen, mapped.qryStart, mapped.qryEnd-mapped.qryStart, mappedRef, mapped.refStart, mapped.refEnd-mapped.refStart
                                        tenotmapped.append(res)
                                    else:
                                        for te in tes:
                                            transductionLength = checkTransduction(mapped, insertEnd, te)
                                            #te positive insertion
                                            res = refChr, a1.refEnd, a1.qryName, insertStart, unaligned, tsdLen, tsdSeq, mapped.qryStart, mapped.qryEnd-mapped.qryStart, mappedRef, mapped.refStart, mapped.refEnd-mapped.refStart, te.repName, te.genoStart, te.genoEnd, te.genoEnd - te.genoStart, transductionLength
                                            temapped.append(res)
    return temapped, tenotmapped, unmapped



def collectTransposon(args):
    print('start reading alignments...')
    readlist = readReadnames(openAndLog(args[1]))
    print(len(readlist))
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
        #tename, testart, teend, transductionLength
        res1 = temap[0:7]
        res2 = temap[7:12] 
        res3 = temap[12:] 
        print(*res1, sep='\t', end='\n', file=outfiletp)
        print(*res2, sep='\t', end='\n', file=outfiletp)
        print(*res3, sep='\t', end='\n', file=outfiletp)
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




