from utilityFunctions import*

def checkTSD(align1, align2):
    tsd_length = align1.refEnd - align2.refStart
    if tsd_length > 0:
        tsd_seq = align2.refSeq[0 : tsd_length]
        tsd_left = align1.refSeq[-1 * (tsd_length + 10) : -1 * tsd_length]
        tsd_right = align2.refSeq[tsd_length : tsd_length + 10]
    else:
        tsd_seq, tsd_left, tsd_right = "N", "N", "N"
    return tsd_length, tsd_seq, tsd_left, tsd_right

def getReferenceChr(alignment):
    return alignment.refChr

# check if mapped here 
def findInsertionFromRedo(alignments):
    logging.info("Checking insertions...\n")
    insertions = collections.defaultdict(list)
    mapped, unmapped = [], []
    tsd = 120
    gap = 150
    #gap = 50
    for readname in alignments.keys():
        alist = alignments[readname]
        if len(alist) >= 2:
            alist.sort(key=operator.attrgetter('refStart'))
            for ref, group in groupby(alist, getReferenceChr):
                group = list(group)
                gsize = len(group)
                if gsize >= 2:
                    for i in range(gsize-1):
                        a1, a2 = group[i], group[i+1]
                        if abs(a2.refStart - a1.refEnd) <= tsd and (a2.qryStart - a1.qryEnd) >= gap:
                            length, seq, left, right = checkTSD(a1, a2)
                            insertion = Insertion(ref, a1.refEnd, readname, a1.qryEnd, a2.qryStart-a1.qryEnd, TSD(length, seq, left, right))
                            pair = isMapped(alist, insertion)
                            if pair:
                                mapped.append(pair)
                            else:
                                unmapped.append(insertion)
    return mapped, unmapped

# check if annotated
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
        insertionstart, insertionend = insertion.insertionStart, insertion.insertionStart+insertion.insertionLength
        donorstart = a.refStart
        donorend = a.refEnd
        refchr = a.refChr
        if refchr in repeats.keys(): 
            cans = getCandidateTE(repeats[refchr], donorstart)
            for te in cans:
                testart, teend = te.genoStart, te.genoEnd
                if donorstart <= te.genoStart and donorend >= te.genoEnd:
                    transduction = checkTransduction_2(a, insertionend, te) 
                    res = a, te, insertion, transduction 
                    #remove duplicates?                  
                    positive_full.append(res)
                    
                elif (te.genoStart <= donorstart and donorend - te.genoEnd >= 20 and te.genoEnd > donorstart):
                    transduction = checkTransduction_2(a, insertionend, te)
                    res = a, te, insertion, transduction
                    #remove duplicates?    
                    positive_partial.append(res)
                      
                else:
                    res = a, insertion
                    #remove duplicates    
                    negative.append(res)
        else:
            res = a, insertion
            negative.append(res)
    return positive_full, positive_partial, negative

def isMapped(alist, ins):
    pair = None
    for a in alist:
        readstart, readend = a.qryStart, a.qryEnd
        insertstart, insertend = ins.insertionStart, ins.insertionStart + ins.insertionLength
        if (abs(insertstart - readstart) <= 50 and abs(insertend - readend) <= 50) :
            pair = a, ins
    return pair


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

def checkTransduction_2(align, insertionend, te):
    length = align.refEnd - te.genoEnd
    if length:
       start = te.genoEnd - align.refStart
       end = start + length
       #transductionSeq = align.refSeq[start:end]
       transductionSeq = align.qrySeq[-1*length:]
       transduction = length, transductionSeq
    else:
        transduction = 0, "None" 
    return transduction

def detectTransposon(args):
    alignments = readAlignments(openAndLog(args[0]))
    repeats = readRepeats(openAndLog(args[1]))
    outputfile = args[2]

    mapped, unmapped = findInsertionFromRedo(alignments)
    #sys.stderr.write("MappedInsertions: " + str(len(mapped)) + '\n')
    #logMapped(outputfile+"-mapped", mapped)
    #logUnmapped(outputfile+"-unmapped", unmapped)
    positive_full, positive_partial, negative = checkCandidates(mapped, repeats)
    #sys.stderr.write("positive-full: " + str(len(positive_full)) + '\n')
    #sys.stderr.write("positive-partial: " + str(len(positive_partial)) + '\n')
    logPositive(outputfile+"-positive-full", positive_full)
    logPositive(outputfile+"-positive-partial", positive_partial)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "alignments.maf repeatsAnnotation outputdir"
    description = "Try to find transposable elements insertion."
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("please give me alignments in MAF format, RepeatAnnotation, and OutputDirectory")
    detectTransposon(args)    




