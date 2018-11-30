from tedet_insertion import*


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
        insertionstart, insertionend = insertion.insertionStart, insertion.insertionEnd
        donorstart = a.refStart
        donorend = a.refEnd
        refchr = a.refChr
        if refchr in repeats.keys(): 
            cans = getCandidateTE(repeats[refchr], donorstart)
            #cans = repeats[refchr]
            for te in cans:
                testart, teend = te.genoStart, te.genoEnd
                if donorstart <= te.genoStart and donorend >= te.genoEnd:
                    transduction = checkTransduction(a, insertionend, te) 
                    res = a, te, insertion, transduction                   
                    positive_full.append(res)
                    
                elif (te.genoStart <= donorstart and te.genoEnd <= donorend and te.genoEnd > donorstart):
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

def parseInsertions(insertionlist):
    insertions = {}
    for insert in insertionlist:
        readname = insert.readName
        ilist = insertions.get(readname)
        if not ilist:
            ilist = []
        #data = align1.refChr, align1.refEnd, align1.qryName, align1.qryEnd, align2.qryStart, tsd[0], tsd[1], tsd[2], tsd[3]
        ilist.append(insert)
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
                    insertstart, insertend = i.insertionStart, i.insertionEnd
                    if (abs(insertstart - readstart) <= 50 and abs(insertend - readend) <= 50) :
                        pair = align, i
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
        target_chr, target_site, donor_chr, donor_start, donor_end = insert.targetChr, insert.targetSite, align.refChr, align.refStart, align.refEnd
        data = donor_chr, donor_start, donor_end, align.qryName, align.qryStart, align.qryEnd
        print(*data, sep='\t', end='\n', file=outfile)
    outfile.close()

def logPositive(filename, result):
    outfile = writeAndLog(filename)
    for r in result:
        align, te, insert, transduction = r[0:4]
        tsd = insert.TSD
        target_chr, target_site, donor_chr, donor_start, donor_end = insert.targetChr, insert.targetSite, align.refChr, align.refStart, align.refEnd
        data1 = target_chr, target_site, donor_chr, donor_start, donor_end, align.qryName, align.qryStart, align.qryEnd, insert.insertionStart, insert.insertionEnd
        data2 = te.repName, te.repClass, te.repFamily, te.genoStart, te.genoEnd
        data3 = tsd.length, tsd.seq, tsd.left_flanking, tsd.right_flanking, transduction[0], transduction[1]
        print(*data1, sep='\t', end='\n', file=outfile)
        print(*data2, sep='\t', end='\n', file=outfile)
        print(*data3, sep='\t', end='\n', file=outfile)
        outfile.write('\n')
    outfile.close()



def findTE(args):
    alignments = readAlignments(args[0])
    repeats = readRepeats(args[1])
    insertions = parseInsertions(findInsertion(alignments))
    sys.stderr.write("insertions: " + str(len(insertions)) + '\n')
    outputfile = args[2]
    selectedPairs = selectedFromAlignments(alignments, insertions)
   # logInsertion(outputfile+"_insertion", selectedPairs)
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
    if len(args) < 3:
        op.error("please give me alignments in MAF format, RepeatAnnnotation, and OutputDirectory")
    findTE(args)




