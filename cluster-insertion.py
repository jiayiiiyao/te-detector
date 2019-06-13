from __future__ import print_function
from utilityFunctions import *

def readInsertionFile(f):
    insertions = collections.defaultdict(list)
    for line in f:
        attributes = line.split()
        targetChr, readname = attributes[0], attributes[2]
        targetSite, readStart, length, og = int(attributes[1]), int(attributes[3]), int(attributes[4]), int(attributes[5])
        ins = targetChr, targetSite, readname, readStart, length, og
        insertions[targetChr].append(ins)
    for targetChr in insertions.keys():
    	insertions[targetChr].sort(key=getTargetSite)
    return insertions

def getTargetSite(ins):
    return ins[1]

def getInsLength(ins):
    return ins[4]

def clusterInsertion(insertions):
    interval = 50
    clusterbysite = collections.defaultdict(list)
    clusters = collections.defaultdict(list)
    for targetChr in insertions.keys():
        tsite = 0
        insertionlist = insertions[targetChr]
        for insert in insertionlist:
            insertlength = insert[4]
            if getTargetSite(insert) - tsite > interval:
                tsite = getTargetSite(insert)
            cls = targetChr, tsite
            clusterbysite[cls].append(insert)
    for cls in clusterbysite.keys():
        clusterbysite[cls].sort(key=getInsLength)
        plength = 0
        for ins in clusterbysite[cls]:
            if abs(getInsLength(ins) - plength) > 50:
                plength = getInsLength(ins)
            targetChr, tsite = cls
            ccls = targetChr, tsite, plength
            clusters[ccls].append(ins)
    return clusters

def logCluster(clusters, f, count):
    for cls in clusters.keys():
        targetChr, tsite, plength = cls
        csize = len(clusters[cls])
        if csize >= count:
            avglength = 0
            for insert in clusters[cls]:
                avglength += getInsLength(insert)
            avglength = int(avglength/csize)
            for insert in clusters[cls]:
                targetChr, targetSite, readname, readStart, length, og = insert 
                res = targetChr, tsite, avglength, csize, targetSite, readname, readStart, length, og
                print(*res, sep='\t', end='\n', file=f)

def logInsertionClusters(args):
	infile, count, outfile = args
	insertions = readInsertionFile(openAndLog(infile))
	logCluster(clusterInsertion(insertions), writeAndLog(outfile), int(count))

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] insertions count outputfile"
    description = "Log insertions clusters"
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 2:
        op.error("Please give me insertions, cluster size limit and outfile")
    logInsertionClusters(args)
