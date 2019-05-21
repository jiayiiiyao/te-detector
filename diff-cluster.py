from __future__ import print_function
from utilityFunctions import *


#key=(chr, targetsite, length)
def readClusters(f):
    clusters = collections.defaultdict(list)
    for line in f:
        attributes = line.split()
        chromosome, clsite, avglength = attributes[0],int(attributes[1]), int(attributes[2])
        cls = chromosome, clsite, avglength
        tsite, readname, startcoordinate, length, og = int(attributes[4]), attributes[5], int(attributes[6]), int(attributes[7]), int(attributes[8])
        ins = tsite, readname, startcoordinate, length, og
        clusters[cls].append(ins)
    return clusters


def isIncluded(p, plist):
    res = False
    tchr1, tsite1 = p
    for p2 in plist:
        tchr2, tsite2 = p2
        if tchr1 == tchr2 and abs(tsite1 - tsite2) <= 30:
            res = True
    return res

def diffCluster(cluster1, cluster2, outfile):
    for cls1 in cluster1:
        chromosome1, clsite1, avglength1 = cls1
        found = False
        for cls2 in cluster2:
            chromosome2, clsite2, avglength2 = cls2
            if (chromosome1 == chromosome2) and (abs(clsite2-cliste1)<= 30) and (abs(avglength1-avglength2)<= 30):
                found = True
                break
        if found == False:
            for insert in clusters[cls1]:
                tsite, readname, startcoordinate, length, og = insert 
                res = chromosome1, clsite1, avglength1, len(clusters[cls1]), tsite, readname, startcoordinate, length, og
                print(*res, sep='\t', end='\n', file=outfile)


def diff_insertions(args): 
    file1, file2, outfile = args[0], args[1], args[2]
    cluster1, cluster2 = readClusters(openAndLog(file1)), readClusters(openAndLog(file2))
    logDiffs(diffCluster(cluster1, cluster2), writeAndLog(outfile))

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] insertions1 insertions2"
    description = "Find unique insertions."
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("Please give me two cluster files and outfile")
    diff_insertions(args)




