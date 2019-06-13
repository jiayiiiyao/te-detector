from __future__ import print_function
from utilityFunctions import *
from operator import itemgetter

MAX = 10E11

def getTargetSite(ins):
    return ins[1]

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


def readTarget(lines):
	targets = {}
	for line in lines:
		attributes = line.split()
		genoName, genoStart, genoEnd = attributes[0], int(attributes[3]), int(attributes[4])
		tlist = targets.get(genoName)
		if not tlist:
			tlist = []
		target = genoStart, genoEnd
		tlist.append(target)
		targets[genoName] = tlist
	for genoName in targets.keys():
		targets[genoName] = sorted(targets[genoName], key = itemgetter(0))
	return targets

def getNewSite(newsitemap, orisite):
	newsite = orisite
	for nmap in newsitemap:
		rangestart, rangeend, diff = nmap
		if orisite < rangestart:
			break
		elif orisite > rangestart and orisite < rangeend:
			newsite = orisite - diff
			break
	return newsite


def alterCoordinatesbyChr(targets, inslist):
	newsitemap = []
	diff = 0
	for i in range(len(targets)-1):
		genoStart, genoEnd = targets[i]
		nextStart, nextEnd = targets[i+1]
		diff = diff + genoEnd - genoStart
		tmap = genoStart, nextStart, diff
		newsitemap.append(tmap)
	fin = targets[len(targets)-1][0], MAX, diff+targets[len(targets)-1][1]-targets[len(targets)-1][0]
	newsitemap.append(fin)
	#print(newsitemap)
	newinslist = []
	for ins in inslist:
	 	targetChr, targetSite, readname, readStart, length, og = ins
	 	newsite = getNewSite(newsitemap, targetSite)
	 	newins = targetChr+'_m', newsite, readname, readStart, length, og
	 	newinslist.append(newins)
	return newinslist





# def alterCoordinatesbyChr(targets, inslist):
# 	newinslist = []
# 	diff = 0
# 	i, j = 0, 0
# 	while i < len(targets)-1:
# 		target = targets[i]
# 		diff = diff + target[1] - target[0]
# 		rangestart = target[0]
# 		rangeend = targets[i+1][0]
# 		print(rangestart)
# 		print(rangeend)
# 		while j < len(inslist):
# 			targetChr, targetSite, readname, readStart, length, og = inslist[j]
# 			if targetSite > rangestart and targetSite < rangeend:
# 				newsite = targetSite - diff
# 				newins = targetChr, newsite, readname, readStart, length, og
# 				newinslist.append(newins)
# 				j = j+1
# 			else:
# 				newinslist.append(inslist[j])
# 				i = i+1
# 				break
# 	return newinslist


def alterCoordinates(args):
	insertions = readInsertionFile(openAndLog(args[0]))
	targets = readTarget(openAndLog(args[1]))
	outfile = writeAndLog(args[2])
	for targetchr in targets.keys():
		altered = alterCoordinatesbyChr(targets[targetchr], insertions[targetchr])
		for a in altered:
		 	print(*a, sep='\t', end='\n', file=outfile)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] insertions targets output"
    description = "AlterCoordinates"
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("Please give me insertions, targets and output")
    alterCoordinates(args)