from __future__ import print_function
from utilityFunctions import *
from operator import itemgetter

MAX = 10E11


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


def alterCoordinatesbyChr(targets, repeats):
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
	newreplist = []
	for r in repeats:
	 	genoName = r.genoName+'_m'
	 	length = r.genoEnd - r.genoStart
	 	newStart = getNewSite(newsitemap, r.genoStart)
	 	newEnd = newStart + length
	 	newrep = genoName, newStart, newEnd, r.strand, r.repName, r.repClass, r.repFamily
	 	newreplist.append(newrep)
	return newreplist


def alterCoordinates(args):
	repeats = readRepeats(openAndLog(args[0]))
	targets = readTarget(openAndLog(args[1]))
	outfile = writeAndLog(args[2])
	for chrs in repeats.keys():
		if chrs in targets.keys():
			altered = alterCoordinatesbyChr(targets[chrs], repeats[chrs])
			for a in altered:
				print(*a, sep='\t', end='\n', file=outfile)
		else:
			for r in repeats[chrs]:
				res = r.genoName, r.genoStart, r.genoEnd, r.strand, r.repName, r.repClass, r.repFamily
				print(*res, sep='\t', end='\n', file=outfile) 


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "%prog [options] annotation targets output"
    description = "AlterCoordinates"
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("Please give me annotation, targets and output")
    alterCoordinates(args)