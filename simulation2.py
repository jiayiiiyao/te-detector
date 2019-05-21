from __future__ import print_function
from operator import itemgetter
from Bio.SeqIO import *
from Bio.Seq import *
from Bio.SeqRecord import SeqRecord


import signal
import sys
import optparse
import operator


def readTarget(lines):
	targets = {}
	for line in lines:
		attributes = line.split()
		genoName, genoStart, genoEnd = attributes[0], int(attributes[2]), int(attributes[3])
		tlist = targets.get(genoName)
		if not tlist:
			tlist = []
		target = genoStart, genoEnd
		tlist.append(target)
		targets[genoName] = tlist
	for genoName in targets.keys():
		targets[genoName] = sorted(targets[genoName], key = itemgetter(0))
	return targets

#cut target transposon
def convertSeq(seq, start, end):
	length = end - start
	mask = 'x'*length
	seqNew = "".join((seq[:start], mask, seq[end:]))
	#seqNew = "".join((seq[:start], seq[end:]))
	return seqNew

def convertReference(file, targets):
	records = []
	chrLen = {}
	for chr_record in parse(file, 'fasta'):
		chrname = chr_record.id
		chrLen[chrname] = len(chr_record.seq)
		if chrname in targets.keys():
			tlist = targets[chrname]
			newSeq = str(chr_record.seq)
			for target in tlist:
				start, end = target[0], target[1]
				newSeq = convertSeq(newSeq, start, end)
			newSeq_m = newSeq.replace('x', '')
			print(str(len(newSeq)) + ',' + str(len(newSeq_m)))
			chrname_m = chrname + '_m'
			rec = SeqRecord(Seq(newSeq_m , None), id=chrname_m)
			records.append(rec)
			#writeRecords(rec, "modified_"+str(chrname))
		else:
			records.append(chr_record)
	return records, chrLen

def writeRecords(records, outfile):
	write(records, outfile, "fasta")

def calculate(telist, chrlen):
	length = []
	origin = []
	modified = []
	count = len(telist)
	for i in range(count):
		if i == 0:
			#print(telist[0])
			length.append(telist[0][0]) 
			origin.append(0)
			modified.append(0)
			#origin[0] = 0
			#modified[0] = 0
		else:
			#length[i] = telist[i][0] - telist[i-1][1]
			length.append(telist[i][0] - telist[i-1][1])
			#origin[i] = telist[i-1][1]
			origin.append(telist[i-1][1])
			#modified[i] = modified[i-1] + lenth[i-1] + 1
			modified.append(modified[i-1] + length[i-1] + 1)
	#length[count] = chrlen - telist[count-1][1]
	length.append(chrlen - telist[count-1][1])
	origin.append(telist[count-1][1])
	modified.append(modified[count-1] + length[count-1] + 1)
	return origin, modified, length

def writeSeg(filename, origin, modified, length, genoName, genoName_m):
	f = open(filename, 'a+')
	count = len(length)
	print(length)
	print(origin)
	print(modified)
	for i in range(count):
		#print(i)
		result = length[i], genoName, origin[i], genoName_m, modified[i]
		print(*result, sep='\t', end='\n', file=f)
	f.close()

def simulate(args):
	inputFile = args[0]
	targets = readTarget(open(args[1]))
	#print (targets)
	records, chrLen = convertReference(inputFile, targets)
	writeRecords(records, args[2])
	seg = args[3]
	for genoName in targets.keys():
		origin, modified, length = calculate(targets[genoName], chrLen[genoName])
		genoName_m = genoName+'_m'
		writeSeg(seg, origin, modified, length, genoName, genoName_m)


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "inFile.fa targets.txt outFile.fa align.seg"
    description = "Delete targeted transposon."
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 4:
        op.error("please give me input file in fasta, targetlist, output fasta and output seg")
    simulate(args)    

 