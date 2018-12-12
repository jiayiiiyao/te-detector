from __future__ import print_function
from Bio.SeqIO import *
from Bio.Seq import *
from Bio.SeqRecord import SeqRecord


import signal
import sys
import optparse
import operator
import signal

def readTarget(lines):
	targets = {}
	for line in lines:
		attributes = line.split()
		genoName, genoStart, genoEnd = attributes[5], int(attributes[6]), int(attributes[7])
		repName, repClass, repFamily = attributes[10:13]
		tlist = targets.get(genoName)
		if not tlist:
			tlist = []
		target = genoStart, genoEnd, repName, repClass, repFamily
		tlist.append(target)
		targets[genoName] = tlist
	return targets

#cut target transposon
def convertSeq(seq, start, end):
	length = end - start
	#mask = 'x'*length
	#seqNew = "".join((seq[:start], mask, seq[end:]))
	seqNew = "".join((seq[:start], seq[end:]))
	return seqNew

def convertReference(file, targets):
	records = []
	for chr_record in parse(file, 'fasta'):
		chrname = chr_record.id
		if chrname in targets.keys():
			tlist = targets[chrname]
			newSeq = str(chr_record.seq)
			for target in tlist:
				start, end = target[0], target[1]
				newSeq = convertSeq(newSeq, start, end)
			rec = SeqRecord(Seq(newSeq, None), id=chrname)
			records.append(rec)
			writeRecords_m(rec, "modified"+str(chrname))
		else:
			records.append(chr_record)
	return records

def writeRecords(records, outfile):
	write(records, outfile, "fasta")

def simulate(args):
	inputFile = args[0]
	targets = readTarget(open(args[1]))
	records= convertReference(inputFile, targets)
	writeRecords(records, args[2])


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "inFile.fa targets.txt outFile.fa"
    description = "Delete targeted transposon."
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 3:
        op.error("please give me input file in fasta, targetlist, and output file name")
    simulate(args)    

