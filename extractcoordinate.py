from __future__ import print_function
from operator import itemgetter
from Bio.SeqIO import *
from Bio.Seq import *
from Bio.SeqRecord import SeqRecord

import signal
import sys
import optparse
import operator

def extractcoordinate(args):
	inputfile, name, start, end, output = args
	records = []
	start = int(start)
	end = int(end)
	for seq_record in parse(inputfile, 'fasta'):
		seqname = seq_record.id
		if name in seqname:
			extractseq = str(seq_record.seq[start:end])
			seqid = seqname+": "+str(start)+" - "+str(end)
			rec = SeqRecord(Seq(extractseq,None), id=seqid)
			records.append(rec)
	write(records, output, 'fasta')

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    usage = "input name start end output"
    description = "Extract fasta file with name:start-end"
    op = optparse.OptionParser(description=description)
    args = op.parse_args()[1]
    if len(args) < 5:
        op.error("please give me input name start end and output")
    extractcoordinate(args) 

