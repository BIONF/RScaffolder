'''
Created on Jul 21, 2015

@author: bardya
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import os, argparse

def parse_args():
    parser = argparse.ArgumentParser(description="From a multi-fasta file extract first x and last x positions of the sequence and print")
    
    parser.add_argument('-i', dest='infilepath', metavar='</path/to/contigfile>', type=str,
                        help='path to an multi fasta file')
    
    parser.add_argument('-o', dest='outfilepath', metavar='</path/to/outfile>', type=str,
                        default='', 
                        help='path to output file prefix, meaning, the suffix ".fasta" will be added automatically'),
    
    parser.add_argument('-min', '--min-length', dest='min_size', metavar='<int>', type=int, 
                        default=1001, 
                        help='specify the minimum size of the sequence, shorter sequences will be ignored'),
    
    parser.add_argument('-l', '--ht-size', dest='ht_max_size', metavar='<int>', type=int,
                        default=1000,
                        help='specify the length of the head and tail subsequence to extract'),
                        
    parser.add_argument('-ms', '--margin-size', dest='margin_size', metavar='<int>', type=int,
                        default=0,
                        help='specify the margin size from the ends of the sequences'),
                        
    parser.add_argument('--version', action='version', version='0.02')
    
    return parser.parse_args()

seqlst = []
#Idea: Make variable sizes of ht
        
def getHtRecords(sequence_records_list, ht_max_size, margin_size):
    """
    For every 
    """
    htrecords = []
    
    for seqrec in sequence_records_list:
        
        if len(seqrec.seq) < ht_max_size * 2:
            ht_size = int((len(seqrec.seq) - 2 * margin_size)/2)
        else:
            ht_size = ht_max_size
                
        rec1 = SeqRecord(
                         seq = seqrec.seq[0 + margin_size : ht_size + margin_size],
                         id = seqrec.id + '_H', name=seqrec.name,
                         #description=seqrec.description)
                         description='')
        rec2 = SeqRecord(
                         seq = seqrec.seq[-ht_size - margin_size : len(seqrec.seq) - margin_size],
                         id = seqrec.id + '_T', name=seqrec.name,
                         #description=seqrec.description)
                         description='')
        
        htrecords.append(rec1)
        htrecords.append(rec2)
        
    return htrecords


def printRecs(recordlst):
    for record in recordlst:
        print(record.format("fasta").rstrip())

def toFile(recordlst,  outpath):
    with open(os.path.join(outpath + '_HT.fasta'), "w") as outf:
            SeqIO.write(recordlst, outf, "fasta")


def parseSeqs(filepath, min_size, ht_size, margin_size):
    
    # Read the file
    input_seq_it = SeqIO.parse(open(filepath, "r"), "fasta")

    # Build a list of short sequences:
    sequence_records_list = [record for record in input_seq_it \
                   if len(record.seq) > min_size]
    
    fname = os.path.basename(filepath.rsplit('.',1)[0])
    
    ##Check if all contig lengths allow to extract head and tail
    for seqrec in sequence_records_list:
        if len(seqrec.seq) < 2* margin_size + 2:
            raise ValueError('Aborted. The contig size of record ' + seqrec.name + ' is smaller than the selected head / tail size')
    
    return sequence_records_list, fname

if __name__ == '__main__':
    args = parse_args()
    #contigs_path = "/share/project/bardya/Acinetobacter/tmp/simulated_assemblies/GCF_000021145.1_ASM2114v1_genomic/assembly.spades.k22-127/contigs.fasta"
    sequence_records_list, fname = parseSeqs(args.infilepath, args.min_size, args.ht_max_size, args.margin_size)
    ht_records = getHtRecords(sequence_records_list, args.ht_max_size, args.margin_size)
    #printRecs(ht_records)
    #outpath='/home/bardya/usr/data/16_02/dev_tmp'
    if args.outfilepath:
        toFile(ht_records, args.outfilepath)
    else:
        printRecs(ht_records)
    
       