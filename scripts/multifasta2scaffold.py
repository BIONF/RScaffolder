'''
Created on Nov 3, 2015

@author: bardya
'''
import os
import argparse
from Bio import SeqIO
import glob
import itertools

def parse_args():
    parser = argparse.ArgumentParser(description='Given an ordered multifasta file, generate it to a scaffold, inserting Ns between them')
    
    parser.add_argument('files', metavar='/path/to/*.file', type=str, nargs='+',
                   help='path(s) to multifasta file(s)')
    
    parser.add_argument('--nlength', dest='nlength', metavar='<INT>', type=int, default=100,
                   help='length of repeated Ns to insert between sequences')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()


def collectSeqStrs(inputfile):
    """Takes an input multi-fasta file and concatenates its sequences"""
    seqs = []
    
    for record in SeqIO.parse(inputfile, "fasta"):
        seqs.append(record.seq)
        
    return [seq.__str__() for seq in seqs]


def mergeSequences(header, seqs, nlength):
    """Merge the sequences by inserting nlength (default:100) Ns inbetween"""
    Nstr = "N" * nlength
    merged_seq = Nstr.join(seqs)
    fasta_str = header + merged_seq
    formatted_fasta_str = fastaLinebreak(fasta_str)
    
    return formatted_fasta_str
 
    
def fastaLinebreak(s, seqrowlength=70):
    """Format the fasta file with linebreaks every 70th character)"""
    fastastring_formatted = []
    for line in s.split("\n"):
        if not line.startswith(">"):
            n = 70
            linesplit = [line[i:i+n] for i in range(0, len(line), n)]
            formattedline = "\n".join(linesplit) + "\n"
            fastastring_formatted.append(formattedline)
        else:
            fastastring_formatted.append(line + "\n")

    return "".join(fastastring_formatted)



if __name__ == '__main__':

    args = parse_args()
    
    file_list= []
   
    for g in args.files:
        file_list.append(glob.glob(g))
    
    input_files = list(itertools.chain.from_iterable(file_list))
    
    print('Output(s):\n')
    for input_file in input_files:
        try:
            inputfile = open(input_file, 'r')
            seqs = collectSeqStrs(inputfile)
            header = '>' + os.path.basename(inputfile.name).rsplit('.',1)[0] + '_length_' + str(sum([len(seq) for seq in seqs])) + '\n'
            if '.' in os.path.basename(inputfile.name): 
                outfile = open(inputfile.name.rsplit('.',1)[0] + '_reformatted.' + inputfile.name.rsplit('.',1)[1], 'w')
            else:
                outfile = open(inputfile.name + '_reformatted', 'w')
            outfile.write(mergeSequences(header, seqs, args.nlength))
            print(' - ' + outfile.name)
        except:
            print('Error occured wth the input file: ' + input_file)
