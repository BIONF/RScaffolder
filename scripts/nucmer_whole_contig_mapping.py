'''
Created on Feb 24, 2017

@author: bardya
'''

import subprocess
import os

def run_full_nucmer(subject, ref, outdir):
    
    base = os.path.basename(ref.rsplit('.', 1)[0])
    outbase = os.path.join(outdir,base)
    
    cmd = ['/home/bardya/usr/bin/MUMmer3.23/nucmer']
    cmd.append("-c 65 -l 20 --maxmatch")
    cmd.append("--prefix={}".format(outbase)) # How to call the output files
    cmd.append(ref)
    cmd.append(subject)
    cmd.append("&&")
    cmd.append('/home/bardya/usr/bin/MUMmer3.23/delta-filter')
    cmd.append('-q {} > {}'.format(outbase + '.delta', outbase + '.filter'))
    cmd.append("&&")
    cmd.append('/home/bardya/usr/bin/MUMmer3.23/show-coords')
    cmd.append('-rcl {} > {}'.format(outbase + '.filter', outbase + '_filter.coords'))
    #cmd.append("&&")
    #cmd.append('/home/bardya/usr/bin/MUMmer3.23/mapview')
    #cmd.append('-n 3 -f pdf -Ir -I -p {} {}'.format(outbase + '_mapviewgraph', outbase + '_filter.coords'))
    #cmd.append("&&")
    #cmd.append('/home/bardya/usr/bin/MUMmer3.23/mummerplot')
    #cmd.append('-postscript -p {} {}'.format(outbase + '_plot',  outbase + '.filter'))
    #cmd.append('/home/bardya/usr/bin/gnuplot {}'.format(outbase + '_plot.gp'))
    cmd.append("&&")
    cmd.append('/home/bardya/usr/bin/MUMmer3.23/dnadiff')
    cmd.append('-d {} -p {}'.format(outbase + '.filter', outbase +  '_filter_dnadiff'))
    cmd.append("&&")
    cmd.append('/home/bardya/usr/bin/MUMmer3.23/show-coords')
    cmd.append('-rcl {} > {} '.format(outbase + '.delta', outbase + '.coords'))
    cmd.append("&&")
    cmd.append('/home/bardya/usr/bin/MUMmer3.23/dnadiff')
    cmd.append('-d {} -p {}'.format(outbase + '.delta', outbase +  '_dnadiff'))
    #cmd.append("&&")
    #cmd.append('/home/bardya/usr/bin/MUMmer3.23/show-tiling')
    #cmd.append('-a -i 90 -l 301 -c -p {} > {}'.format(outbase + 'pseudo', outbase + '.delta', outbase + '_tiling_path.txt'))

    return cmd


if __name__ == '__main__':
    
    subject = '/home/bardya/usr/data/Reference_Assisted_Scaffolding/Gold_Standards/Acinetobacter_baumannii_ACICU/contigs.fasta'
    ref = '/home/bardya/usr/data/Reference_Assisted_Scaffolding/Gold_Standards/Acinetobacter_baumannii_ACICU/original_dna_and_gff/Acinetobacter_baumannii_ACICU.fna'
    outdir = '/home/bardya/usr/data/Reference_Assisted_Scaffolding/Gold_Standards/Acinetobacter_baumannii_ACICU/mapped_on_self'
    
    cmd = " ".join(run_full_nucmer(subject, ref, outdir))
    #cmd = "\n".join().split('&&')
    print(cmd)
    subprocess.call(cmd)
    