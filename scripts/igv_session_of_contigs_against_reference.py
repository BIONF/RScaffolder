'''
Created on May 31, 2016

@author: Ba1
'''
#!/usr/bin/python3

import os, subprocess, csv

cwd = os.path.dirname(__file__)

GCPRCNT_SCRIPT_PATH = os.path.join(cwd, 'gc_track_for_igv.sh')
GCPRCNT_WINDOW_SIZE = 1000
GCSKEW_SCRIPT_PATH = os.path.join(cwd, 'gcskew_track_for_igv.sh')
GCSKEW_WINDOW_SIZE = 5000

COORDS2GFF_SCRIPT_PATH = os.path.join(cwd, 'coords2gff3.py')
LENGTH_DISTR_SCRIPT_PATH = '/home/bardya/usr/scripts/bash_scripts/collect_seqtype_length_distribution_from_gff.sh'
QUAST_PATH = '/home/bardya/usr/bin/quast-4.4/quast.py'
IGVTOOLS_CMD = ['/share/applications/java/java8/bin/java', '-Xmx1500m', '-jar ~/usr/bin/IGVTools/igvtools.jar']

barebone_session_XML = open(os.path.join(cwd, 'assets' , 'xml_igv_session_template.txt'), 'r').read()

def getIgvFileMinMax(igvfilepath, colnum=4):
    """extracts the minimum and maximum of a IGV track file"""
    with open(igvfilepath) as infile:
        try:
            col3 = [float(row[4]) for row in list(csv.reader(infile, delimiter='\t'))]
        except:
            print("the generated igv file " + igvfilepath + " has at least one not float-convertable entry in col number " + str(colnum))
        
        return min(col3), max(col3)

def getIgvFileFirstValue(igvfilepath):
    """extracts the minimum and maximum of a IGV track file"""
    with open(igvfilepath) as infile:
        try:
            row0 = infile.readline()
            value0 = row0.strip().split('\t')[-1]
        except:
            print("the generated igv file " + igvfilepath + " has no float value in its last column of row 0 ")
        
        return str(value0)

#get the current directories subdirectories

def mappedContigCoverageTrack(mapped_features_gff_filepath, per_base_coverage_txt_filepath, outpath):
    mcontig_map = {}
    with open(mapped_features_gff_filepath) as mcontgff:
        for line in mcontgff:
            if not line.startswith('#'):
                lnlst = line.strip().split('\t')
                reference = lnlst[0]
                name = lnlst[-1].split(';',2)[1].rsplit('=')[1]
                start = lnlst[3]
                end = lnlst[4]
                mcontig_map[name]= (reference, start, end)
    
    #print(mcontig_map.keys())
        
    with open(per_base_coverage_txt_filepath, 'r') as covfile:
        with open(outpath, 'w') as outf:
            for line in covfile:
                seqname, position, coverage = line.strip().split('\t')
                if seqname in mcontig_map:
                    reference, start, end = mcontig_map[seqname]
                    corrected_start = int(start) + int(position) - 1    #in igv format a single position ranges from 0 to 1 (end is excluded)
                    corrected_end = int(start) + int(position)
                    output_line = '\t'.join([reference, str(corrected_start), str(corrected_end), 'Coverage', coverage]) + '\n'
                    outf.write(output_line)
                    


if __name__ == '__main__':
    
    #declare the working directory
    START_WD = '/home/bardya/usr/data/Reference_Assisted_Scaffolding/Gold_Standards'
    os.chdir(START_WD)
    subjectlist = next(os.walk('.'))[1] # os.walk gives you a tuple of root,dirs,files for the specified path
    if 'reference_genomes' in subjectlist:
        subjectlist.remove('reference_genomes')
    
    for subject in subjectlist:
        igv = {}        
        new_wd = os.path.join(os.getcwd(), subject, 'original_dna_and_gff')
        os.chdir(new_wd)
        igv["session_genome_id"] = subject
        igv["session_genome_path"] = os.path.join(new_wd, subject + '.fna')
        gff_filepath = os.path.join(new_wd, subject + '.gff')    

        #calculate and save path of the GC statistics tracks and get the min and max of each track        
        exit_code = subprocess.call(['/bin/bash', GCPRCNT_SCRIPT_PATH, igv["session_genome_path"], str(GCPRCNT_WINDOW_SIZE)])
        exit_code = subprocess.call(['/bin/bash', GCSKEW_SCRIPT_PATH, igv["session_genome_path"], str(GCSKEW_WINDOW_SIZE)])
        igv["resource_gcprcnt_path"] = os.path.join(new_wd, 'gcpercent', subject + '.fna_{}bps.igv'.format(GCPRCNT_WINDOW_SIZE))
        igv["resource_gcprcnt_min"], igv["resource_gcprcnt_max"] = getIgvFileMinMax(igv["resource_gcprcnt_path"])
        igv["resource_gcprcnt_id"] = igv["resource_gcprcnt_path"] + '_' + getIgvFileFirstValue(igv["resource_gcprcnt_path"])
        
        igv["resource_gcskew_path"] = os.path.join(new_wd, 'gcpercent', subject + '.fna_{}bps_gcskew.igv'.format(GCSKEW_WINDOW_SIZE))
        igv["resource_gcskew_min"], igv["resource_gcskew_max"] = getIgvFileMinMax(igv["resource_gcskew_path"])
        igv["resource_gcskew_id"] = igv["resource_gcskew_path"] + '_' + getIgvFileFirstValue(igv["resource_gcskew_path"])

        
        #filter RNA and repeats
        seqtype = 'repeat_region'
        outpath1 = os.path.join(new_wd, subject + '_{}.gff3'.format(seqtype))
        cmd1 = ['python3', COORDS2GFF_SCRIPT_PATH, '-gi', gff_filepath, '-t', seqtype,
            '>', outpath1]
        
        seqtype = 'rRNA tRNA'
        outpath2 = os.path.join(new_wd, subject + '_{}.gff3'.format('RNA'))
        cmd2 = ['python3', COORDS2GFF_SCRIPT_PATH, '-gi', gff_filepath, '-t', seqtype,
            '>', outpath2] 
        
        exit_code = subprocess.call(" ".join(cmd1), shell=True)
        exit_code = subprocess.call(" ".join(cmd2), shell=True)
        igv["resource_Repeats_features"] = outpath1
        igv["resource_RNA_features"] = outpath2
        
        #map contigs and create features file of contigs and one of the contig break regions
        os.chdir(os.path.join(START_WD, subject))
        new_wd = os.path.join(os.getcwd(), 'mapped_on_self')
        #coordsfile_path = os.path.join(new_wd, subject.replace('genomic', 'chromosome_filter_dnadiff' + '.mcoords'))
        #coordsfile_path = os.path.join(new_wd, subject.replace('genomic', 'chromosome_dnadiff' + '.mcoords'))
        coordsfile_path = os.path.join(new_wd, subject + '_dnadiff.mcoords')
        outpath1 = os.path.join(new_wd, 'mcoords_mapped_contigs')
        outpath2 = os.path.join(new_wd, 'mcoords_not_covered_regions')
        
        cmd1 = ['python3', COORDS2GFF_SCRIPT_PATH, '-i', coordsfile_path, '-f', '-c', '80', '-l', '250', '-u', '-o', outpath1]
        cmd2 = ['python3', COORDS2GFF_SCRIPT_PATH, '-i', coordsfile_path, '-f', '-br', '-c', '80', '-l', '250', '-u', '-o', outpath2]
        
        #create a file containing the length distribution
        outpath3 = os.path.join(new_wd, 'uncovered_seqs_lengths.txt')
        cmd3 = [LENGTH_DISTR_SCRIPT_PATH, "intercontig_region", outpath2 + '.gff3', outpath3]
        
        #create 2 gffs listing one collectiing all overlapping annotated features of breaks 
        #(i.e.)
        outpath4 = os.path.join(new_wd, 'breaks_in_annotated_regions')
        outpath5 = os.path.join(new_wd, 'breaks_in_unannotated_regions.gff3')
        cmd4 = ['python3', COORDS2GFF_SCRIPT_PATH, '-gi', outpath2 + '.gff3', '-gm', gff_filepath, '-f', '-go', outpath4, '>', outpath5]
        
        exit_code = subprocess.call(cmd1)
        exit_code = subprocess.call(cmd2)
        exit_code = subprocess.call(cmd3)
        exit_code = subprocess.call(" ".join(cmd4), shell=True)
        igv["resource_mapped_contigs"] = outpath1 + '.gff3'
        igv["resource_not_covered"] = outpath2 + '.gff3'
        
        
        #Parsing a per-base coverage file of the contigs and fit its positions to the genome positions of the mapped contigs
        new_wd = os.path.join(START_WD, subject)
        os.chdir(os.path.join(new_wd, START_WD))
        outpath6 = os.path.join(new_wd, 'mapped_coverage.igv')
        mappedContigCoverageTrack(outpath1 + '.gff3', os.path.join(new_wd, 'per_base_coverage.txt'), outpath6)
        igv["resource_coverage_min"], igv["resource_coverage_max"] = getIgvFileMinMax(outpath6)
        
        outpath7 = os.path.join(new_wd, 'mapped_coverage_sorted.tdf')
        exit_code = subprocess.call(" ".join(IGVTOOLS_CMD + ['sort', outpath6, os.path.join(new_wd, 'mapped_coverage_sorted.igv')]), shell=True)
        exit_code = subprocess.call(" ".join(IGVTOOLS_CMD + ['toTDF', os.path.join(new_wd, 'mapped_coverage_sorted.igv'), outpath7, 
                                     igv["session_genome_path"] + '.fai']), shell=True)
        
        igv["resource_coverage_path"] = outpath7
        igv["resource_coverage_id"] = os.path.basename(igv["resource_coverage_path"])
        
        exit_code = subprocess.call(["rm", outpath6, os.path.join(new_wd, 'mapped_coverage_sorted.igv') ])
        
        
        #FINALIZE
        igv["session_xml_path"] = os.path.join(START_WD, subject, 'igv_session.xml')
        os.chdir(START_WD)
        
        with open(outpath1 + '.gff3', 'r') as mapped:
            for line in mapped:
                if line.startswith("##sequence-region"):
                    linlst = line.strip().split(' ')
                    igv.update({"session_locus_name": linlst[1], "session_locus_start": linlst[2], "session_locus_end": linlst[3]})
        
        with open(igv["session_xml_path"], 'w') as outf:
            outf.write(barebone_session_XML.format(**igv))
        
        #### Doing QUAST analysis #####
#         new_wd = os.path.join(START_WD, subject)
#         os.chdir(new_wd)
#         env = os.environ.copy()
#         env["PYTHONPATH"]= "/usr/lib/python2.7"
#         ref_genome_path = os.path.join(new_wd, 'original_dna_and_gff', subject + '.fna')
#         ref_genes_path = os.path.join(new_wd, 'original_dna_and_gff', subject + '.gff')
#         scaffolds_filepath = os.path.join(new_wd, 'contigs.fasta')
#         cmd6 = ['python2', QUAST_PATH, scaffolds_filepath, '-R', ref_genome_path,
#                 '-G', ref_genes_path, '--gage']
#         exit_code = subprocess.call(cmd6, env=env)
        os.chdir(START_WD)

