'''
Created on Jun 22, 2016

@author: bardya
'''
import os, subprocess, argparse

barebone_script = """#!/bin/bash
#
#$-N JOB_subjectname_-_{toolname}
#$-S /bin/bash/
USER=bardya;
MY_DATE=$(/bin/date +%Y%m%d%H%M%S)
MY_HOST=$(hostname)
#$-cwd
#$-o {outdir}/'"$JOB_ID"'_subjectname_{timestamp}.out
#$-e {outdir}/'"$JOB_ID"'_subjectname_{timestamp}.err
#$-S /bin/bash
#-M djahanschiri@bio.uni-frankfurt.de
#$-m a
#$-V
#$-sync y
#$-q all.q
#-l h_vmem=2G,s_vmem=1.9G,h_rt=0:15:00
#-pe mpi XX

##added to allow for self deletion of script and outputs
currentscript="$0"; 
function finish { 
    echo "Deleting stdout and stderr files produced by ${currentscript}"; 
    rm {outdir}/${JOB_ID}_subjectname_{timestamp}.out; 
    rm {outdir}/${JOB_ID}_subjectname_{timestamp}.err;
    echo "Securely shredding ${currentscript}"; shred -u ${currentscript};
}


"""
# 
# WDIR=SCRIPTS_DIR/toolname/subjectname/
# mkdir -p $WDIR
# if [ ! -d $WDIR ]
# then
#   echo $WDIR not created
#   exit
# fi
# cd $WDIR

def parse_args():
    parser = argparse.ArgumentParser(description="A multi-fasta map will get mapped against all genomes in the given reference dir of the working directory ")
    
    parser.add_argument('-i', dest='infilepath', metavar='</path/to/multifastafile>', type=str,
                        help='path to an multi fasta file')
    
    parser.add_argument('-o', dest='outdir', metavar='</path/to/outdir>', type=str,
                        default='', 
                        help='path to output directory'),
                        
    parser.add_argument('-rd', dest='refdir', metavar='</path/to/refdir>', type=str,
                        default='', 
                        help='path to references directory with containing the fna sequence files to map against'),
                        
    parser.add_argument('--version', action='version', version='0.02')
    
    return parser.parse_args()


def generate_prescript(toolname, subjectname, outdir, mpi_slots=False):
    '''
    Based on 3 parameters and the barebone_script configs this generates a
    str including the sge settings and args for the cluster +
    the commands are added in this step.
    :param SCRIPTS_DIR:
    :param toolname:
    :param subjectname:
    '''
    import datetime
    now = datetime.datetime.now()
    strtime = now.strftime("%Y%m%d%H%M%S")
    
    res = barebone_script
    res = res.replace("subjectname",str(subjectname))
    res = res.replace("{toolname}",str(toolname))
    res = res.replace("{outdir}",str(outdir))
    res = res.replace("{timestamp}", strtime)
    if mpi_slots:
        res = res.replace('#-pe mpi XX','#$-pe mpi ' + str(mpi_slots))
    return res

options = "--gene-finding --debug" 
TOOLNAME = "Nucmer"
TOOLPATH = "/home/bardya/usr/bin/MUMmer3.23/"
TOOLPATH2 = "/bin/grep"
TOOLPATH3 = "/home/bardya/usr/bin/gnuplot"

def RunNucmer(contig, outdir, refdir):

    commands = []
    os.chdir(outdir)
    references = {}
    if refdir:
        references_dir = refdir
    else:
        references_dir = os.path.join(os.path.dirname(outdir), "references")  #find the references folder in the PROJECT_DIR (2 up)
    #Get all reference genomes
    for filename in os.listdir(references_dir):
        if filename.endswith(".fasta") or filename.endswith(".fna"):
            references[filename.rsplit('.',1)[0]] = os.path.join(references_dir, filename)
      
    for ref, val in references.items():
        prefix = os.path.join(outdir, ref)
                                                                    #ref, #query
        commands.append([TOOLPATH + 'nucmer', '-c', '55' ,'-l', '18', '-g', '120', '--maxmatch', '--prefix='+ prefix, val , contig])
        #commands.append([TOOLPATH + 'delta-filter', '-q', prefix + '.delta', '>', prefix + '.filter' ])
        #commands.append([TOOLPATH + 'show-coords', '-rcl', prefix + '.filter', '>', prefix + '_filter' + '.coords'])
        #commands.append([TOOLPATH + 'mapview', '-n', '3', '-f', 'pdf', '-Ir', '-I', prefix + '.coords', '-p', prefix + '_mapviewgraph' ])
        #commands.append([TOOLPATH + 'mapview', '-n', '3', '-f', 'pdf', '-Ir', '-I', '-p', ref + '_mapviewgraph',   prefix + '_filter' + '.coords'])
        #commands.append([TOOLPATH + 'mummerplot', '-postscript', '-p', prefix +'_plot', prefix + '.filter'])
        #commands.append([TOOLPATH3, prefix + '_plot' + '.gp'])
        #commands.append([TOOLPATH + 'dnadiff', '-d', prefix + '.filter', '-p', prefix + '_filter'+ '_dnadiff']) #executed on filtered

        commands.append([TOOLPATH + 'show-coords', '-rcl', prefix + '.delta', '>', prefix + '.coords'])
        commands.append([TOOLPATH + 'dnadiff', '-d', prefix + '.delta', '-p', prefix + '_dnadiff']) #executed on non filtered
        commands.append([TOOLPATH + 'show-tiling', '-a' , '-i', '90', '-l', '301', '-c', '-p', prefix + '_pseudo', prefix + '.delta', '>' , prefix + '_tiling_path.txt' ])
        
    if isinstance(commands, list):
        command = " && ".join([" ".join(cmd) for cmd in commands])
    else:
        command = commands
        
    script = generate_prescript(TOOLNAME, 'mapping', outdir)
    script += command
    scrpt_ending_name = "_" + TOOLNAME + '.sh'
    scriptfile = open(os.path.join(outdir, 'mapping' + scrpt_ending_name), "w")
    scriptfile.write(script)
    
    import stat
    st = os.stat(scriptfile.name)
    os.chmod(scriptfile.name, st.st_mode | stat.S_IEXEC)
    scriptfile.close()
    
    exit_code = subprocess.call(["qsub", scriptfile.name])
    
    return exit_code
    
if __name__ == '__main__':
    args = parse_args()
    #infilepath = "/home/bardya/usr/data/16_02/dev_tmp/contigs_HT.fasta"
    
    RunNucmer(args.infilepath, args.outdir, args.refdir)
    
    
# --maxmatch     Use all anchor matches regardless of their uniqueness

# -b int
# --breaklen     Distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)

# -c int
# --mincluster     Minimum cluster length (default 65)

# --[no]delta     Toggle the creation of the delta file. Setting --nodelta prevents the alignment extension step and only outputs the match clusters (default --delta)

# --depend     Print the dependency information and exit
# -d float

# --diagfactor     Maximum diagonal difference factor for clustering, i.e. diagonal difference / match separation (default 0.12)

# --[no]extend     Toggle the outward extension of alignments from their anchoring clusters. Setting --noextend will prevent alignment extensions but still align the DNA between clustered matches and create the .delta file (default --extend)

# -f
# --forward     Align only the forward strands of each sequence

# -g int
# --maxgap     Maximum gap between two adjacent matches in a cluster (default 90)

# -l int
# --minmatch     Minimum length of an maximal exact match (default 20)

#Decreasing the values of the -mincluster and --minmatch options will increase the sensitivity of the alignment but may produce less reliable alignments.
#In addition, significantly raising the value of the --maxgap value (say to 1000) can be crucial in producing alignments for more divergent genomes.