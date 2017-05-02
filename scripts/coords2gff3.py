'''
Created on Dec 15, 2015

This program parses a coords file which is the output of a nucmer run into a valid gff3 file
with the coordinate system of the reference which can then be used in IGV for display for ex-
ample.

@author: bardya
'''
import os
import argparse
import sys
from collections import OrderedDict

GFF3_FORMAT_ORDER = ["seqid", "source", "seqtype", "start", "end", "score", "strand", "phase", "seqattributes"]
GFF3_SEP = "\t"
GFF3_END = "\n"
GFF3_UNDEF = "."
GFF3_FIRST_LINE = "##gff-version 3\n"
#"#!gff-spec-version 1.21"
#"#!processor NCBI annotwriter"
#"#!genome-build ASM18720v4"
#"#!genome-build-accession {}"
GFF3_SECOND_LINE = "##sequence-region {} {} {}\n"
GFF3_LAST_LINE = "###"

def parse_args():
    parser = argparse.ArgumentParser(description="Convert coords files produced by MUMMER or nucmer into GFF file format")
    
    parser.add_argument('-i', dest='infilepath', metavar='</path/to/coordsfile>', type=argparse.FileType('rt'),
                        help='path to an coords file')
    
    parser.add_argument('-o', dest='outfilepath', metavar='</path/to/outfile>', type=str,
                        default='', 
                        help='path to output file prefix, meaning, the suffix ".gff3" will be added automatically'),
    
    parser.add_argument('-source', dest='source', metavar='<source>', type=str, 
                        default='contig', 
                        help='specify what should be written in the "source" field of the GFF output'),
    
    parser.add_argument('-c', '--coverage-threshold', dest='coverage_threshold', metavar='<float>', type=float,
                        default=0.0,
                        help='specify a coverage threshold. All features with less coverage will not be imported as features'),
    
    parser.add_argument('-l', '--length-threshold', dest='length_threshold', metavar='<int>', type=int,
                        default=0,
                        help='specify a length threshold. All features with less length will not be imported as features'),
    
    parser.add_argument('-ht', '--ht-flag', dest='htflag', action='store_true',
                        help='if flag is set, it assumed that the contig names end with "_H" or "_T" which will be concatenated to the derived NODE ids'),
                        
    parser.add_argument('-u', '--only-unique', dest='unique_requirement', action='store_true', 
                        help='if specified, requires a feature to be unique or its parts are overlapping, other thresholds are not applied on overlapping'),
                        
    parser.add_argument('-br', '--breaks-only', dest='breaks_only', action='store_true',
                        help='run with this flag to generate a gff file which has the breaks and the uncovered region as features'),
                         
    parser.add_argument('-gi', '--gff-in', dest='gffinfilepath', metavar='</path/to/gff-infile>', type=argparse.FileType('rt'),
                        default=None,
                        help='specify a path to the gff file to parse'),
                   
    parser.add_argument('-go', '--gff-out', dest='gffoutfilepath', metavar='</path/to/gff-outfile>', type=str, 
                        default='',
                        help='path to output file (just the prefix)'),
                    
    parser.add_argument('-t', '--gff-type-include', dest='gff_type_include', nargs='+', metavar='<type1 type2 ...>',
                        default=None,
                        help='give a list of space-separated gff seqtypes to include in the gff output'),
                    
    parser.add_argument('-T', '--gff-type-exclude', dest='gff_type_exclude', nargs='+', metavar='<type1 type2 ...>',
                        default=None,
                        help='give a list of space-separated gff seqtypes to exclude in the gff output (overwrites -t)'),
                        
    parser.add_argument('-a', '--gff-attribute-include', dest='gff_attr_include', nargs='+', metavar='<str1 str2 ...>',
                        default=None,
                        help='if any of these strings was found anywhere in the gff column "attributes" the feature will be included in the gff output'),
                    
    parser.add_argument('-A', '--gff-attribute-exclude', dest='gff_attr_exclude', nargs='+', metavar='<str1 str2 ...>',
                        default=None,
                        help='if any of these strings was found anywhere in the gff column "attributes" the feature will be excluded from the gff output (overwrites -a)'),
                        
    parser.add_argument('-gm', '--gff-map-in', dest='mappingGFFpath', metavar='</path/to/gff-infile>', type=argparse.FileType('rt'),
                        default=None,
                        help='specify a path to a gff file we would like to map our first gff file against'), 
    
    parser.add_argument('-x', '--only-genes-in-aligned-regions', dest='only_aligned', action='store_true',
                        help='if set, only genes in aligned regions of the contig mapping will be displayed'),
                        
    parser.add_argument('-f', '--force-overwrite', dest='force_ow', action='store_true',
                        help='if set output files will be forced to overwrite existing files.'),
                        
    parser.add_argument('--version', action='version', version='0.13')
    
    return parser.parse_args()


class Feature:
    
    def __init__(self, **kwargs):
        """For every possible feature:value pair in the given dictionary, create a 
        a instance attribute for easier access"""
        
        for k,val in kwargs.items():      #go through all parameter variable names
            setattr(self,k,val) #add attributes to class and #add them into the field
    
    
    def __str__(self, style='gff3'):
        """Create a valid GFF format string representation of this Feature instance
        """
        
        if style == 'gff3':
            orddict = OrderedDict()
            for elem in GFF3_FORMAT_ORDER:
                try:
                    orddict[elem] = self.__dict__[elem]
                except:
                    orddict[elem] = GFF3_UNDEF

            outpt = GFF3_SEP.join(orddict.values()) + GFF3_END
            return outpt
    
    def __repr__(self):
        try:
            assert(self.node_id)
            return self.node_id
        except:
            return super(Feature, self).__repr__()
            #return object.__repr__(self)
    
class FeatureLst:
    """A class which is instantiated with a list containing Feature objects. It automatically
    creates some commentary files and harbours the methods to print them into a gff format file
    or stdout"""
    
    def __init__(self, featurelst=[]):
        
        self.featurelst = featurelst
        self.featuredct = {}
        self.seqregion_dict = {}     
        
        for feature in featurelst:
            refid = feature.__dict__["refid"]
            
            if not refid in self.seqregion_dict.keys():
                self.seqregion_dict[refid] = [feature.__dict__["start"], feature.__dict__["end"]]
            else:
                start = feature.__dict__["start"]
                end = feature.__dict__["end"]
                if  int(start) < int(self.seqregion_dict[refid][0]):
                    self.seqregion_dict[refid][0] = start
                if  int(end) > int(self.seqregion_dict[refid][1]):
                    self.seqregion_dict[refid][1] = end
    
    def generateFeatureDict(self, gff=False):
        """Creates a dict of type contig:feat which lets us access features via the name of the
        contig (node) it is located on."""
        import collections

        featuredct = collections.defaultdict(list) 
        for feat in self.featurelst:
            if "node_id" in feat.__dict__.keys():
                featuredct[feat.__dict__["node_id"]] = feat
        
        if gff:
            for feat in self.featurelst:
                if "seqid" in feat.__dict__.keys():
                    featuredct[feat.__dict__["seqid"]] = feat
        
        return featuredct 
    
    def toFile(self, prefix="./output", suffix="gff3", force_overwrite=None):
        """Creates a gff file by opening a user-specified output path and writing the
        object's instance content into the file line by line. Also writes the constant
        string (defined in the top part of this file) in the beginning and ending lines
        of the file."""
         
        try:    #if file already exists, change name to ~_1 or ~_2 etc. until a non-existing file name found (only if f-flag not set
            i=1
            if not force_overwrite:
                while os.path.isfile(prefix + "." + suffix):
                    if i==1:
                        prefix += "_0"
                    prefix = os.path.join(os.path.dirname(prefix),
                                        os.path.basename(prefix)[0:-(1 + len(str(i-1)))] + '_' + str(i))
                    i+=1
            filepath = prefix + "." + suffix
        except IOError:
            sys.stderr.write("Error while searching for already existing filenames.")
            
        with open(filepath, 'w', encoding='utf-8') as fh:
            fh.write(GFF3_FIRST_LINE)
            
            for refid, vl in self.seqregion_dict.items():   #creates the commentary lines like ##sequence-region NC_017847.1 1 3964912
                fh.write(GFF3_SECOND_LINE.format(refid, vl[0], vl[1]))
                
            for feat in self.featurelst:
                fh.write(feat.__str__(style=suffix))
            fh.write(GFF3_LAST_LINE)
            
    def __str__(self):
        """returns a gff file formatted string of the object's instance"""
        
        outpt  = [GFF3_FIRST_LINE]
                
        for refid, vl in self.seqregion_dict.items():   #creates the lines like ##sequence-region NC_017847.1 1 3964912
            outpt += GFF3_SECOND_LINE.format(refid, vl[0], vl[1])
                        
        outpt += [str(feat) for feat in self.featurelst]
        outpt += [GFF3_LAST_LINE]
        
        return "".join(outpt)
        
   
def coordsFileParser(coordsfile, args):
    """Takes a nucmer or mummer produced *.coords file goes through it line 
    by line. The non-commment lines will be split into list and given to the
    fitting FeatureParser. The output of the latter is then collected into a
    list which finally will be returned as an object of type FeatureLst"""
    
    flst = []                                        #list all Features in a list
    featcount = 0
    coordsfiletype=coordsfile.name.rsplit('.',1)[1]    # auto-detect file ending from coords filepath
    #print(coordsfiletype)
    for linenum, line in enumerate(coordsfile):
        featcount += 1
        if line.startswith('#') or line.isspace():
            featcount -= 1
            continue
        elif coordsfiletype == "1coords" or coordsfiletype == "mcoords":
            flst.append(coordsLineToFeatureParser(line.rstrip().split('\t'), args.source, featcount, args.htflag))
        elif coordsfiletype == "coords":
            valuelst = line.rstrip().split()
            valuelst = list(filter(lambda x: x!='|', valuelst)) #delete list elements which only containing a pipe sign
            flst.append(coordsLineToFeatureParser(valuelst, args.source, featcount, args.htflag))
    
    flst = postProcessingCoords(flst, args)          
    return FeatureLst(featurelst=list(filter(None, flst)))


def coordsLineToFeatureParser(valuelst, source, featcount, htflag, refid=""):
    """Takes a number of values in a list and assigns them their field. The created dictionary
    is transferred to create and return a Feature object. This method therefore needs to by 
    modified if the structure of the input file or GFF format has changed.
    Moreover, this method accepts a featurecount variable for naming as well as threshold
    values which would prohibit the Feature object to be created if the value does not meet
    these values. In these cases None is returned"""
    import re
    spades_pat = re.compile('^NODE_[\d]{0,4}_')
    
    feat = {}
    try:
        if refid == "":
            if "|" in valuelst[11] and len(valuelst[11].split('|')) > 3:
                refid = valuelst[11].split('|')[3]
            else:
                refid = valuelst[11]    
        
        if len(valuelst) == 13:
            feat["seqid"] = refid
            feat["source"] = "coords-file"
            feat["seqtype"] = source
            if int(valuelst[3]) >= int(valuelst[2]):
                feat["start"] = valuelst[0]
                feat["end"] = valuelst[1]
                feat["strand"] = '+'
                feat["qstart"] = valuelst[3]
                feat["qend"] = valuelst[2]
            else:
                feat["start"] = valuelst[0]
                feat["end"] = valuelst[1]
                feat["strand"] = '-'
                feat["qstart"] = valuelst[2]
                feat["qend"] = valuelst[3]
            #feat["score"] = 
            #feat["phase"] = 
            
            feat["alength"] = valuelst[4]
            feat["alength_q"] = valuelst[5]
            feat["total_length_q"] = valuelst[8]
            feat["cov"] =  valuelst[10]
            feat["refid"] = valuelst[11]
            feat["qid"] = valuelst[12]
            
            if spades_pat.match(valuelst[12]): #verify that its a spades standard contig naming scheme
                namesplit = valuelst[12].split('_')
                if htflag:
                    feat["node_id"] = '_'.join(namesplit[0:2]) + '_' + namesplit[-1]
                    feat["ht"] = namesplit[-1]
                else:
                    feat["node_id"] = '_'.join(namesplit[0:2]) 
            else:
                ## if not spades scheme search for first occurence of numbers which is assumed to be the id and delete leading zeros
                k = re.findall(r'\d+', valuelst[12])
                node_id = str(int(k[0]))
                
                if htflag:
                    feat["node_id"] = "CONTIG_" + node_id + '_' + valuelst[12].split('_')[-1]
                    #str(featcount) + '_' + namesplit[-1] #else give auto-increment name
                    feat["ht"] = valuelst[12].split('_')[-1]
                else:
                    feat["node_id"] = "contig_" + str(featcount)
                
            feat["seqattributes"] = "ID={};Name={};idy={};cov={};lenq={};lena={};Is_circular={};mol_type=genomic{};".format(
                                     feat["node_id"],valuelst[12],valuelst[6],valuelst[10],valuelst[8],valuelst[5],"true","DNA")
            
            
            return Feature(**feat)
        
        else:
            raise("Pseudo-Error")
    except:
        if len(valuelst) != 13:
            sys.stderr.write("Error parsing line. The valid coords-file-line needs to have 13 tab separated columns. Skipping.\n")
        else:
            sys.stderr.write("Error parsing line.\n")
            import traceback
            traceback.print_exc()
        return None


def gffFileParser(gfffile, args , featuredct = None):
    flst = []
    featcount = 0
    gfffile.seek(0)
    for linenum, line in enumerate(gfffile):
        featcount += 1
        if line.startswith('#') or line.isspace():
            featcount -= 1
            if line.startswith("##FASTA"):
                break
            else:
                continue
        
        elif len(line.rstrip().split('\t')) == 9:
            flst.append(gffLineToFeatureParser(line.rstrip().split('\t'), args, 
                                               featuredct))
    
    return FeatureLst(featurelst=list(filter(None, flst)))


def checkPresence(gfflst, args):
    """Returns False if the specified filters detected an unwanted feature"""
    if args.gff_type_include:
        if not gfflst[2] in  args.gff_type_include:
            return False
    
    if args.gff_type_exclude:
        if gfflst[2] in args.gff_type_exclude:
            return False
    
    
    if args.gff_attr_include:
        found = False
        for textbit in args.gff_attr_include:
            if textbit in gfflst[8]:
                found = True
                break
        if not found:
            return False
    
    if args.gff_attr_exclude:
        for textbit in args.gff_attr_exclude:
            if textbit in gfflst[8]:
                return False

    return True
    
def gffLineToFeatureParser(gfflst, args, featuredct):
    feat = {}
    
    if featuredct:
        feat["source"] = gfflst[1]
        feat["seqtype"] = gfflst[2]
    
        feat["strand"] = gfflst[6]
        feat["score"] = gfflst[5]
        feat["phase"] = gfflst[7]  
        feat["seqattributes"] = gfflst[8]
    
        node_id = gfflst[0]
        
        if not node_id in featuredct.keys():        #cannot add gff lines if original node id is not found in the FeatureLst dict given
            return None
        
        if not checkPresence(gfflst, args):
            return None
        
        feat["seqid"] = featuredct[node_id].__dict__["seqid"]
        feat["refid"] = feat["seqid"]
        feat["start"] = str(int(featuredct[node_id].__dict__["start"]) + int(gfflst[3]))
        feat["end"] = str(int(featuredct[node_id].__dict__["start"]) + int(gfflst[4]))
          
        if args.only_aligned:    #looks at the aligned range of the contig node, if the Feature lies completely outside this range
                            #it will be skipped
            if int(feat["start"]) > int(featuredct[node_id].__dict__["end"]):
                return None
            if int(feat["end"]) < int(featuredct[node_id].__dict__["start"]):
                return None
    else:
        feat["seqid"] =  gfflst[0]
        feat["source"] = gfflst[1]
        feat["seqtype"] = gfflst[2]
        feat["refid"] = feat["seqid"]
        feat["start"] = gfflst[3]
        feat["end"] = gfflst[4]
        feat["score"] = gfflst[5]
        feat["strand"] = gfflst[6]
        feat["phase"] = gfflst[7]  
        feat["seqattributes"] = gfflst[8]
        
        if not checkPresence(gfflst, args):
            return None
        
        if args.only_aligned:
            raise("Option-Error, can't use only-aligned option when no coords file is specified") #TODO put into argparse as conditional option
           
    return Feature(**feat)


def postProcessingCoords(pre_flst, args):
    """Checks whether a contig is unique, if is not it checks whether
    its parts do overlap and their combined size is not larger than 150%
    of the original. If it is unique it applies the thresholds"""
    
    unique_feat_qids = {}
    pro_flst = []
    
    def isOverlapping(featOLD, featNEW):
        ofs = int(featOLD.start)
        ofe = int(featOLD.end)
        
        if ( (ofs <= int(feat.start) <= ofe) and ( ofe < int(feat.end)) ) or \
           ( (ofs <= int(feat.end)   <= ofe) and ( int(feat.start) < ofs) ):
            return True
    
    def isSameAsFirst(featNEW, featlst):
        
        if feat.node_id == featlst[0].node_id: 
            return True
    
    def isContained(featA, featB):
        ofs = int(featA.start)
        ofe = int(featA.end)
        
        if ( int(featB.start) <= ofs <= ofe <= int(featB.end) or \
                    int(featB.start) <= ofe <= ofs <= int(featB.end) or \
                    int(featB.start) >= ofs >= ofe >= int(featB.end) ):
            return True
    
    def isMoreCovered(featA, featB, factor=20):
        ofcov = float(featA.cov)
        if ofcov >= float(featB.cov) * factor:
            return True
    
    def isQuerySeqContained(featA, featB):
        ofqs = int(featA.qstart)
        ofqe = int(featA.qend)
        
        if ( int(featB.qstart) >=  ofqs >= ofqe >= int(featB.qend) or \
                    int(featB.qstart) <=  ofqs <= ofqe <= int(featB.qend) ):
            return True
    
    if args.unique_requirement:
        for feat in pre_flst:
            
            feat.overlapflag = False
            
            if feat.qid in unique_feat_qids:
                f = unique_feat_qids[feat.qid]
                
                if isContained(feat, f):
                    pro_flst.append(None)
                    continue               
                
                elif isOverlapping(f, feat) or isSameAsFirst(feat, pro_flst):
                    feat.overlapflag = True
                    f.overlapflag = True
                    former_id = f.seqattributes.split(';',1)[0].split('=')[1]
                    new_id = former_id + '*'
                
                elif isContained(f, feat):
                    new_id = f.node_id
                    feat.overlapflag = f.overlapflag
                    
                    for i, o in enumerate(pro_flst[::-1]):
                        if o and o.qid == feat.qid:
                            #del pro_flst[-i-1]
                            pro_flst[-i-1] = None
                            break
                
                else:
                    former_id = f.seqattributes.split(';',1)[0].split('=')[1]
                    new_id = former_id + '^'
                
                feat.seqattributes = ';'.join(['ID=' + new_id , feat.seqattributes.split(';',1)[1]])
                feat.node_id = new_id
            
            pro_flst.append(feat)
            unique_feat_qids[feat.qid] = feat
    
                    
#     if args.unique_requirement:
#         for feat in pre_flst:
#             #check uniqueness and criteria to keep although not unique
#             if feat.qid in unique_feat_qids:
#                 f = unique_feat_qids[feat.qid]
# 
#                 ofs = int(f.start)
#                 ofe = int(f.end)
#                 ofqs = int(f.qstart)
#                 ofqe = int(f.qend)
#                 
#                 
#                 #only overlaps no contaiment allowed, assume circularity, therefore allow non-unique at the beginning and end
#                 if ( (ofs <= int(feat.start) <= ofe) and ( ofe < int(feat.end)) ) or \
#                    ( (ofs <= int(feat.end)   <= ofe) and ( int(feat.start) < ofs) ) or \
#                     feat.node_id == pro_flst[0].node_id or \
#                     float(feat.cov) + float(f.cov) > 2* args.coverage_threshold:
#                     
#                     feat.overlapflag = True
#                     f.overlapflag = True
#                     #print("overlap1")
#                     former_id = f.seqattributes.split(';',1)[0].split('=')[1]
#                     new_id = former_id + '*'
#                     feat.seqattributes = ';'.join(['ID=' + new_id , feat.seqattributes.split(';',1)[1]])
#                     feat.node_id = new_id
#                     pro_flst.append(feat)
#                 
#                     unique_feat_qids[feat.qid] = feat
#                 
#                 #if the old one is contained within the new alignment, then delete old (replace with 'None')
#                 elif ( int(feat.start) <= ofs <= ofe <= int(feat.end) or \
#                     int(feat.start) <= ofe <= ofs <= int(feat.end) or \
#                     int(feat.start) >= ofs >= ofe >= int(feat.end) ) and float(f.cov) < args.coverage_threshold:
#                     feat.node_id = f.node_id
#                     feat.overlapflag = f.overlapflag
#                     
#                     for i, o in enumerate(pro_flst[::-1]):
#                         if o and o.qid == feat.qid:
#                             #del pro_flst[-i-1]
#                             pro_flst[-i-1] = None
#                             break
#                     
#                     pro_flst.append(feat)
#                     unique_feat_qids[feat.qid] = feat         
#                    
#                 elif ( ofs <= int(feat.start) <= int(feat.end) <= ofe or \
#                     ofs <= int(feat.end) <= int(feat.start) <= ofe or \
#                     ofs >= int(feat.start) >= int(feat.end) >= ofe ) and float(feat.cov) < args.coverage_threshold:
#                     pro_flst.append(None)
#                     print("The new feature was deleted,", feat.node_id, "its start and end are in between the one already existing")
#                 
#                 # if the mapped part of the same contig is contained in the new one and the old one's maplength is only 5% or less, then delete the old one 
#                 elif ( int(feat.qstart) >=  ofqs >= ofqe >= int(feat.qend) or \
#                     int(feat.qstart) <=  ofqs <= ofqe <= int(feat.qend) ):
#                         if not abs(ofqs-ofqe) <= (abs(int(feat.qstart) - int(feat.qend)) * 0.05):
#                             feat.overlapflag = False
#                             f.overlapflag = False
#                             #print("overlap1")
#                             former_id = f.seqattributes.split(';',1)[0].split('=')[1]
#                             new_id = former_id + '*'
#                             feat.seqattributes = ';'.join(['ID=' + new_id , feat.seqattributes.split(';',1)[1]])
#                             feat.node_id = new_id
#                             
#                             pro_flst.append(feat)
#                             unique_feat_qids[feat.qid] = feat   
#                             
#                         else:
#                         
#                             feat.overlapflag = False
#                             f.overlapflag = False
#                             #print("overlap1")
#                             former_id = f.seqattributes.split(';',1)[0].split('=')[1]
#                             new_id = former_id
#                             feat.seqattributes = ';'.join(['ID=' + new_id , feat.seqattributes.split(';',1)[1]])
#                             feat.node_id = new_id
#                             
#                             for i, o in enumerate(pro_flst[::-1]):
#                                 if o and o.qid == feat.qid:
#                                     pro_flst[-i-1] = None
#                                     break
#                             
#                             pro_flst.append(feat)
#                             unique_feat_qids[feat.qid] = feat   
#                 
#                 elif ( int(ofqs) >=  int(feat.qstart) >= int(feat.qend) >= ofqe or \
#                     int(ofqs) <=  int(feat.qstart) <= int(feat.qend) <= ofqe ) and \
#                     abs(int(feat.qstart) - int(feat.qend)) <= (abs(ofqs-ofqe)* 0.05):
#                 
#                     pro_flst.append(None)
#                     #print("The new feature was deleted, it smaller 5% than the one already in the list and its start and end are already covered in the bigger one")
#                 
#                 else:
#                     #print(feat.qstart, ofqs, ofqe, feat.qend)
#                     #print(abs(int(feat.qstart) - int(feat.qend)), '<=', (abs(ofqs-ofqe)* 0.05))
#                     print('Strange case, revisit the mccords file and compare')
#                      
#             else:
#                 #print('NO_OVER',feat.qid)
#                 feat.overlapflag = False
#                 unique_feat_qids[feat.qid] = feat #in case an overlap occurs the new one replaces the old entry (to allow multiple splits in the contig
#                 pro_flst.append(feat)
#                 #print(pro_flst)
    else:
        for feat in pre_flst:
            setattr(feat, 'overlapflag', False)
        pro_flst = pre_flst
    
    #print('PRE',pre_flst)
    #print('PRO',pro_flst)
    #print()
    final_flst = []
    for feat in pro_flst:
        #if feat:
        #    print(feat.node_id, feat.cov, feat.alength, feat.overlapflag)
        if feat and not feat.overlapflag \
                and args.coverage_threshold <= float(feat.cov) \
                and args.length_threshold <= int(feat.alength) :
                #and int(feat.alength_q) / int(feat.total_length_q) >= args.length_ratio_threshold:
                
            final_flst.append(feat)
        elif feat and feat.overlapflag:
            final_flst.append(feat)
        else:
            #if feat:
            #    print('filtered', feat.node_id)
            final_flst.append(None)
    #print('FIN',final_flst)    
    return final_flst


def parseContigBreaks(featLst):
    """Expects an ordered FeatureLst of Contigs for a single seqid,finds the uncovered regions
    and saves theses sequences into a feature file (gff)"""
    prelst = []
    ID = 1
    
    for i, feat in enumerate(featLst.featurelst):
        
        if i == 0:
            #a pseudo feature from which the uncovered region at the start of the sequence starts
            start_pseudo_feature = {"start": "0", "end": "1" }
            feat0 = Feature(**start_pseudo_feature)
            tmp_breaklength = min(int(feat.end), int(feat.start)) - max(int(feat0.start), int(feat0.end)) + 1
            if tmp_breaklength < 0: #skip if genome is covered right from position 1
                feat0 = feat
                continue
        
        breaklength = min(int(feat.end), int(feat.start)) - max(int(feat0.start), int(feat0.end)) + 1
        
        if breaklength < 0:
            name_short = "OvBreak "
            ns = "obr"
            breaklength *=-1
            end = str(max(int(feat0.start), int(feat0.end)))
            start = str(min(int(feat.end), int(feat.start)))
            custom_seqtype = "intercontig_overl_region"

        else:
            name_short = "Break "
            ns = "br"
            start = str(max(int(feat0.start), int(feat0.end)) + 1)
            end = str(min(int(feat.end), int(feat.start)) - 1)
            custom_seqtype = "intercontig_region"
        
        break_feat = {"seqid": feat.seqid,
                      "refid": feat.refid,
                      "source": feat.source, 
                      "seqtype": custom_seqtype,
                      "start": start, 
                      "end": end,
                      "strand": '+',
                      "score": '.',
                      "phase": '.',
                      "seqattributes": "ID={};Name={};idy={};cov={};lenq={};lena={};Is_circular={};mol_type=genomic{};".format(
                                     ns + str(ID), name_short + str(ID), "NA", "NA", breaklength, breaklength, "true", "DNA")}
             
        ID += 1              
        prelst.append(Feature(**break_feat))
        feat0 = feat
    
    #if circular features are allowed
    if prelst[0].start == '0':
        prelst[0].start = featLst[-1].end #make the first break (if it starts at pos 0)
        
    #TODO: How about the break at the end of the featurelst, i.e. the uncovered part of the sequence
    #starting from the end of the last feature to the end of the whole sequence?
    #Need to parameterize the total seq length somehow...            
    
    return FeatureLst(featurelst=list(filter(None, prelst)))


if __name__ == '__main__':

    args = parse_args()
    # With a coords file go ahead and parse it into a list of features
    if args.infilepath:
        featLst = coordsFileParser(args.infilepath, args) #args.source, args.coverage_threshold, args.length_threshold, args.unique_requirement

        if args.breaks_only:
            featLst = parseContigBreaks(featLst)
            
        if args.outfilepath == "":
            print(featLst)
        else:
            featLst.toFile(prefix=args.outfilepath, force_overwrite=args.force_ow)
    
    # With an additional or independent gff file go ahead and parse it
    if args.gffinfilepath:
        
        if args.infilepath:
            gff_adjust_featLst = gffFileParser(args.gffinfilepath, args, featuredct = featLst.generateFeatureDict())
            
            if args.gffoutfilepath == "":
                print(gff_adjust_featLst)
            else:
                gff_adjust_featLst.toFile(prefix=args.gffoutfilepath, force_overwrite=args.force_ow)
                        
        else:
            featLst = gffFileParser(args.gffinfilepath, args)
            
            if args.mappingGFFpath:
                overlapping_features = []
                nonoverlapping_features = []
                featLst2 = gffFileParser(args.mappingGFFpath, args)
                for f in featLst.featurelst:
                    flag = False
                    #print(f.start, '..', f.end)
                    for mf in featLst2.featurelst:
                        #print('\t', mf.start, '..', mf.end)
                        #{ ... [[ ... } ... ]] right overlap
                        #[[ ... { ... ]] ... } left overlap
                        #[[ ... { ... } ... ]] inclusion (other case of inclusion is already handled by first two) 
                        if (int(f.start) <= int(mf.start) <= int(f.end)) or \
                           (int(f.start) <= int(mf.end) <= int(f.end)) or \
                           ((int(mf.start) <= int(f.start) <= int(mf.end)) and (int(mf.start) <= int(f.end) <= int(mf.end))):          
                            overlapping_features.append(mf)
                            #print(f.start, "<=", mf.start, "<=", f.end)
                            #print(f.start, "<=", mf.end, "<=", f.end)
                            flag = True
                        else:
                            if flag:
                                break
                    if not flag:
                        nonoverlapping_features.append(f)
                mappedFeats = FeatureLst(featurelst=overlapping_features)
                unmappedFeats = FeatureLst(featurelst=nonoverlapping_features)
                #print(mappedFeats)
                print(unmappedFeats)
                featLst = mappedFeats
                
            if args.gffoutfilepath == "":
                print(featLst)
            else:
                featLst.toFile(prefix=args.gffoutfilepath, force_overwrite=args.force_ow)
    
# Description of the Format
# 
# GFF3 files are nine-column, tab-delimited, plain text files. Literal use of tab, newline, carriage return, the percent (%) sign,
# and control characters must be encoded using RFC 3986 Percent-Encoding; no other characters may be encoded.
# Backslash and other ad-hoc escaping conventions that have been added to the GFF format are not allowed.
# The file contents may include any character in the set supported by the operating environment, although for portability with
# other systems, use of Latin-1 or Unicode are recommended.
# 
#     tab (%09)
#     newline (%0A)
#     carriage return (%0D)
#     % percent (%25)
#     control characters (%00 through %1F, %7F) 
# 
# In addition, the following characters have reserved meanings in column 9 and must be escaped when used in other contexts:
# 
#     ; semicolon (%3B)
#     = equals (%3D)
#     & ampersand (%26)
#     , comma (%2C) 
# 
# Note that unescaped spaces are allowed within fields, meaning that parsers must split on tabs, not spaces. 
# Use of the "+" (plus) character to encode spaces is depracated from early versions of the spec and is no longer allowed.
# 
# Undefined fields are replaced with the "." character, as described in the original GFF spec.
# 
# Column 1: "seqid"
#     The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain any characters, 
#      but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain unescaped 
#      whitespace and must not begin with an unescaped ">". 
# Column 2: "source"
#     The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type creating a new composite type that is a subclass of the type in the type column. 
# 
# Column 3: "type"
#     The type of the feature (previously called the "method"). This is constrained to be either: 
#     (a)a term from the "lite" version of the Sequence Ontology - SOFA,
#     (b)a term from the full Sequence Ontology - it must be an is_a child of sequence_feature (SO:0000110) or 
#     (c)a SOFA or SO accession number. The latter alternative is distinguished using the syntax SO:000000. 
# 
#     (d) Look up at: https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo
# 
# Columns 4 & 5: "start" and "end"
# 
#     The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative to the landmark given in column one. Start is always less than or equal to end. For features that cross the origin of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to be less than or equal to end is satisfied by making end = the position of the end + the length of the landmark feature.
# 
#     For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the indicated base in the direction of the landmark.
# Column 6: "score"
#     The score of the feature, a floating point number. As in earlier versions of the format, the semantics of the score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and that P-values be used for ab initio gene prediction features. 
# Column 7: "strand"
#     The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown. 
# Column 8: "phase"
# 
#     For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon. In other words, a phase of "0" indicates that the next codon begins at the first base of the region described by the current line, a phase of "1" indicates that the next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins at the third base of this region. This is NOT to be confused with the frame, which is simply start modulo 3.
# 
#     For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted from the end field.
# 
#     The phase is REQUIRED for all CDS features.
# Column 9: "attributes"
# 
#     A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape. Attribute values do not need to be and should not be quoted. The quotes should be included as part of the value by parsers and not stripped.
# 
#     These tags have predefined meanings:
# 
#     ID
#         Indicates the ID of the feature. IDs for each feature must be unique within the scope of the GFF file. In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines. All lines that share an ID collectively represent a single feature. 
#     Name
#         Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file. 
#     Alias
#         A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file. 
#     Parent
#         Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, an so forth. A feature may have multiple parents. Parent can *only* be used to indicate a partof relationship. 
#     Target
#         Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20. 
#     Gap
#         The alignment of the feature to the target if the two are not collinear (e.g. contain gaps). The alignment format is taken from the CIGAR format described in the Exonerate documentation. See "THE GAP ATTRIBUTE" for a description of this format. 
#     Derives_from
#         Used to disambiguate the relationship between one feature and another when the relationship is a temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See "PATHOLOGICAL CASES" for further discussion. 
#     Note
#         A free text note. 
#     Dbxref
#         A database cross reference. See the section "Ontology Associations and Db Cross References" for details on the format. 
#     Ontology_term
#         A cross reference to an ontology term. See the section "Ontology Associations and Db Cross References" for details. 
#     Is_circular
#         A flag to indicate whether a feature is circular. See extended discussion below. 
# 
#     Multiple attributes of the same type are indicated by separating the values with the comma "," character, as in:
# 
#     Parent=AF2312,AB2812,abc-3
# 
#     In addition to Parent, the Alias, Note, Dbxref and Ontology_term attributes can have multiple values.
# 
#     Note that attribute names are case sensitive. "Parent" is not the same as "parent".
# 
#     All attributes that begin with an uppercase letter are reserved for later use. Attributes that begin with a lowercase letter can be used freely by applications.
