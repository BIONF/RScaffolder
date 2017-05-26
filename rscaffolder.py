'''
Created on Feb 21, 2017

@author: Ba1
'''
import argparse
import sys
from os import listdir, getcwd, makedirs
from os.path import join, isfile, isdir, basename, dirname, abspath, expanduser
import subprocess
from Bio import SeqIO
import pandas as pd
import re
from scripts.pyGFA import parseGFA
from scripts.coords2gff3 import Feature
import networkx as nx
import graphviz as gv
import functools


def parse_args(args):
    parser = argparse.ArgumentParser(description="Multi-Reference-Genome-Assisted Scaffolding of target assembly")
    
    parser.add_argument('-i', '--contigs_filepath', dest='contigs_path', metavar='</path/to/contigfile>', type=argparse.FileType('rt'), required=True,
                        help='path to an contigs fasta file')

    parser.add_argument('-rg', '--reference-genomes-directory', dest='refgenomesdir', metavar='<path>', type=str,
                        default='./references', 
                        help='specify a directory containing the reference genomes (multi-) fasta nucleotide sequence files'),
    
    parser.add_argument('-spdir', '--assembly-graph-dir', dest='assemblygraph', metavar='<path>', type=str, #type=argparse.FileType('rt'), 
                        help='specify the spades output directory to parse the assembly graph file in GFA 1.0 format'),

    parser.add_argument('-o', '--outdir', dest='outdirpath', metavar='</path/to/out/>', type=str,
                        default=getcwd(), 
                        help='path to output file prefix'),
                        
    parser.add_argument('-DEBUGexclude', dest='exclude',  metavar='<str>',type=str, default=False, nargs='+',
                   help='References with sequence headers containing any of the specified string will get deleted.')
    
    parser.add_argument('-f', '--force', dest='force_flag', metavar='<force overwrite flag>', action='store_const', const=1, default=0,
                   help='if set, the files in the specified output directory will be overwritten')

    return parser, parser.parse_args()


def getRefFilePaths(refgenomesdir):
    """Check reference directory path and return list of genome files
    """
    assert isdir(refgenomesdir), "{} is not a directory"
    
    accepted_suffices = (".fna", ".fa", ".fasta")
    reference_genomes = []
    
    for f in listdir(refgenomesdir):
        if isfile(join(refgenomesdir, f)) and '.' + f.rsplit('.', 1)[-1] in accepted_suffices:
            reference_genomes.append(join(refgenomesdir,f))
            
    assert len(reference_genomes) > 0, "{} does not contain any sequence files in FASTA Format. Accepted files suffixes: {}".format(refgenomesdir, accepted_suffices)
    
    return reference_genomes


def parseContigs(contigs_filepath, per_base_coverages=False, pattern=False):
    """Parse the information of the nodes, create a statistics file with coverage information
    indicate which contigs are considered for the scaffolding process"""

    seqs = SeqIO.parse(open(contigs_filepath, "r"), "fasta")
    
    seq_vallist = []
    index_list = []
    recordlist = []
        
    for record in seqs:
        recordlist.append(record)       #Biopython Seq object
        m = re.findall(pattern, record.id)
        if m:
            if m[0][0] in index_list:
                continue
            index_list.append(m[0][0])
            value_dict = {'length': int(m[0][1]), 'cov': float(m[0][2])}
            seq_vallist.append(value_dict)
        else:
            m = re.findall(r'\d+', record.id)
            m = str(int(m[0][0]))
            if m in index_list:
                continue
            index_list.append(m)
            import numpy as np
            value_dict = {'length': int(len(record)), 'cov': np.nan}
            seq_vallist.append(value_dict)
            
            
        df = pd.DataFrame(seq_vallist, index=index_list)
    
    bools = df['length'] > 0
    df = df.assign(used=bools.values)

    if 'cov' in df.columns:
        try:
            cov_median = df['cov'].median()
            mad = abs(df['cov'] - cov_median).median()
 
            df.loc[:,'multipl'] = pd.Series(df['cov']/df['cov'].median(), index=df.index)
        
            for index, row in df.iterrows():
                if abs(round(row['multipl']) - row['multipl']) < 2* mad:
                    df.set_value(index,'multipl',round(row['multipl']))
                else:
                    df.set_value(index,'multipl', np.nan)
        except:
            df['multipl'] = np.nan
    
    
    return recordlist, df


def filterAndStoreContigs(recordlist, df, pattern=False, key='lr'): 
    record_dict = dict()
    
    for seqrec in recordlist:
        m = re.findall(pattern, seqrec.id)
        if m:
            record_dict[seqrec.id.split('_')[1]] = seqrec
        else:
            k = re.findall(r'\d+', seqrec.id)
            record_dict[str(int(k[0]))] = seqrec 
    
    # OUTPUT #
    filtered_contigs_path = join(args.outdirpath,'contigs_filtered.fasta')
    
    ##Filter by ###
    if isinstance(key, (int, float)): # absolute value
        print("Using absolute value as threshold for minimum contig size", key)
        filtered_records_dict = {rid:record for rid,record in record_dict.items() if len(record.seq) >= key}
    
    if key == 'lr': # longest repetitve 
        repeats_df = df[df["multipl"] > 1.0]
        longest_repeat_length = repeats_df.max().loc['length']
        filtered_records_dict = {rid:record for rid,record in record_dict.items() if len(record.seq) >= longest_repeat_length + 1}

    #######################################################################################
            
    SeqIO.write(filtered_records_dict.values(), filtered_contigs_path, "fasta")    
    
    return record_dict, filtered_records_dict, filtered_contigs_path


def parseSpadesAssemblyGraph(spadesdir, records_dict):
        
    import scripts.pyGFA
    
    G = scripts.pyGFA.parseGFA(join(spadesdir,'assembly_graph.gfa'))
    
    edge2seg_map = dict()
    edges_list = []
    seqs = SeqIO.parse(open(join(spadesdir, 'assembly_graph.fastg'), "r"), "fasta")    
    
    for i, record in enumerate(seqs):
        if i % 2 == 0: #skip over the inverse edges
            edges_list.append(record)       #Biopython Seq object

    for segid, seq in G.segments.items():
        for record in edges_list:
            if seq == record.seq:
                pattern = re.compile('^EDGE_([\d]+)_length_([\d]+)_cov_([\d*\.\d*]+).*')
                m = re.findall(pattern, record.id)
                edge2seg_map[m[0][0]] = segid
                break
    
    def otransl(b):
        if isinstance(b, bool):
            if b:
                return '+'
            else:
                return '-'
        if b in ("+","-"):
            if b == "+":
                return True
            else:
                return False
        
    current = False # only added for pydev to stfu
    
    with open(join(spadesdir,'contigs.paths'), "r") as cpf:
        flag = False
        for __, line in enumerate(cpf):
            if flag == True:
                path = []
                for x in line.strip().split(','):
                    path.append((edge2seg_map[x[:-1]], otransl(x[-1])))

                records_dict[current].__dict__["path"] = path
                
                #if path[0][1]:
                records_dict[current].__dict__["head"] = (path[0][0], not path[0][1])
                #else:
                #    records_dict[current].__dict__["head"] = (path[0][0], path[0][1])
                
                #if path[-1][1]:
                records_dict[current].__dict__["tail"] = (path[-1][0],path[-1][1])
                #else:
                #    records_dict[current].__dict__["tail"] = (path[-1][0], not path[-1][1])

            if line.startswith('NODE_') and not line.strip().endswith("'"):
                pattern = re.compile('^NODE_([\d]+)_length_([\d]+)_cov_([\d*\.\d*]+).*')
                m = re.findall(pattern, line.strip())
                current = m[0][0]
                flag = True
                continue
            flag = False
    
    ordering = [str(i) for i in sorted([int(rid) for rid in records_dict])]
    
    with open(join(args.outdirpath, "gfa_contigs.paths"), 'w') as cpaths:
        for rid in ordering:
            rec = records_dict[rid]
            cpaths.write(rec.name + '\n' + str(rec.path) + '\n')
    
    
    seg2edge_map = {v: k for k, v in edge2seg_map.items()}
    node_edges_still_uncon = {}
    
    #print(seg2edge_map)
    #print('###########')
    
    def parseLinksToList(itself, link_dict, records_dict):
        """
        For each link in the link_dict associated with an edge
        go through the dictionary of contigs and collect those where
        this link links to a contigs' edge paths' head or tail edge
        """
        
        linklist = []
        for link in link_dict: #ignoring the CIGAR value of each link key
            for rid, rec in records_dict.items():
                if not len(rec.path) == 1:
                    if link == rec.path[0]:
                        linklist.append([itself , 'NODE_' + rid + '_H'])
                        
                    if link == rec.path[-1]:
                        linklist.append([itself , 'NODE_' + rid + '_T'])
                        
                else:
                    if link == rec.path[0]:
                        linklist.append([itself , 'NODE_' + rid + '_H'])
                        #continue
                    if link == (rec.path[-1][0], not rec.path[-1][1]) :
                        linklist.append([itself , 'NODE_' + rid + '_T'])
                        #continue
                    
                    
        return linklist
    
    cons = []
    unconnected = []

    for rid in ordering:
        rec = records_dict[rid]
        h = rec.head
        t = rec.tail
        not_h = (rec.head[0], not rec.head[1])
        not_t = (rec.tail[0], not rec.tail[1])
        
        node_edges_still_uncon[rid] = [h,t,not_h,not_t]
        
        for link, connectedlinks in G.links.copy().items():
            if link == h:
                node_edges_still_uncon[rid][0] = False
                #print(rid, 'HEAD', link, '[', seg2edge_map[link[0]] +  otransl(link[1]), ']', '->' , connectedlinks)
                if not connectedlinks:
                    unconnected.append(rid + '_H ')
                else:
                    cons.append(parseLinksToList('NODE_' + rid + '_H', connectedlinks, records_dict))
                    #print("JUST APPENDED,", parseLinksToList('NODE_' + rid + '_H', connectedlinks, records_dict))
            elif link == t:
                node_edges_still_uncon[rid][1] = False
                #print(rid, 'TAIL', link, '[', seg2edge_map[link[0]] +  otransl(link[1]), ']', '->' , connectedlinks)
                if not connectedlinks:
                    unconnected.append(rid + '_T ')
                else:
                    cons.append(parseLinksToList('NODE_' + rid + '_T', connectedlinks, records_dict))
                    #print("JUST APPENDED,", parseLinksToList('NODE_' + rid + '_T', connectedlinks, records_dict))
            elif link == not_h:
                node_edges_still_uncon[rid][2] = False
                #print(rid, 'CONTIG_DIRECTION', 'HEAD', link, '[', seg2edge_map[link[0]] +  otransl(link[1]) ,']', '->' , connectedlinks)
            elif link == not_t:
                node_edges_still_uncon[rid][3] = False
                #print(rid, 'CONTIG_DIRECTION', 'TAIL', link, '[', seg2edge_map[link[0]] +  otransl(link[1]) ,']', '->' , connectedlinks)
#     for rid, edges in node_edges_still_uncon.items():
#         if any(edges):
#             print(rid, ':', edges)

    con_tuples_dict = dict()

    con_tuples = []
    
    for i in cons:
        con_tuples.extend([item for sublist in i for item in sublist])
    
    con_tuples = list(zip(*2*[iter(con_tuples)]))
    
    i = 0
    for tup in con_tuples:
        con_tuples_dict[i] = tup
        i += 1
    
    return con_tuples_dict, unconnected


def filterConnections(contig_connections, unconnected, filtered_records_dict):
    
    
    ccon = contig_connections.copy()
    
    valid_contigs = []
    for cid in filtered_records_dict:
        valid_contigs.extend(['NODE_{}_{}'.format(cid, 'H'),'NODE_{}_{}'.format(cid, 'T')])
    
    for i, con in contig_connections.items():
        if not all(c in valid_contigs for c in con) or con[0]==con[1]:
            del ccon[i]
    
    uncon = [u for u in unconnected if 'NODE_' + u.strip() in valid_contigs]
    
    return ccon, uncon


def extractHeadTail(contigs_path, outdirpath, ht_max_length=2000, margin_length=150):
    
    cmd = [PYTHON_PATH, join(dirname(__file__), 'scripts' ,'fasta_headtail_extract.py'), 
           '-i', contigs_path,
            '-min', str(2 + 2 * int(margin_length)), 
            '-l', str(ht_max_length), 
            '-ms', str(int(margin_length)), 
            '-o', join(outdirpath, basename(contigs_path).split('.')[0])]
    
    exit_code = subprocess.call(' '.join(cmd), shell=True)
    
    if exit_code > 0 :
        print('Error while trying to extract the head and tail sequences from fasta file:', contigs_path)
        raise
    else:
        outpath = join(outdirpath, basename(contigs_path).split('.')[0] + '_HT.fasta')
    
    return outpath
    

def nucmerMapping(contigs_path, refgenomesdir, outdir):
    cmd = [PYTHON_PATH, join(dirname(__file__), 'scripts' , 'nucmer_mapping.py'), 
           '-i', contigs_path, '-o', outdir,
           '-rd', refgenomesdir]

    #print(' '.join(cmd))
    
    exit_code = subprocess.call(' '.join(cmd), shell=True)
    
    if exit_code > 0 :
        print('Error while trying to map', contigs_path,'to reference genome sequences in', refgenomesdir)
        raise IOError()
    else:
        return


def getMaporder(refmapdir, cov = 70, alen_thresh = 1300, unique = True, htflag = False):
    '''get the mapping order of the contigs which meet the 
    thresholds for each reference in the reference directory'''
        
    file_ending_glob = '*_dnadiff.mcoords'
    
    import scripts.coords2gff3 as c2g
    import glob
    
    args2 = argparse.Namespace()
    args2.__dict__.update({'coverage_threshold': cov,
                           'length_threshold' : alen_thresh,
                           'unique_requirement' : unique,
                           'source' : 'contig',
                           'htflag' : htflag})
    
    perms = []
        
    for ref_map_path in glob.glob(join(refmapdir, file_ending_glob)):
        args2.__dict__['infilepath'] = open(ref_map_path ,'r') 
        feats = c2g.coordsFileParser(args2.__dict__['infilepath'], args2)
        
        perm = dict()
        
        for feat in feats.featurelst:
            if not feat.seqid in perm:
                perm[feat.seqid] = []
            #if vis_only:
            #    if htflag:
            #        feat_repr = feat.strand + feat.__repr__().split('_')[1]
            #    else:
            #        feat_repr = feat.strand + feat.__repr__().split('_')[1] + feat.__repr__().split('_')[-1]
            #    perm[feat.seqid].append(feat_repr)
            #else:
            perm[feat.seqid].append(feat)
        perms.append(perm)
    
    return perms


def filterPerms(perms, filtered_records_dict):
    
    valid_contigs = []
    for cid in filtered_records_dict:
        valid_contigs.extend(['NODE_{}_{}'.format(cid, 'H'),'NODE_{}_{}'.format(cid, 'T')])
    
    filtered_perms = []
    for refgenome in perms:
        seqdict = {}
        for seq, perm in refgenome.items():
            seqdict[seq] = [c for c in perm if c.__repr__() in valid_contigs]
        filtered_perms.append(seqdict)
    
    return filtered_perms


def initGraph(cons=False):
    
    def getShortLabel(label):
        return label.replace('NODE_','') + ' '
    
    def getOppositeLabel(label):
        if label.endswith('H '):
            return label.replace('_H ', '_T ')
        else:
            return label.replace('_T ', '_H ')
    
    G = nx.Graph()
    
    if not cons:
        print("Information from assembly graph not used")
        return G
    
    for con in cons.values():
        G.add_edge(getShortLabel(con[0]), getShortLabel(con[1]), weight=1, ref=('ASSEMBLY',), edgetype=1)
        
    for node in G.nodes():
        if not G.has_edge(node, getOppositeLabel(node)):
            G.add_edge(node, getOppositeLabel(node), weight=100, ref=('DEFAULT',), edgetype=1)

    return G


def toGraph(perms, dist_thresh, G=False, htflag = False):
    
    if not G:
        G = nx.Graph()
        
    from itertools import cycle
    
    def getNodeLabel(feat, htflag):
        """Gives out a 2-tuple, with the long version label plus orientation, 
        id and head/tail info on the contig, as well as a shorter version
        leaving out the head and tail info"""
        if htflag:
            node_label = feat.node_id
            node_label = '_'.join(node_label.split('_')[1:]) + ' '  #replace NODE_ or CONTIG_ in front of contig label
            lab = node_label.rsplit('_')[0]
        else:
            node_label = feat.strand + '_' + feat.node_id
            node_label = '_'.join(node_label.split('_')[1:]) + ' '
            lab = node_label = node_label.replace('_', '')
        
        return node_label, lab
    
    def getInvLabel(label):
        if label.startswith('+'):
            label = '-' + label[1:]
        else:
            label = '+' + label[1:]
            
        return label
       
    adj_lists = {}
        
    for permdic in perms:
        for k, perm in permdic.items():
            adj_lists[k] = []
            
            all_connected = cycle(perm)
            next(all_connected) #to start with element two
            for mcon in perm:
                inverse = False
                #if isinstance(mcon, Feature):
                
                mcon_label, lab = getNodeLabel(mcon, htflag=htflag)
                
                if not htflag and G.has_node(getInvLabel(mcon_label)):
                        mcon_label = lab = getInvLabel(mcon_label)
                        inverse = True

                if not G.has_node(mcon_label):
                    G.add_node(mcon_label)
               
                nextcon = next(all_connected)
                nextcon_label, nlab = getNodeLabel(nextcon, htflag=htflag )
                       
                if inverse:
                    nextcon_label = nlab = getInvLabel(nextcon_label)
                    
                distance = (mcon_label, nextcon_label, abs(int(nextcon.start) - int(mcon.end)))                                      
                adj_lists[k].append(distance)
                                
                if nlab == lab:
                    if mcon_label == nextcon_label: # if this is something of the style +1_H +1_H +1_H +1_H +1_T , 
                        # ignore the first mappings of the head
                        continue
                    #insert here total length threshold for distance between head and tail in comparison to total
                    if G.has_edge(mcon_label, nextcon_label):
                        G[mcon_label][nextcon_label]["weight"] += 1
                        G[mcon_label][nextcon_label]["ref"] += (k,)
                    else:
                        G.add_edge(mcon_label, nextcon_label, weight=1, ref=(k,), edgetype=1)
                
                elif (abs(int(nextcon.start) - int(mcon.end)) < dist_thresh):# or nextcon.overlapflag: #threshold + margin
                    if G.has_edge(mcon_label, nextcon_label):
                        G[mcon_label][nextcon_label]["weight"] += 1.0
                        G[mcon_label][nextcon_label]["ref"] += (k,)
                    else:
                        G.add_edge(mcon_label, nextcon_label, weight=1.0, ref=(k,), edgetype=1)
#                 else:
#                     print('No edge drawn. Too far apart: ', mcon_label, mcon.end, 'and', nextcon_label, nextcon.start, \
#                          abs(int(nextcon.start) - int(mcon.end)), 'bp')
#     toremove = []
#     for e in G.edges_iter(data=True):
#         e[2]['weight'] -= 1
#         if e[2]['weight'] < 1:
#             toremove.append((e[0],e[1]))
#     for e in toremove:
#         G.remove_edge(*e)

    
    for permdic in perms:                           #second iteration, now going for contigs spanning two
        for k, perm in permdic.items():
            #print(k)
            
            all_connected = cycle(perm)
            all_connected2 = cycle(perm)
            all_connected3 = cycle(perm)
            
            [next(all_connected) for _ in range(1)] #to later start with element two
            [next(all_connected2) for _ in range(2)] #to later start with element three
            [next(all_connected3) for _ in range(3)] #to later start with element four
            
            for mcon in perm:
                mcon_label, lab = getNodeLabel(mcon, htflag= htflag)
                
                nextcon = next(all_connected)
                nextcon2 = next(all_connected2)
                nextcon3 = next(all_connected3)
                
                nextcon_label, nlab = getNodeLabel(nextcon, htflag= htflag)
                nextcon2_label, nlab2 = getNodeLabel(nextcon2, htflag= htflag)
                nextcon3_label, nlab3 = getNodeLabel(nextcon3, htflag= htflag)
                
                if not htflag and not G.has_node(mcon_label):
                    mcon_label = lab = getInvLabel(mcon_label)
                    nextcon_label = nlab = getInvLabel(nextcon_label)
                    nextcon2_label = nlab2 = getInvLabel(nextcon2_label)
                    nextcon3_label = nlab3 = getInvLabel(nextcon3_label)
                 
                
                if G.has_edge(mcon_label, nextcon_label) and k in G[mcon_label][nextcon_label]["ref"]:
                    if G.has_edge(nextcon_label, nextcon2_label) and k in G[nextcon_label][nextcon2_label]["ref"]:
                        if G.has_edge(nextcon2_label, nextcon3_label) and k in G[nextcon2_label][nextcon3_label]["ref"]:
                            if len(set([lab,nlab,nlab2,nlab3])) > 2: 
                                if "color" in G[mcon_label] and "color" in G[nextcon3_label]:
                                    G.node[mcon_label]["color"].update(nextcon3_label)
                                    G.node[nextcon3_label]["color"].update(mcon_label)
                                else:
                                    G.node[mcon_label]["color"] = set([nextcon3_label,])
                                    G.node[nextcon3_label]["color"] = set([mcon_label,])
#                 else:
#                     print('Not of same color: ', mcon_label, mcon.end, 'and', nextcon3_label, nextcon3.start, \
#                          abs(int(nextcon3.start) - int(mcon.end)), 'bp')
    

    return G, adj_lists


def visualizeConnectivityGraph(G, outfilepath, smallest_contig_length='"default" --> smallest considered contig size in'):
    
    graph = functools.partial(gv.Graph, format='svg')
    G1 = graph()
    
    #digraph = functools.partial(gv.Digraph, format='svg')
    #G1 = digraph()
    
    
    styles = {
        'graph': {
            'label': 'Mapped Contigs Graph, threshold_distance: ' + str(smallest_contig_length) + ' bp',
            'fontsize': '12',
            'fontcolor': 'black',
            'bgcolor': '#FFFFFF',
            'rankdir': 'BT',
            'size' : '7.75,10.25',
        },
        'nodes': {
            'fontname': 'Helvetica',
            'shape': 'rarrow',
            'fontcolor': 'black',
            'color': 'black',
            'style': 'filled',
            'fillcolor': '#e8e8e8',
        },
        'edges': {
            'style': 'solid',
            'color': 'black',
            'fontname': 'Courier',
            'fontsize': '12',
            'fontcolor': 'black',
            #'arrowhead' :'vee',
            #'dir':'both',
        }
    }
        
    def apply_styles(graph, styles):
        graph.graph_attr.update(
            ('graph' in styles and styles['graph']) or {}
        )
        graph.node_attr.update(
            ('nodes' in styles and styles['nodes']) or {}
        )
        graph.edge_attr.update(
            ('edges' in styles and styles['edges']) or {}
        )
        return graph
    
    G1 = apply_styles(G1, styles)
    
    
    def add_nodes(graph, nodes):
        for n in nodes:
            if isinstance(n, tuple):
                graph.node(n[0], **n[1])
            else:
                graph.node(n)
        return graph
    
    def add_edges(graph, edges):
        for e in edges:
            if isinstance(e[0], tuple):
                graph.edge(*e[0], **e[1])
            else:
                graph.edge(*e)
        return graph
    
    #G1 = add_nodes(G1, G.nodes())
    #G1 = add_edges(G1, G.edges())
    
    for __, n in enumerate(G.nodes()):
        
        temp = 0
        if n.startswith('+'):
            #print(n, type(n))
            temp+=1
        if '_H' in n:
            temp+=2
            
        if temp == 0: #-T
            G1.node(n, shape='rect', fillcolor='#779ECB')
        elif temp == 1: #+T
            G1.node(n, shape='larrow', fillcolor='#779ECB')
        elif temp == 2: #-H
            G1.node(n, shape='rect', fillcolor='#DEA5A4')
        elif temp == 3: #+H
            G1.node(n, shape='larrow', fillcolor='#DEA5A4')
    
    import math   
    
    for __, e in enumerate(G.edges()):
        edge_weight = G[e[0]][e[1]]['weight']
        
        if edge_weight >= 15:
            G1.edge(*e, label= ' ' + str(edge_weight), weight=str(edge_weight),
                   color='#CB6D51', penwidth = str(math.log(edge_weight, 2) + 1 ))
        else:
            G1.edge(*e, label= ' ' + str(edge_weight), weight=str(edge_weight),
                   color='#CB99C9', penwidth = str(math.log(edge_weight, 2) + 1 ))
           
    for __,n in enumerate(G.nodes()):
        if "color" in G.node[n]:
            dist_indegree = len(G.node[n]["color"])
            if dist_indegree != 1:
                G1.node(n, shape='rarrow', fillcolor='#FFFFFF')
            else:
                (con_node,) = G.node[n]["color"]
                e = (con_node, n)
                if not G.has_edge(n, con_node):
                    G.add_edge(*e, edgetype=2)
                    G1.edge( *e, color='#E2E2E2')
    
    
    
    #G1.node(n, shape='rarrow', fillcolor='#E2E2E2')
        
    #print(G1.node)
    
    #print(G1.source)
    
    filename = G1.render(filename=outfilepath)
    
    return filename
    
def buildNodeNeighborDict(G):
    diff_neighbors = {}
    
    for adj_tuple in G.adjacency_iter():
        n = adj_tuple[0]
        neighbordict = adj_tuple[1]
        
        pink_adj = []
        grey_adj = []
        
        for neighbor, edge_attrib_dict in neighbordict.items():
            if edge_attrib_dict['edgetype'] == 1:
                pink_adj.append(neighbor)
            else:
                grey_adj.append(neighbor)
            
        diff_neighbors[n] = {'pink_neighbors': pink_adj, 
                             'grey_neighbors': grey_adj}

    
    return diff_neighbors


def getCovUsedContigDict(df):
    
    cov_dict = {k +'_H ':v for k,v in df['multipl'].to_dict().items()}
    cov_dict.update({k +'_T ':v for k,v in df['multipl'].to_dict().items()})
    
    used_dict = {k:0 for k, v in cov_dict.copy().items()}
    
    return cov_dict, used_dict


def findPath(n, diff_neighbors, current_path, cov_dict, used_dict, MAX=sys.maxsize, G=False, htflag=False, guides=False):
    
    new_vale = used_dict[n] + 1
    
    if guides and new_vale > cov_dict[n]:
        print('COVERAGE EXCEEDED')
        return current_path
    else:
        used_dict[n] += 1
        current_path.append(n)
    
    if len(current_path) >= MAX or len(set(current_path)) * 3 < len(current_path) : #Protection against infinite loopS
        print('INFINITE LOOP STOP')
        return current_path 
    #print("Currently: ",current_path)
    
    
    if htflag:
        n_id, n_ht = n.split('_')
        
        if n_ht == 'H ':
            cn = '_'.join([n_id, 'T '])
        else:
            cn = '_'.join([n_id, 'H '])
    
    pn = diff_neighbors[n]['pink_neighbors'][:]
    
    if len(current_path) > 1:
        prev = current_path[-2]
        #print("Previous Node: ",prev,"Pink Neighbor Nodes: ", pn, "of current Node", n)
        #THE PROBLEM: THERE IS NO EDGE BETWEEN 28
        if prev in pn:
            pn.remove(prev)
        
    if len (current_path) > 2:
        pre_prev = current_path[-3] #print("Preprevious Node: ",pre_prev)
    else:
        pre_prev = None
    
    ### Algorithm ###
        
    if len(current_path) > 1 and n == current_path[0]:
        print("When path returns to its starting node, break and return")
        return current_path[:-1]
    
               
    if len(pn) == 0:
        print('Reached a dead end', current_path)
        return current_path
    
    elif len(pn) == 1:
        if guides and pn[0] == guides[0]:
            guides = guides[1:]
        return findPath(pn[0], diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)
    
    elif len(pn) > 1:
        if htflag:
            if cn in pn: #directly go to its opposite if connected
                if guides and cn == guides[0]:
                    guides = guides[1:]
                return findPath(cn, diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)
            
            assembly_pns = [p for p in pn if "ASSEMBLY" in G[n][p]["ref"]] #Is there a definitive vote from the assembly?
            highvote_pns = [p for p in pn if not G[n][p]["weight"] == 1]   #only those neighbors with higher vote than 1
    
            if pre_prev and diff_neighbors[pre_prev]['grey_neighbors']:
                
                gn = diff_neighbors[pre_prev]['grey_neighbors']
                print("OK", n, pn, gn, guides)
                print("Current Path:", current_path)
                if not guides and len(gn) == 1 and gn[0] in pn: # If there is exactly one grey neigbor, continue
                    return findPath(gn[0], diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G,htflag=htflag, guides=guides)
                
                if len(assembly_pns) == 1 and guides and not assembly_pns[0] == guides[0]:
                    print("Assembly graph supports an intermediate node")
                    return findPath(assembly_pns[0], diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)
                                   
                if guides and ((guides[0] in gn and guides[0] in pn) or (not gn and guides[0] in pn) or (guides[0] in pn)):
                    guidenode, guides = guides[0], guides[1:]
                    print("Guide's next element same as target of grey edge from pre-previous node, or using guide because no grey edge")
                    return findPath(guidenode, diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G,htflag=htflag, guides=guides)
            
    
            if len(assembly_pns) == 1 and guides and not assembly_pns[0] == guides[0]:
                    return findPath(assembly_pns[0], diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)
                
            if len(highvote_pns) == 1:
                    return findPath(highvote_pns[0], diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)           
                
            elif guides and guides[0] in pn:
                
#                 if len(assembly_pns) == 1 and not assembly_pns[0] == guides[0]:
#                     return findPath(assembly_pns[0], diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)
#                  
#                 if len(highvote_pns) == 1:
#                     return findPath(highvote_pns[0], diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)
#                  
#                 else:
                guidenode, guides = guides[0], guides[1:]
                print("Using guides to prolong path to", guidenode)
                return findPath(guidenode, diff_neighbors, current_path, cov_dict, used_dict, MAX=MAX, G=G, htflag=htflag, guides=guides)
                
            else:
                if guides:
                    if pre_prev:
                        gn = diff_neighbors[pre_prev]['grey_neighbors'] 
                        print('AMBIGUOUS PINK AND NO OR MULTIPLE GREY EDGES 2', n, pn, gn, cn, current_path, guides)
                    else: print('AMBIGUOUS PINK AND NO GREY EDGES', n, pn, cn, current_path, guides)
                
                print("No test postive, stop path elongation",current_path,n, cn,  pn, pre_prev, guides)            
                return current_path
#             else:
#                 
#                 print("the reason")

def findPaths(G, df, startnodes=False, htflag=False):
    """Start at every node and run through with the developed algorithm"""
    
    #print("startnodes",startnodes)
    
    diff_neighbors = buildNodeNeighborDict(G)

    cov_dict, used_dict = getCovUsedContigDict(df)
    
    sorted_nodes = sorted(G.nodes(), key=lambda x: len(diff_neighbors[x]['pink_neighbors']))
    
    MAX = len(sorted_nodes) * 2
    
    if not startnodes:
        startnodes = sorted_nodes
    
    all_paths = []
    for n in startnodes:
        current_path = []
        if G.has_node(n):
            print("Start from", n)
            all_paths.append(findPath(n, diff_neighbors, current_path, cov_dict, used_dict, G=G, MAX=MAX, htflag=htflag))
    
    for n in sorted_nodes:
        if used_dict[n] == 0:
            current_path = []
            all_paths.append(findPath(n, diff_neighbors, current_path, cov_dict, used_dict, G=G, htflag=htflag))
    
    all_paths = [p for p in all_paths if p] # filter out 'None' 
    all_paths = [p for p in all_paths if len(p) > 1]
    
    ### try to delete duplicates and subsets
    nr_paths = []
    
    for p1 in sorted(all_paths, key=lambda x: len(x), reverse=True):
        
        p1_str = ''.join(p1).strip()
        p1_str_r = ''.join(p1[::-1]).strip()
        
        current_paths_str = [''.join(p).strip() for p in nr_paths]
        
        flag = True
        for pstr in current_paths_str:
            circ_pstr = pstr + ' ' + pstr

            if not p1_str in pstr \
             and not p1_str_r in pstr \
             and not p1_str in circ_pstr \
             and not p1_str_r in circ_pstr:
                continue
            else:
                flag = False
                break
        if flag:
            nr_paths.append(p1)
    
    ### try to merge the ends of paths
    nr_paths = sorted(nr_paths, key=lambda x: len(x), reverse=True)
    
    from itertools import takewhile

    def longestPrefix(a, b):
        return [x for (x, _) in takewhile(lambda x: x[0] == x[1], zip(a, b))]

    def longestSuffix(a, b):
        return longestPrefix(reversed(a), reversed(b))[::-1]
    
    from functools import reduce

    def getAllSuffixes(a):
        return [a[i-len(a):] for i in reversed(range(len(a)))]
    
    #def getAllPrefixes(a):
    #    return [x[::-1] for x in getAllSuffixes(a[::-1])]
    
    def longestPrefixThatIsSuffix(a,b):
        lcp = []
        for suf in getAllSuffixes(b):
            if len(longestPrefix(a, suf)) == len(suf):
                lcp = suf
        return lcp
    
    def mergeSeqs(a, b, c = None, merge_mode = None):
        #print('MMM', "MODE", merge_mode, '\nA: ', a, '\nB: ', b, '\n   overlap: ', c)
        if merge_mode == 'lcp':
            return b[:-len(c)] + a
        if merge_mode == 'inv_lcp':
            return a + b[len(c):]
        if merge_mode == 'lcs':
            return a + list(reversed(b))[len(c):]
        if merge_mode == 'inv_lcs':
            return list(reversed(a))+ b[len(c):]
    
    def findAndMerge(a, b):
        d = {'lcp': longestPrefixThatIsSuffix(a, b), 
             'inv_lcp': longestPrefixThatIsSuffix(b, a), 
             'lcs': longestPrefixThatIsSuffix(list(reversed(a)), b),
             'inv_lcs': longestPrefixThatIsSuffix(b, list(reversed(a)))}
        
        m = max(d, key=lambda k: len(d[k]))

        if len(d[m]) > 2:
            a = mergeSeqs(a, b, c=d[m], merge_mode=m)
        return a
    
    def getContigID(contigend):
        for e in contigend.split('_'):
            if e.isdigit():
                return e
    
    def getInvLabel(label):
        if label.endswith('_H '):
            return label.replace('_H ', '_T ')
        else:
            return label.replace('_T ', '_H ')
    
    
    def correctScaffoldEnds(nr_paths, used_dict):
        
        for i, p in enumerate(nr_paths[:]):
                    
            if p[0] != getInvLabel(p[1]):
                nr_paths[i] = [getInvLabel(p[0])] + nr_paths[i]
                used_dict[nr_paths[i][0]] += 1
            if p[-2] != getInvLabel(p[-1]):
                nr_paths[i] += [getInvLabel(p[-1])]
                used_dict[nr_paths[i][-1]] += 1
            
            for j, c in enumerate(p):
                used_dict[c] += 1
            
            
        marked = []
        for i, p in enumerate(nr_paths):
            for j, c in enumerate(p):
                if used_dict[c] > 1:
                    if c in marked:
                        nr_paths[i][j] = False
                    else:
                        marked.append(c)
            nr_paths[i] = [x for x in nr_paths[i] if x]
        
        return nr_paths, used_dict
    
    paths = nr_paths[:]
    used_dict = {k:0 for k, v in cov_dict.copy().items()}

    nr_paths, used_dict = correctScaffoldEnds(list(sorted(nr_paths, key= lambda x: len(x), reverse=True)), used_dict)
        
    pos_to_delete = []
    
    for i, p in enumerate(nr_paths):
        a = p
        if i in pos_to_delete:
            continue
        #print('i:',i)
        for j, p1 in enumerate(nr_paths[:]):
            if j in pos_to_delete or i >= j :
                continue
            #print('j:',j)
            new_a = findAndMerge(a, p1)
            if len(new_a) > len(a):
                a = new_a
                pos_to_delete.append(j)
        paths[i] = a
    
    df1 = pd.DataFrame.from_dict(used_dict, orient='index')
    df1.columns = ['count']
    df = df.join(df1)
    
    
    return all_paths, nr_paths, sorted([p for i, p in enumerate(paths) if not i in pos_to_delete], key=lambda x: len(x), reverse=True), df


def adjlistsToFile(adj_lists, outfilepath):
    
    for k, li in adj_lists.items():
        for i, elem in enumerate(li):
            adj_lists[k][i] = (adj_lists[k][i][0], adj_lists[k][i][1], abs(adj_lists[k][i][2] - 2 * READ_LENGTH))
    
    with open(outfilepath, 'w') as of:
        for k, v in adj_lists.items():
            of.write('{}\t{}\n'.format(k, v))


def checkAndConvertHTpaths(non_redundant_paths, ):
    """First checks whether the ordering in each path
    shows consistent pairs of H and T of the same
    contig ID. Second: converts the notion of the nodes
    to a notion containing + and - for strandedness info"""
    #from itertools import tee, islice, chain, zip
    
    def grouped(s, n):
        """s transformed to (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), 
        (s2n,s2n+1,s2n+2,...s3n-1), ..."""
        return zip(*[iter(s)]*n)
    
    converted_paths = []
    
    
    for j, p in enumerate(non_redundant_paths):

        current_id = []
        current_ht = []
        orientations = []
        
        # here check whether the first contig is incomplete and if yes infer orientation
        # and let the rest of the path be iterated starting from the second element to make
        # the modulo 2 criterion work
                
        gr_path = list(grouped(p, 2))
        print(p)
        print(list(gr_path))
        print("""####\n####\n####\n####""")
        
        flag = 0
        #gr = next(gr_path)
        
        if not gr[0].split('_')[0] == gr[1].split('_')[0]:
            current_id.append(gr[0].split('_')[0])
            if  gr[0].split('_')[1].startswith('T'):
                orientations.append('+')
            else:
                orientations.append('-')
            flag = 1
            
        for i, e in enumerate(p[flag:]):
            contig_id, ht = e.split('_')
            
            if i%2 == 0:
                current_ht.append(ht)
                current_id.append(contig_id)
                if ht.startswith('T'):
                    orientations.append('-') 
                else:
                    orientations.append('+')
            else:
                if not contig_id == current_id[-1] or (contig_id == current_id[-1] and ht == current_ht[-1]):
                    print('Warning. Inconsistency in path number', j ,', at position', i)
                    current_ht.append(ht)
                    current_id.append(contig_id)
                    if ht.startswith('T'):
                        orientations.append('-') 
                    else:
                        orientations.append('+')
        
        converted_paths.append(['{}{}'.format(orientations[i], current_id[i]) for i,__ in enumerate(current_id)])
    
    return converted_paths
  

def checkAndCountOccurences(non_redundant_paths, df):
    from collections import Counter
    
    counts = Counter(non_redundant_paths[0])
    for p in non_redundant_paths[1:]:
        counts += Counter(p)
    
    count_dict = {k.split('_')[0]:int(v) for k,v in counts.items()}
    
    cov_dict = df['multipl'].to_dict()
    
    df1 = pd.DataFrame.from_dict(count_dict, orient='index')
    df1.columns = ['count']
    df = df.join(df1)

    for k,v in cov_dict.items():
        if k in count_dict:
            if not v == count_dict[k]:
                print('Warning: Coverage_criterion not fulfilled: contig', k, ' #predicted:', v, '#observed:', count_dict[k] )
        else:
            pass
    
    return df
    
def produceScaffoldRecords(final_paths, record_dict, mode='100N'):
    """ Constructs the actual scaffolded contigs and writes em in a fasta file """
    
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    used_ids = []
    
#     print('\nRECORD DICT\n')
#     for k,v in record_dict.items():
#         print(k ,':', v)
#     print('END')
    
    scaffolds_records = []
        
    for i, path in enumerate(final_paths):
        scaf_id = ['SCAFFOLD', str(i+1), 'path', '['] #contruct fasta scaffold header indicating underlying contig path
        scaf_seq = Seq("")
        
        for j, e in enumerate(path):
            scaf_id.append(e)
            used_ids.append(e[1:])
            
            if e.startswith('-'):
                contig_seq = record_dict[e[1:]].seq.reverse_complement()
            else:
                contig_seq = record_dict[e[1:]].seq
        
            if mode == '100N':
                scaf_seq += Seq('N'*100)
                scaf_seq += contig_seq
            elif mode == 'K127': ##############HERE NEEDS IMPLEMENTATION OF CHECK BREAKS, IF OVERLAP BREAK THEN DO THIS ELSE N MODE
                scaf_seq += contig_seq[127:]
            
        scaffolds_records.append( SeqRecord( scaf_seq, id = '_'.join(scaf_id) + '_]' + '_length_' + str(len(scaf_seq)), description=""))
    
    for k,v in record_dict.items():
        if not k in used_ids:
            scaffolds_records.append(v)
    
    return scaffolds_records




if __name__ == '__main__':
    
    ################################### TEST INPUT ########################################
    
#     parser, args = parse_args(['-i','/home/bardya/usr/data/Reference_Assisted_Scaffolding/Assemblies/GCF_000018445.1_ASM1844v1_genomic/assembly.spades/contigs.fasta',
#                        '-rg', '/home/bardya/usr/data/Reference_Assisted_Scaffolding/References',
#                        '-spdir', '/home/bardya/usr/data/Reference_Assisted_Scaffolding/Assemblies/GCF_000018445.1_ASM1844v1_genomic/assembly.spades',
#                        '-o', '/home/bardya/usr/data/Reference_Assisted_Scaffolding/Scaffolding',
#                         '-f'])
    
    parser, args = parse_args([])
    
    LENGTH_THRESHOLD = 601
    READ_LENGTH = 150
    PYTHON_PATH = sys.executable
    
    #spades contig name pattern
    PATTERN = re.compile('^NODE_([\d]+)_length_([\d]+)_cov_([\d*\.\d*]+).*')
    
    #######################################################################################
    
    ################################### PREPARE OUTPUT DIR ################################
    
    try:
        assert isdir(args.outdirpath), "Output directory does not exist"
    except:
        makedirs(args.outdirpath, 0o777)
    
    if not args.force_flag:
        assert len(listdir(args.outdirpath)) == 0, "Output directory not empty, try the '-f' option to force overwriting"
    
    args.outdirpath = abspath(expanduser(args.outdirpath))
    
    #######################################################################################
    
    ################################### GET & PARSE INPUT FILES ###########################
    
    #a1) contigs file
    recordlist, df = parseContigs(args.contigs_path.name, pattern=PATTERN)
    
    #a2) store and filter contigs
    records_dict, filtered_records_dict, filtered_contigs_filepath = filterAndStoreContigs(recordlist, df, pattern=PATTERN, key='lr')

    #b) reference genomes
    reference_genomes = getRefFilePaths(args.refgenomesdir)
    
    #c) (optional) spades assembly graph (needs the directory, as several files need to be parsed)
    contig_connections, unconnected = parseSpadesAssemblyGraph(args.assemblygraph, records_dict)
    
    with open(join(args.outdirpath, "gfa_contig_connections.txt"), 'w') as ccon:
        for i in sorted(contig_connections.keys(), key=lambda x: int(x)):
            ccon.write("{}\n".format(contig_connections[i]))
    
    #######################################################################################
    
    ################################### Filter Contig Connections #########################
    
    filtered_contig_connections, filtered_unconnected = filterConnections(contig_connections, unconnected, filtered_records_dict)
    
    #######################################################################################
        
    ############################ Prepare the Contig Heads and Tails #######################
    
    head_tail_filepath = extractHeadTail(filtered_contigs_filepath, args.outdirpath, 
                                         ht_max_length = 2000, margin_length=READ_LENGTH*1.5)
    
    #######################################################################################
    
    ############################## Perform Mapping using Nucmer ###########################
    
    map_outdir = join(args.outdirpath, 'mapped')
    if not isdir(map_outdir):
        makedirs(map_outdir, 0o777)
    #nucmerMapping(head_tail_filepath, args.refgenomesdir, map_outdir)
    
    #######################################################################################

    #### DETERMINE MINIMAL PHYSICAL DISTANCE OF TWO CONTIG ENDS TO BE CALLED CONNECTED ####   
          
    smallest_contig_length = min([len(rec.seq) for __,rec in filtered_records_dict.items() if len(rec.seq) >= LENGTH_THRESHOLD ])
    
    #######################################################################################
    
    ############################## Extract Contig Permutations ############################
    
#TODO: Iterative approach to also include the smaller contigs, with length smaller than alen_thresh
    perms_ht = getMaporder(join(args.outdirpath, 'mapped'), cov = 90, 
                           alen_thresh = 1300 + 2 * READ_LENGTH * 1.5,
                           #alen_thresh = smallest_contig_length + 2 * READ_LENGTH * 1.5, 
                           unique = False, htflag = True)
    
    #######################################################################################
    
    ############################## Filter Perms according to Records Dict #################
    
    filtered_perms_ht = filterPerms(perms_ht, filtered_records_dict)
    
    #######################################################################################
    
    ############################## Build the Graph ########################################
    
    G = initGraph(filtered_contig_connections)

    G, adj_lists = toGraph(filtered_perms_ht, G=G, 
                           dist_thresh = smallest_contig_length + 2 * READ_LENGTH * 1.5, 
                           htflag = True)
    
    #######################################################################################    
    
    #### SAVE TO FILE THE DISTANCES OF MAPPED HEAD AND TAIL REGIONS IN ALL REFERENCES ##### 
     
    adjlistsToFile(adj_lists, join(args.outdirpath, 'ht_mapping_distances.tsv'))
    
    #######################################################################################
    
    ############################## VISUALIZE THE GRAPH ####################################
    
    visualizeConnectivityGraph(G, join(args.outdirpath, 'ht_connectivity_graph'), 
                               smallest_contig_length=smallest_contig_length)
    
    #######################################################################################
    
    ################################### FIND PATHS ########################################
    filtered_df = df[df["length"] >= smallest_contig_length]

    all_paths, overlap_paths, non_redundant_paths, filtered_df = findPaths(G, filtered_df, startnodes=filtered_unconnected, 
                                                                  htflag=True)
        
    for p in non_redundant_paths:
        print(p)
    
    #######################################################################################
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print("Starting Refinement...")
    
    #######################################################################################
    
    ## DETERMINE NEW MINIMAL PHYSICAL DISTANCE OF TWO CONTIG ENDS TO BE CALLED CONNECTED ##
          
    smallest_contig_length = min([len(rec.seq) for __,rec in records_dict.items() if len(rec.seq) >= LENGTH_THRESHOLD ])
    
    #######################################################################################
    
    ############################## Build 2nd Graph ########################################
    
    G2 = initGraph(contig_connections)

    G2, adj_lists = toGraph(perms_ht, G=G2, 
                           dist_thresh = smallest_contig_length + 2 * READ_LENGTH * 1.5, 
                           htflag = True)
    
    #######################################################################################  
    
    ############################## VISUALIZE 2nd GRAPH ####################################
    
    visualizeConnectivityGraph(G2, join(args.outdirpath, 'ht_connectivity_graph2'), 
                               smallest_contig_length=smallest_contig_length)
    
    
    ################################### FIND 2nd PATHS ########################################
    extended_paths = []
    diff_neighbors = buildNodeNeighborDict(G2)
    cov_dict, used_dict = getCovUsedContigDict(df)
    
    for guide_path in sorted(non_redundant_paths, key=lambda x: len(x), reverse=True):
        
        startnode = guide_path[1]
        current_path = [guide_path[0]]
        if G2.has_node(startnode):
            used_dict[guide_path[0]] += 1
            #print(guide_path)
            extended_paths.append(findPath(startnode, diff_neighbors, current_path, cov_dict, used_dict, MAX=len(guide_path*3), G=G2, htflag=True, guides=guide_path[2:]))
    
    for n in sorted(G2.nodes()):
        if used_dict[n] == 0:
            current_path = []
            extended_paths.append(findPath(n, diff_neighbors, current_path, cov_dict, used_dict, G=G2, htflag=True))
            
    for p in extended_paths:
        print(p)
    
    #######################################################################################

    
    ##### Finalize and produce output files ####################
    
    df = checkAndCountOccurences(extended_paths, df)
        
    df.to_csv(join(args.outdirpath, 'contigs_stats.tsv'), sep= '\t')
    
    final_paths = checkAndConvertHTpaths(extended_paths)
    
    print("Final_paths:")
    for p in final_paths:
        print(p)
    
    scaffold_records = produceScaffoldRecords(final_paths, records_dict)
    
    outfile = open(join(args.outdirpath, 'scaffolds.fasta'), 'w')
    
    for record in scaffold_records:
        SeqIO.write(record, outfile, "fasta")
    outfile.close()
    
    
    ###################### PRINT CONTIG FOLD COVERAGES & STATS ############################
    
    df.to_csv(join(args.outdirpath, 'contigs_stats.tsv'), sep= '\t')
    