'''
Created on Jul 25, 2016

@author: bardya
'''
import argparse
import glob
import os
import subprocess
from Bio import SeqIO
import pandas as pd

PYTHON_PATH = '/home/bardya/anaconda3/envs/py3k/bin/python'
SCRIPTS_PATH = os.path.join(os.path.dirname(__file__), 'scripts')

def parse_args():
    parser = argparse.ArgumentParser(description="Multi-Reference-Genome-Assisted Scaffolding of target assembly")
    
    parser.add_argument('-i', '--contigs_filepath', dest='contigs_path', metavar='</path/to/coordsfile>', type=argparse.FileType('rt'),
                        help='path to an contigs fasta file')
    
    parser.add_argument('-o', '--outdir', dest='outdirpath', metavar='</path/to/out/>', type=str,
                        default=os.getcwd(), 
                        help='path to output file prefix'),
    
    parser.add_argument('-rg', '--reference-genomes-directory', dest='refgenomesdir', metavar='<path>', type=str, 
                        default='./references', 
                        help='specify a directory containing the reference genomes (multi-) fasta nucleotide sequence files'),
    
    parser.add_argument('-rm_full', '--reference_mappings_dir_full', dest='refmapdir', metavar='<path>', type=str, 
                        default='./assembly.ref_mapping', 
                        help='specify a directory containing the nucmer mappings of the contigs against each reference genome'),
    
    parser.add_argument('-rm_ht', '--reference_mappings_dir_ht', dest='ht_refmapdir', metavar='<path>', type=str, 
                        default= None, 
                        help='specify a directory containing the nucmer mappings of the contigs\' ht against each reference genome'),
                        
    parser.add_argument('-in_ht', '--contigs_ht_filepath', dest='contigs_path_ht', metavar='<path>', type=str, 
                        default= None, 
                        help='specify the path to a contigs_ht fasta file'),
    
    parser.add_argument('-mlen', '--margin_length', dest='margin_length', metavar='<int>', type=int,
                        default=150,
                        help='for head/tail mapping, define a margin size at both ends of the contig'),
                        
    parser.add_argument('-htlen', '--ht_length', dest='ht_length', metavar='<int>', type=int,
                        default=2000,
                        help='for head/tail mapping, define a margin size at both ends of the contig'),
    
    parser.add_argument('-cov_thresh', '--coverage_threshold', dest='cov_thresh', metavar='<int>', type=int,
                        default=40,
                        help='min % alignment coverage indicating a successfull mapping of a head or tail sequence'),
    
    parser.add_argument('--skip_ht', dest='skip_ht', action='store_true',
                        help='if this flag is set, the script will not perform a scaffolding based on head and tail sequences'),
                        
    parser.add_argument('--skip_full', dest='skip_full', action='store_true',
                        help='if this flag is set, the script will not perform a scaffolding with the full contigs'),
                        
    parser.add_argument('--manual_scaffolds', dest='manual_scaffolds_filepath', metavar='<path>', type=str, 
                        default= None, 
                        help='if given, the program will take the contig file and produce a scaffold file based on the information in this file.')
                    
    parser.add_argument('--version', action='version', version='0.14')
    
    return parser, parser.parse_args()


def getMaporder(contigs_path, refmapdir, cov = 70, alen_thresh = 1400, ht_len = 2000, unique = True, htflag = False):
    '''get the mapping order of the contigs which meet the 
    thresholds for each reference in the reference directory'''
        
    file_ending_glob = '*_dnadiff.mcoords'
    
    import coords2gff3
    
    args2 = argparse.Namespace()
    args2.__dict__.update({'coverage_threshold': cov,
                           'length_threshold' : alen_thresh,
                           'unique_requirement' : unique,
                           'source' : 'contig',
                           'htflag' : htflag})
    
    perms = []

    for ref_map_path in glob.glob(os.path.join(refmapdir, file_ending_glob)):
        args2.__dict__['infilepath'] = open(ref_map_path ,'r') 
        feats = coords2gff3.coordsFileParser(args2.__dict__['infilepath'], args2)
        
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


def extractHeadTail(contigs_path, outdirpath, ht_length = 2000, margin_length= 150):
    
    cmd = [PYTHON_PATH, os.path.join(SCRIPTS_PATH, 'fasta_headtail_extract.py'), 
           '-i', contigs_path,
            '-min', str(ht_length+2*margin_length), 
            '-l', str(ht_length), 
            '-ms', str(margin_length), 
            '-o', os.path.join(outdirpath, os.path.basename(contigs_path).split('.')[0])]
    
    exit_code = subprocess.call(' '.join(cmd), shell=True)
    
    if exit_code > 0 :
        print('Error while trying to extract the head and tail sequences from fasta file:', contigs_path)
        raise
    else:
        outpath = os.path.join(outdirpath, os.path.basename(contigs_path).split('.')[0] + '_HT.fasta')
    
    return outpath
    

def nucmerMapping(contigs_path, refgenomesdir, outdir):
    cmd = [PYTHON_PATH, os.path.join(SCRIPTS_PATH, 'nucmer_mapping.py'), 
           '-i', contigs_path, '-o', outdir,
           '-rd', refgenomesdir]

    #print(' '.join(cmd))
    
    exit_code = subprocess.call(' '.join(cmd), shell=True)
    
    if exit_code > 0 :
        print('Error while trying to map', contigs_path,'to reference genome sequences in', refgenomesdir)
        raise IOError()
    else:
        return


def toGraph(perms, dist_thresh=1300, htflag = False):
    import networkx as nx
    
    #G = nx.DiGraph()
    G = nx.Graph()
        
    from itertools import cycle
    
    def getNodeLabel(feat, htflag):
        """Gives out a 2-tuple, with the long label with orientation, 
        id and and head/tail info on the contig, and a shorter version
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
            #print(k)
            adj_lists[k] = []
            
            all_connected = cycle(perm)
            next(all_connected) #to start with element two
            for mcon in perm:
                inverse = False
                mcon_label, lab = getNodeLabel(mcon,  htflag= htflag)
                
                if not htflag and G.has_node(getInvLabel(mcon_label)):
                        mcon_label = lab = getInvLabel(mcon_label)
                        inverse = True
                
                G.add_node(mcon_label)
               
                nextcon = next(all_connected)
                nextcon_label, nlab = getNodeLabel(nextcon, htflag= htflag )
                        
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
                else:
                    print('No edge drawn. Too far apart: ', mcon_label, mcon.end, 'and', nextcon_label, nextcon.start, \
                         abs(int(nextcon.start) - int(mcon.end)), 'bp')
            
    #second iteration, now going for contigs spanning two
    for permdic in perms:
        for k, perm in permdic.items():
            #print(k)
            
            all_connected = cycle(perm)
            all_connected2 = cycle(perm)
            all_connected3 = cycle(perm)
            
            [next(all_connected) for i in range(1)] #to later start with element two
            [next(all_connected2) for i in range(2)] #to later start with element three
            [next(all_connected3) for i in range(3)] #to later start with element four
            
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
                else:
                    print('Not of same color: ', mcon_label, mcon.end, 'and', nextcon3_label, nextcon3.start, \
                         abs(int(nextcon3.start) - int(mcon.end)), 'bp')
    
        
    return G, adj_lists
    

def visualizeConnectivityGraph(G, outfilepath, smallest_contig_length='"default" --> smallest considered contig size in'):
    
    import graphviz as gv
    import functools
    
    graph = functools.partial(gv.Graph, format='svg')
    
    digraph = functools.partial(gv.Digraph, format='svg')
    
    G1 = graph()
    
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
    
    
    for i, n in enumerate(G.nodes()):
        
        temp = 0
        if n.startswith('+'):
            #print(n, type(n))
            temp+=1
        if '_H' in n:
            temp+=2
            
        if temp == 0: #-T
            G1.node(n, shape='rarrow', fillcolor='#779ECB')
        elif temp == 1: #+T
            G1.node(n, shape='larrow', fillcolor='#779ECB')
        elif temp == 2: #-H
            G1.node(n, shape='rarrow', fillcolor='#DEA5A4')
        elif temp == 3: #+H
            G1.node(n, shape='larrow', fillcolor='#DEA5A4')
            
    import math   
    for i, e in enumerate(G.edges()):
        edge_weight = G[e[0]][e[1]]['weight']
        
        if edge_weight >= 15:
            G1.edge(*e, label= ' ' + str(edge_weight), weight=str(edge_weight),
                   color='#CB6D51', penwidth = str(math.log(edge_weight, 2) + 1 ))
        else:
            G1.edge(*e, label= ' ' + str(edge_weight), weight=str(edge_weight),
                   color='#CB99C9', penwidth = str(math.log(edge_weight, 2) + 1 ))
    
            
    for i,n in enumerate(G.nodes()):
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
    
    filename = G1.render(filename=os.path.join(outfilepath))


def findPaths(G, htflag=False):
    """Start at every node and run through with the developed algorithm"""
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
    
    sorted_nodes = sorted(G.nodes(), key=lambda x: len(diff_neighbors[x]['pink_neighbors']))
    
    def findPath(n, diff_neighbors, current_path):
        current_path.append(n)
        if len(current_path) >= MAX or len(set(current_path))*3 < len(current_path) : #Protection against infinite loopS
            return current_path 
        #print("Currently: ",current_path)
        
        if htflag:
            n_id, n_ht = n.split('_')
            
            if n_ht == 'H ':
                cn = '_'.join([n_id, 'T '])
            else:
                cn = '_'.join([n_id, 'H '])
        
        pn = diff_neighbors[n]['pink_neighbors'][:]
        
        #print(n,pn)
        
        pre_prev = None
        
        if len(current_path) > 1:
            prev = current_path[-2]
            #print("Previous Node: ",prev,"Pink Neighbor Nodes: ", pn, "of current Node", n)
            #THE PROBLEM: THERE IS NO EDGE BETWEEN 28
            
            pn.remove(prev)
        if len (current_path) > 2:
            pre_prev = current_path[-3]
            #print("Preprevious Node: ",pre_prev)
                
        if len(pn) == 0:
            #print('DONE')
            return current_path
        
        if len(current_path) > 1 and n == current_path[0]:
            #print("When path returns to its starting node, break and return")
            return current_path[:-1]
        
        elif len(pn) == 1:
            return findPath(pn[0], diff_neighbors, current_path)
    
        elif len(pn) > 1:
            if htflag:
                if cn in pn:
                    return findPath(cn, diff_neighbors, current_path)
                elif pre_prev:
                    gn = diff_neighbors[pre_prev]['grey_neighbors']
                    if len(gn) == 1:
                        if gn[0] in pn:
                            return findPath(gn[0], diff_neighbors, current_path)
                        else:
                            return current_path
                    else:
                        return current_path
            else:
                if pre_prev:
                    gn = diff_neighbors[pre_prev]['grey_neighbors']
                    if len(gn) == 1:
                        return findPath(gn[0], diff_neighbors, current_path)
                    else:
                        return current_path
                
        
    all_paths = []
    
    MAX = len(sorted_nodes) * 2
    
    for n in sorted_nodes:
        current_path = []
        prev = ''
        all_paths.append(findPath(n, diff_neighbors, current_path))
    
    all_paths = [p for p in all_paths if p] # filter out 'None' 
    all_paths = [p for p in all_paths if len(p) > 1]
    
    #print(all_paths)

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
    
    #print(nr_paths)
    
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
    
    paths = nr_paths[:]
    
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
    
    return all_paths, nr_paths, sorted([p for i, p in enumerate(paths) if not i in pos_to_delete], key=lambda x: len(x), reverse=True)
    #h = nx.adjacency_matrix(G, nodelist=None, weight='weight')
    #print(h)


def mappingOrderToFile(perms, outfilepath):
    with open(outfilepath, 'w') as of:
        for permdic in perms:            
            for k,v in permdic.items():
                of.write('{}\t{}\n'.format(k, v))


def adjlistsToFile(adj_lists, outfilepath):
    with open(outfilepath, 'w') as of:
        for k, v in adj_lists.items():
            of.write('{}\t{}\n'.format(k, v))


def getCoverages(contigs_path, per_base_coverages = False):
    
    if per_base_coverages:
        pass
    else:
        import re
        
        #spades contig name pattern
        pattern = re.compile('^NODE_([\d]+)_length_([\d]+)_cov_([\d*\.\d*]+).*')
        # Read the file
        input_seq_it = SeqIO.parse(open(contigs_path, "r"), "fasta")
        
        seq_vallist = []
        index_list = []
        recordlist = []
        
        for record in input_seq_it:
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
        
        record_dict = {}
        for seqrec in recordlist:
            m = re.findall(pattern, seqrec.id)
            if m:
                record_dict[seqrec.id.split('_')[1]] = seqrec
            else:
                k = re.findall(r'\d+', seqrec.id)
                record_dict[str(int(k[0]))] = seqrec 
        
        #print(record_dict.items())
        
        return df, record_dict


def checkAndConvertHTpaths(non_redundant_paths, ):
    """First checks whether the ordering in each path
    shows consistent pairs of H and T of the same
    contig ID. Second: converts the notion of the nodes
    to a notion containing + and - for strandedness info"""
    #from itertools import tee, islice, chain, zip
    
    def grouped(iterable, n):
        """s transformed to (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), 
        (s2n,s2n+1,s2n+2,...s3n-1), ..."""
        return zip(*[iter(iterable)]*n)
    
    converted_paths = []
    
    
    for j, p in enumerate(non_redundant_paths):

        current_id = []
        current_ht = []
        orientations = []

        # here check whether the first element is incomplete and if yes infer orientation
        # and let the rest of the path be iterated starting from the second element to make
        # the modulo 2 criterion work
                
        gr_path = grouped(p, 2)
        flag = 0
        gr = next(gr_path)
        
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
    
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    used_ids = []
    
#     print('\nRECORD DICT\n')
#     for k,v in record_dict.items():
#         print(k ,':', v)
#     print('END')
    
    scaffolds_records = []
        
    for i, path in enumerate(final_paths):
        scaf_id = ['SCAFFOLD', str(i+1), 'path', '[']
        scaf_seq = Seq("")
        
        for j, e in enumerate(path):
            scaf_id.append(e)
            used_ids.append(e[1:])
            
            if e.startswith('-'):
                scaf_seq += record_dict[e[1:]].seq.reverse_complement()
            else:
                scaf_seq += record_dict[e[1:]].seq
        
            if not j == len(path)-1:
                if mode == '100N':
                    scaf_seq += Seq('N'*100)
            
        scaffolds_records.append( SeqRecord( scaf_seq, id = '_'.join(scaf_id) + '_]' + '_length_' + str(len(scaf_seq)), description=""))
    
    for k,v in record_dict.items():
        if not k in used_ids:
            scaffolds_records.append(v)
    
    return scaffolds_records


if __name__ == '__main__':
    
    parser, args = parse_args()
#     args = parser.parse_args(['-i','/share/project/bardya/Acinetobacter/tmp/simulated_assemblies/GCF_000018445.1_ASM1844v1_genomic/assembly.spades/contigs.fasta',
#                        '-o','/share/project/bardya/Acinetobacter/tmp/simulated_assemblies/GCF_000018445.1_ASM1844v1_genomic/scaffolding_reffi',
#                        '-rm_full','/share/project/bardya/Acinetobacter/tmp/simulated_assemblies/GCF_000018445.1_ASM1844v1_genomic/assembly.ref_mapping',
#                        '-rg', '/share/project/bardya/Acinetobacter/tmp/simulated_assemblies/references',
#                        '-rm_ht', '/share/project/bardya/Acinetobacter/tmp/simulated_assemblies/GCF_000018445.1_ASM1844v1_genomic/scaffolding_reffi/ht_length_2000_margin_length_0',
#                        '-in_ht', '/share/project/bardya/Acinetobacter/tmp/simulated_assemblies/GCF_000018445.1_ASM1844v1_genomic/scaffolding_reffi/ht_length_2000_margin_length_0/contigs_HT.fasta',
#                        '--skip_full'])
    
    ###############################Using the full length of the contig###################################
    if not args.skip_full and not args.manual_scaffolds_filepath:
        if not args.refmapdir:
            args.refmapdir =  os.path.join(args.outdirpath, 'full_ref_mapping')
            
            if not os.path.exists(args.refmapdir):
                os.makedirs(args.refmapdir)
            
            nucmerMapping(args.contigs_path.name, args.refgenomesdir, args.refmapdir)
            
        perms = getMaporder(args.contigs_path.name, args.refmapdir, unique = False)
        
        mappingOrderToFile(perms, os.path.join(args.outdirpath, 'full_map_order.tsv'))
        
        G, adj_lists = toGraph(perms, dist_thresh = 1000)
        
        adjlistsToFile(adj_lists, os.path.join(args.outdirpath, 'full_adj_mapping_distances.tsv'))

        visualizeConnectivityGraph(G, os.path.join(args.outdirpath, 'full_connectivity_graph'))
        
        non_redundant_paths = findPaths(G, htflag=False)
        
        #print(non_redundant_paths)
    
    ###########################Using the head and tail of the contig#####################################
    
    if not args.skip_ht and not args.manual_scaffolds_filepath:
        if not args.ht_refmapdir:
            dirname = 'ht_length_' + str(args.ht_length) + '_margin_length_' + str(args.margin_length)
            args.ht_refmapdir =  os.path.join(args.outdirpath, dirname)
            
            if not os.path.exists(args.ht_refmapdir):
                os.makedirs(args.ht_refmapdir)
            
            if not args.contigs_path_ht:
                args.contigs_path_ht = extractHeadTail(args.contigs_path.name, args.ht_refmapdir, ht_length = args.ht_length, margin_length = args.margin_length)
    
            nucmerMapping(args.contigs_path_ht, args.refgenomesdir, args.ht_refmapdir)
        
        #GET THE ORDER OF MAPPED CONTIGS
        perms_ht = getMaporder(args.contigs_path_ht, args.ht_refmapdir, cov = args.cov_thresh, alen_thresh = args.ht_length * (args.cov_thresh/100), ht_len = 2000, unique = False, htflag = True)
                
        mappingOrderToFile(perms_ht, os.path.join(args.outdirpath, 'ht_map_order.tsv'))
        
        df, record_dict = getCoverages(args.contigs_path.name)
        
        smallest_contig_length = min([len(rec.seq) for __,rec in record_dict.items() if len(rec.seq) > 2*(args.ht_length + args.margin_length) ])
        
        smallest_contig_length = 1400
        
        print("Smallest Contig Size / Max gap between mapped contigs to create a link : ", smallest_contig_length)
        
        G, adj_lists = toGraph(perms_ht, dist_thresh = smallest_contig_length + 2 * args.margin_length, htflag = True)
        
        for k, li in adj_lists.items():
            for i, elem in enumerate(li):
                adj_lists[k][i] = (adj_lists[k][i][0], adj_lists[k][i][1], abs(adj_lists[k][i][2] - 2 * args.margin_length))
         
        adjlistsToFile(adj_lists, os.path.join(args.outdirpath, 'ht_adj_mapping_distances.tsv'))
        
        visualizeConnectivityGraph(G, os.path.join(args.outdirpath, 'ht_connectivity_graph'), smallest_contig_length=smallest_contig_length)
        
        all_paths, overlap_paths, non_redundant_paths = findPaths(G, htflag=True)
        
        with open(os.path.join(args.outdirpath, 'all_paths.txt'), 'w') as all_paths_out:
            for path in sorted(all_paths, key=lambda x: len(x)):
                all_paths_out.write("[" + ",".join(path) + "]\n")
        
        if not non_redundant_paths:
            print("No paths found. No scaffolding performed as output file looks like input file. \nSuggest checking thresholds for mapping of contigs /-HT.")
            exit()
        
        df = checkAndCountOccurences(non_redundant_paths, df)
        
        df.to_csv(os.path.join(args.outdirpath, 'contigs_stats.tsv'), sep= '\t')
        
        final_paths = checkAndConvertHTpaths(non_redundant_paths)
        
        print("Final_paths:")
        for p in final_paths:
            print(p)
        
        scaffold_records = produceScaffoldRecords(final_paths, record_dict)
        
        outfile = open(os.path.join(args.outdirpath, 'scaffolds.fasta'), 'w')
        
        for record in scaffold_records:
            SeqIO.write(record, outfile, "fasta")
        outfile.close()
    
    if args.manual_scaffolds_filepath:
        df, record_dict = getCoverages(args.contigs_path.name)
        
        with open(args.manual_scaffolds_filepath, 'r') as f:
            final_paths = f.read().splitlines() 
        
        try:
            import ast
            final_paths = [ast.literal_eval(x.strip()) for x in final_paths]
        except:
            final_paths = [x.strip().split(',') for x in final_paths]
        
        df.to_csv(os.path.join(args.outdirpath, 'contigs_stats.tsv'), sep= '\t')
        scaffold_records = produceScaffoldRecords(final_paths, record_dict)
        
        outfile = open(os.path.join(args.outdirpath, 'scaffolds.fasta'), 'w')
        
        scaffold_records = list(scaffold_records)
        sorted_records = sorted(scaffold_records, key=lambda x: len(x), reverse=True)
        
        for record in sorted_records:
            SeqIO.write(record, outfile, "fasta")
        outfile.close()