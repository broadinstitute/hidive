from collections import defaultdict
import sys
sys.path.append("../AnchorbasedGraphicalGenome/")
import CCGG_extension as CCGG
import numpy
import json
import argparse
import AGG
from Levenshtein import distance

################input#################
# ref = 'MHC-CHM13'
# source_seq = graph.nodes[sanchor]['seq'] ### input
# sink_seq = graph.nodes[eanchor]['seq']  ### input
# source_pos = int(graph.nodes[sanchor][ref])-1
# sink_pos = int(graph.nodes[eanchor][ref])-1
# print(source_pos,sink_pos )
# ref_HLAC = contig[source_pos:sink_pos] ### input


# filename = './HLA-A/HLA-A_anchorgraph.aln.fa'





#################Functions#####################

def loadFasta(filename):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into a list of headers 
        and fragment sequences for each sequence contained.
        The resulting sequences are 0-indexed! """
    if (filename.endswith(".gz")):
        fp = gzip.open(filename, 'rb')
    else:
        fp = open(filename, 'rb')
    # split at headers
    data = fp.read().decode().split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    headers = []
    sequences = []
    for sequence in data:
        lines = sequence.split('\n')
        headers.append(lines.pop(0))
        sequences.append(''.join(lines))
    return (headers, sequences)

def find_reference_and_source_sink_anchors(reference, genestart, geneend):
    header, sequences = loadFasta(reference)
    chromosome = 6
    chr6 = sequences[chromosome-1].upper()
    ref_HLA = chr6[genestart:geneend]
    return ref_HLA
    
def get_reference_kmer_profile(ref_HLA, k):
    ref_kmer_profile = {}
    for i in range(len(ref_HLA) - k + 1):
        kmer = ref_HLA[i:i+k]
        ref_kmer_profile[kmer] = ref_kmer_profile.get(kmer,0) + 1

    unique_kmerlist = []
    for kmer, count in ref_kmer_profile.items():
        if count == 1:
            unique_kmerlist.append(kmer)
    
    return unique_kmerlist

def classify_sequences_by_samples(headers, sequences):
    Sample_reads = {}
    for i, h in enumerate(headers):
        seq = sequences[i]
        sample = h.split("|")[-1]
        Sample_reads[sample] = Sample_reads.get(sample, []) +[h]
    return Sample_reads

def reverse_complement(kmer):
    return ''.join([{'A':'T','C':'G','G':'C','T':'A', 'N':'N'}[base] for base in reversed(kmer)])

def find_sequences_between_sanchor_eanchor(headers, sequences,ref_HLA):
    count = 0
    HLA_samples = set()

    HLA_seq = {}

    for i, seq in enumerate(sequences):
        h = headers[i]
        seq = seq.upper()
        HLA_samples.add(h.split("|")[-1])
        HLA_seq[h] = seq
   
    HLA_seq['MHC-CHM13'] = ref_HLA
    HLA_samples.add('MHC-CHM13')

    return count, HLA_samples, HLA_seq

def map_reference_unique_kmers_to_seq(unique_kmerlist, HLA_seq, k):
    SampleDict = defaultdict(list)
    PositionDict = defaultdict(dict)
    for kmer in unique_kmerlist:
        kmer_rev = reverse_complement(kmer)
        SampleDict[kmer]
        SampleDict[kmer_rev]
        PositionDict[kmer]
        PositionDict[kmer_rev]

    for read, seq in HLA_seq.items():
        for i in range(1, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer in SampleDict:
                SampleDict[kmer] = SampleDict.get(kmer, []) + [read.split("|")[-1]]
                PositionDict[kmer] = PositionDict.get(kmer, {})
                PositionDict[kmer][read] = PositionDict[kmer].get(read, []) + [i]
    return SampleDict, PositionDict

def get_anchor_information(SampleDict,HLA_samples, PositionDict, k):
    anchorlist = []
    for kmer, samplelist in SampleDict.items():
        if len(set(samplelist)) == len(HLA_samples):
            D = PositionDict[kmer]
            if len(D)>1 and len(D["MHC-CHM13"]) == 1:
                anchorlist.append(kmer)

    AnchorInfo = {}
    for kmer in anchorlist:
        pos = PositionDict[kmer]['MHC-CHM13'][0]
        anchorname = "A%06d" %  (pos//k + 1)
        AnchorInfo[anchorname] = {}
        AnchorInfo[anchorname]['seq'] = kmer
        AnchorInfo[anchorname]['pos'] = pos

    anchornames = sorted(AnchorInfo.keys())
    anchor_unadjacent_list = []
    index = 0
    sanchor = anchornames[index]
    while sanchor < anchornames[-1]:
        for danchor in anchornames[index+1:]:
            if AnchorInfo[danchor]['pos'] - AnchorInfo[sanchor]['pos'] > k+1:
                break
            index += 1
        anchor_unadjacent_list += [sanchor, danchor]
        sanchor = danchor
    anchor_unadjacent_list = sorted(set(anchor_unadjacent_list))[:-1]
    Final_anchor = {anchor:AnchorInfo[anchor] for anchor in anchor_unadjacent_list}

    return Final_anchor

def mapping_info(AnchorInfo, contig, k):
    anchor_list = list(AnchorInfo.keys())
    Anchorseq = {AnchorInfo[anchor]['seq']:anchor for anchor in anchor_list}

    seqlist = Anchorseq.keys()
    PositionDict = defaultdict(list)
    for anchor_seq in seqlist:
        anchor_rev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(anchor_seq)])
        PositionDict[anchor_seq]
        PositionDict[anchor_rev]

    for i in range(1, len(contig) - k + 1):
        kmer = contig[i:i+k]
        if kmer in PositionDict:
            PositionDict[kmer] = PositionDict.get(kmer, []) + [i]
            
    A = {}
    SVs = {}
    for anchor, D in AnchorInfo.items():
        anchor_seq = D['seq']
        #anchor_rev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(anchor_seq)])
        poslist = PositionDict[anchor_seq] #  + PositionDict[anchor_rev]
        if len(poslist) == 1:
            A[anchor] = poslist[0]
        else:
            SVs[anchor] = poslist
            
    return A, SVs

def find_edge_info(src_pos, dst_pos, k, contig, contigname, sample, Anchorseq):
    E = {} # edgeinfo
    # get source infomation
    if src_pos == 0:
        src = "SOURCE"
        src_seq = ""
        pr = False
    else:
        src_seq = contig[src_pos:src_pos + k]
        try:
            src = Anchorseq[src_seq]
            pr = False
        except:
            src_seq = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(src_seq)])
            src = Anchorseq[src_seq]
            pr = True

    dst_seq = contig[dst_pos:dst_pos+k]
    
    if dst_pos == len(contig): # fix sink edge issue 08/29/23
        dst = "SINK"
        dst_seq = ""
        sr = True
    else:
        try:
            dst = Anchorseq[dst_seq]
            sr = False
        except:
            dst_seq = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(dst_seq)])
            dst = Anchorseq[dst_seq]
            sr = True
            
    if src_pos == 0:
        edge_seq = contig[src_pos:dst_pos] # first edge fix bug 08/28/2023
    else:
        edge_seq = contig[src_pos+k:dst_pos] # fix bug 08/23/2023
        
    if pr and sr:
        edge_seq = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(edge_seq)])
        node = src
        src = dst
        dst = node 

    
    E = {}
    E['seq'] = edge_seq
    E['src'] = src
    E['dst'] = dst
    E['reads'] = [contigname]
    E['strain'] = [sample]


    return E


# loop through all reads
def create_edgefile(HLA_seq, AnchorInfo, k):
    Edge_info = {}
    Outgoing = {}
    
    edgenum_perread = []
    anchor_list = list(AnchorInfo.keys())
    Anchorseq = {AnchorInfo[anchor]['seq']:anchor for anchor in anchor_list}
    contig_index = 0
    
    for contig_name, contig in HLA_seq.items():
        contig_index += 1
        sample_name = contig_name.split('|')[-1]
        A, SVs = mapping_info(AnchorInfo, contig, k)
        splitposlist = sorted(A.values())
        edgeindex = 0 # 
        src_pos = 0
        for dst_pos in splitposlist:
            E = find_edge_info(src_pos, dst_pos, k, contig, contig_name, sample_name, Anchorseq)
            src = E['src']
            edgelist = Outgoing.get(src, [])
            for edge in edgelist:
                if Edge_info[edge]['dst'] != E['dst']:
                    continue
                if Edge_info[edge]['seq'] == E['seq']:
                    Edge_info[edge]['reads'] += E['reads']
                    Edge_info[edge]['strain'] = sorted(set(Edge_info[edge]['strain']) | set(E['strain']))
                    break
            else:
                edgename = "E%05d.%04d" % (contig_index, edgeindex)
                Edge_info[edgename] = E
                Outgoing[src] = Outgoing.get(src,[]) + [edgename]
                edgeindex += 1
            # update
            src_pos = dst_pos
        
        dst_pos = len(contig) # fix bug 08/29/2023 
        E = find_edge_info(src_pos, dst_pos, k, contig, contig_name, sample_name, Anchorseq)
        src = E['src']
        edgelist = Outgoing.get(src, [])
        for edge in edgelist:
            if Edge_info[edge]['dst'] != E['dst']:
                continue
            if Edge_info[edge]['seq'] == E['seq']:
                Edge_info[edge]['reads'] += E['reads']
                Edge_info[edge]['strain'] = sorted(set(Edge_info[edge]['strain']) | set(E['strain']))
                break
        else:
            edgename = "E%05d.%04d" % (contig_index, edgeindex)
            Edge_info[edgename] = E
            Outgoing[src] = Outgoing.get(src,[]) + [edgename]
            edgeindex += 1
        
        edgenum_perread.append(edgeindex)   
        
    return Edge_info, Outgoing


def write_gfa(AnchorInfo, Edge_info, outputfilename):
    header = ['H\tVN:Z:1.0\n']
    # add json annotation to edge and anchor sequences
    AnchorS = []
    for anchor,D in AnchorInfo.items():
        seq = D.pop('seq')
        json_string = json.dumps(D)
        AnchorS += ['S\t%s\t%s\t%s\n' % (anchor, seq, "PG:J:" + json_string)]
   
    EdgeS = []   
    Link = [] 
    for edge,edge_dict in Edge_info.items():
        seq = edge_dict.pop('seq')
        src = edge_dict.pop('src')
        dst = edge_dict.pop('dst')
        
        json_string = json.dumps(edge_dict)
        if "reads" in edge_dict:
            EdgeS += ['S\t%s\t%s\t%s\t%s\n' % (edge, seq, "PG:J:" + json_string, "RC:i:" + str(len(edge_dict['reads'])))]
        else:
            EdgeS += ['S\t%s\t%s\n' % (edge, seq)]
        Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (src, "+", edge, "+", "0M"))
        Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (edge, "+", dst, "+", "0M"))     

    with open(outputfilename, 'w') as fp:
        for h in header:
            fp.write(h)
        for s in AnchorS:
            fp.write(s)
        for s in EdgeS:
            fp.write(s)
        for l in Link:
            fp.write(l)

def write_gfa(AnchorInfo, Edge_info, outputfilename):
    header = ['H\tVN:Z:1.0\n']
    # add json annotation to edge and anchor sequences
    AnchorS = []
    for anchor,D in AnchorInfo.items():
        seq = D.pop('seq')
        json_string = json.dumps(D)
        AnchorS += ['S\t%s\t%s\t%s\n' % (anchor, seq, "PG:J:" + json_string)]
   
    EdgeS = []   
    Link = [] 
    for edge,edge_dict in Edge_info.items():
        seq = edge_dict.pop('seq')
        src = edge_dict.pop('src')
        dst = edge_dict.pop('dst')
        
        json_string = json.dumps(edge_dict)
        if "reads" in edge_dict:
            EdgeS += ['S\t%s\t%s\t%s\t%s\n' % (edge, seq, "PG:J:" + json_string, "RC:i:" + str(len(edge_dict['reads'])))]
        else:
            EdgeS += ['S\t%s\t%s\n' % (edge, seq)]
        Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (src, "+", edge, "+", "0M"))
        Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (edge, "+", dst, "+", "0M"))     

    with open(outputfilename, 'w') as fp:
        for h in header:
            fp.write(h)
        for s in AnchorS:
            fp.write(s)
        for s in EdgeS:
            fp.write(s)
        for l in Link:
            fp.write(l)

def write_graph_fromgraph(filename, graph):
    for edge in graph.edges:
        graph.edges[edge]['src'] = graph.incoming[edge][0]
        graph.edges[edge]['dst'] = graph.outgoing[edge][0]

    write_gfa(graph.anchor, graph.edges, filename)

# error correction -- underconstruction

def error_correction(graph, ref = "MHC-CHM13", threshold = 4):
    deleted_edge = []
    ###############################
    for edge in graph.edges:
        if ref in graph.edges[edge]['reads']:
            continue
        if len(graph.edges[edge]['reads']) < threshold:
            deleted_edge.append(edge)

    for edge in deleted_edge:
        edge_seq = graph.edges[edge]['seq']
        edge_readsets = graph.edges[edge]['reads']
        edge_samples = graph.edges[edge]['strain']
        ######### parallel path ###############
        src = graph.incoming[edge][0]
        dst = graph.outgoing[edge][0]
        read_sets = AGG.find_all_reads(graph)
        paths = AGG.Find_all_Path_between_anchors(graph, src, dst, read_sets)
        score = []
        for p, readsets in paths.subpath:
            if edge in p:
                score.append(99999999)
                continue
            seq = AGG.reconstruct_path_seq(graph, p)
            score.append(distance(edge_seq, seq))

        index = numpy.argmin(score)
        p, read_sets = paths.subpath[index]
        # print(p)
        for item in p:
            if item.startswith('E'):
                graph.edges[edge]['reads'] = list(set(graph.edges[edge]['reads'] + edge_readsets))
                graph.edges[edge]['samples'] = list(set(graph.edges[edge]['reads'] + edge_samples))

        ############## delete edge #################
    for edge in deleted_edge:    
        del graph.edges[edge]


def trim_graph(graph,ref = "MHC-CHM13", threshold = 4):
    ################### delete source and sink edges ############
    edgelist = graph.incoming.get("SINK", []) + graph.outgoing.get('SOURCE', []) + graph.incoming.get("SINK", []) + graph.incoming.get("SOURCE", [])
    for edge in edgelist:
        if edge in graph.edges:
            del graph.edges[edge]
            src = graph.incoming[edge][0]
            dst = graph.outgoing[edge][0]
            graph.outgoing[src].remove(edge)
            graph.incoming[dst].remove(edge)

    error_correction(graph, threshold)

class Get_Series_Parallel_Graph:
    
    def __init__(self, graph):
        self.initial_set = self.find_all_reads(graph)
        self.nodelist = self.series_parallel_graph_nodelist(graph)
        #print(self.nodelist)
        self.anchor, self.edges, self.outgoing, self.incoming = self.series_parallel_graph(self.nodelist, graph)

    def find_all_reads(self, graph):
        read_sets = set()
        edgelist = graph.edges.keys()
        for item in edgelist:
            readlist = graph.edges[item]['reads']
            for read in readlist:
                read_sets.add(read)
        return read_sets
    def find_furthest_node(self, node_candidate, subgraph, start_node):
        max_distance = -1
        node = ""
        for n in node_candidate:
            d = numpy.absolute(subgraph.anchor[n]['pos'] - subgraph.anchor[start_node]['pos'])
            if d > max_distance:
                node = n
        return node

    def series_parallel_graph_nodelist(self, subgraph): 
        #  need to consider the loop i.e.the last anchor to the first anchor

        start_node = sorted(subgraph.anchor.keys())[0]
        end_node = sorted(subgraph.anchor.keys())[-1]
        print(start_node, end_node)
        Nodelist = [start_node]

        edgelist = subgraph.outgoing[start_node]
        node_candidate = []
        for edge in edgelist:
            nodelist = subgraph.outgoing[edge]
            node_candidate += nodelist
            if nodelist[0] not in subgraph.anchor:
                    continue
        node_candidate = sorted(node_candidate)
        #node = node_candidate[-1] ## node should be selected by the largest distance instead of the sorting order
        node = self.find_furthest_node(node_candidate, subgraph, start_node)

        # find the furthurest anchor
        Nodelist.append(node) # append the furthest node
        


        while True:
            if node == end_node or node == '':
                break
            edgelist = subgraph.outgoing[node]

            node_candidate = []
            for edge in edgelist:
                nodelist = subgraph.outgoing[edge]
                # exclude deadend
                if "SINK" in nodelist:
                    continue
                if nodelist[0] not in subgraph.anchor:
                    continue
                if nodelist[0] not in subgraph.outgoing:
                    continue
                node_candidate += nodelist

            node_candidate = sorted(node_candidate)
            node = self.find_furthest_node(node_candidate, subgraph, node)
            if node in set(Nodelist):
                Nodelist.append(node)
                break
            Nodelist.append(node) # append the furthest node  
                
        return Nodelist
    
    def series_parallel_graph(self, Nodelist, subgraph):
        Node_dict = {}
        Edge_dict = {}
        Outgoing_dict = {}
        Incoming_dict = {}
        for i, node in enumerate(Nodelist[:-1]):
            start_node = node
            end_node = Nodelist[i+1]
            Node_dict[start_node] = subgraph.anchor[start_node]
            # print(start_node, end_node)
            path = AGG.Find_all_Path_between_anchors(subgraph, start_node, end_node, self.initial_set)
            #print(start_node,end_node, len(path.subpath))
            index = 0
            for p, rs in path.subpath:
                edgename = 'E%05d.%04d' % (int(start_node[1:]), index)
                seq = AGG.reconstruct_path_seq(subgraph, p[1:-1])
                Edge_dict[edgename] = {}
                Edge_dict[edgename]['seq'] = seq
                Edge_dict[edgename]['src'] = start_node
                Edge_dict[edgename]['dst'] = end_node
                Edge_dict[edgename]['reads'] = list(rs)
                Edge_dict[edgename]['strain'] = list(set([item.split("|")[-1] for item in list(rs)]))

                Outgoing_dict[start_node] = Outgoing_dict.get(start_node, []) + [edgename]
                Outgoing_dict[edgename] = [end_node]

                Incoming_dict[end_node] = Incoming_dict.get(end_node, []) + [edgename]
                Incoming_dict[edgename] = [start_node]
                index += 1
        return Node_dict, Edge_dict, Outgoing_dict, Incoming_dict
    

######################### execution ##############################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Merged and aligned fastas, file name")
    parser.add_argument("-r", "--reference", type=str, required=True,
                        help="reference sequence fasta of T2T")
    
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="output gfa prefix")
    parser.add_argument("-k", "--kmersize", type=int, required= True,
                        help="integer for kmer size")
    parser.add_argument("-s", "--genestart", type=int, required= True,
                        help="start coordinate of gene in T2T coordinates")
    parser.add_argument("-e", "--geneend", type=int, required= True,
                        help="end coordinate of gene in T2T coordinates")
    parser.add_argument("-t", "--threshold", type=int, required= True,
                        help="minimal number of reads supporting each path")
    
    args = parser.parse_args()
    filename = args.input
    reference = args.reference
    k = args.kmersize
    outputfilename = args.output
    gene_start = args.genestart
    gene_end = args.geneend
    threshold = args.threshold

    ref_HLA = find_reference_and_source_sink_anchors(reference, gene_start, gene_end)
    headers, sequences = loadFasta(filename)
    unique_kmerlist = get_reference_kmer_profile(ref_HLA, k)
    # Sample_reads = classify_sequences_by_samples(headers, sequences)
    count,HLA_samples,HLA_seq = find_sequences_between_sanchor_eanchor(headers, sequences, ref_HLA)
    SampleDict, PositionDict = map_reference_unique_kmers_to_seq(unique_kmerlist, HLA_seq, k)
    Final_anchor = get_anchor_information(SampleDict,HLA_samples, PositionDict, k)
    Edge_info, Outgoing = create_edgefile(HLA_seq, Final_anchor, k)
    
    write_gfa(Final_anchor, Edge_info, outputfilename + '.raw.gfa')
    graph = AGG.GraphicalGenome(outputfilename + '.raw.gfa')
    trim_graph(graph, threshold)
    write_graph_fromgraph(outputfilename + '.trimmed.gfa', graph)
    graph = AGG.GraphicalGenome(outputfilename + '.trimmed.gfa')
    series_parallele_graph = Get_Series_Parallel_Graph(graph)
    write_graph_fromgraph(filename=outputfilename + '.sp.gfa', graph=series_parallele_graph)


if __name__ == "__main__":
    main()