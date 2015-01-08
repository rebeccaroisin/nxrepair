#!/usr/bin/env python

# Copyright (c) 2014, Illumina
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import pysam
import collections
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from intervalNode import IntervalNode

def meansd(frq):

    """

    Function to calculate mean and standard deviation from a dictionary of frequencies.

    Return a tuple of (mean, std).

    Arguments:
    frq: a dictionary of frequencies: key = insert size, value = frequency

    """

    keys = frq.keys()
    keys.sort()
    w = np.empty(len(keys),np.float)
    for i,k in enumerate(keys):
        w[i] = frq[k]

    x = np.abs(keys)
    xbar = np.average(x,weights=w)
    xsd = np.average(np.sqrt(np.power((x - xbar),2)),weights=w)
    return (xbar,xsd)

def MAD(frq):

    """

    Function to calculate median and median absolute deviation (MAD) from a dictionary of frequencies.

    Return a tuple of (mean, MAD).

    Arguments:
    frq: a dictionary of frequencies: key = insert size, value = frequency
     
    """

    all_lengths = []
    for k, v in frq.iteritems():
        new_vals = [k] * int(v)
        all_lengths.extend(new_vals)
    all_lengths = np.array(sorted(all_lengths))
    mid = len(all_lengths)/2
    median = all_lengths[mid]
    residuals = sorted(abs(all_lengths - median)) # difference between val and median
    MAD = residuals[mid] # median of residuals
    print MAD
    return median, MAD

def find_intersections(tree, s, e):

    """

    Function to find inserts that bridge a region of a contig.

    Arguments:
    tree: interval tree of start/end positions of mate pair inserts
    s: start position of interval on contig 
    e: end position of interval on contig
     
    Return a list of nodes from tree whose start and end positions span the interval given by s and e. 
    
    """

    # find all reads that bridge a gap
    intersections = []
    tree.intersect(s, e, lambda x: intersections.append(x)) # see interval node for implementation
    return intersections

def get_insertlengths(reads):

    """

    Function to calculate interval sizes of a set of mate pairs.

    Arguments:
    reads: a list of mate pairs as interal tree node objects

    Return two numpy arrays: an array of insert sizes (integers) and an array of strand alignments (boolean)  
    """

    distances = []
    strands = []
    for read in reads:
        distances.append(read.end - read.start) # insert length
        strands.append(read.other[1]) # boolean: correct alignment
    return np.array(distances), np.array(strands)

def probability_of_readlength(read_length, mu, sigma, pi1, L):

    """

    Function to calculate the probability that mate pair insert sizes are not anomalous.
    Return an array of probabilites.

    Arguments:
    read_length: a numpy array of insert sizes
    mu: mean insert size (float)
    sigma: insert size standard deviation (float)
    pi1: prior probability of being anomalous (float)
    L: length of contig to which the reads in read_length are aligned

    """
    p_0 = pi1 * (1 / float(L)) # anomaly
    # probability of drawing from a gaussian with mean mu and std sigma
    p_1 = (1 - pi1) * stats.norm.pdf(read_length, loc=mu, scale=sigma) 

    p_total = p_1 / (p_0 + p_1)
    return p_total

class aligned_assembly:

    """

    Class to hold a set of mate pair or paired end reads aligned to the scaffolded genome assembly

    """    

    def __init__(self, bamfile, fastafile, min_size, threshold, step, window, minmapq, maxinsert, fraction, prior):

        """

        Initialiser function for the aligned assembly class.

        Arguments:
        bamfile: a sorted bam file of reads aligned to the scaffolded assembly
        fastafile: the scaffolded assembly in fasta format
        min_size: the minimum size contig to consider for analysis (integer)
        threshold: threshold in Z score below which a misassembly is called (float)
        step: step size to walk contigs (integer)
        window: width of window around each position from which mate pair insert sizes are fetched (integer)
        minmapq: the minimum mapq value for which paired reads are evaluated (float)
        maxinsert: the maximum insert size for which genome population statistics are calculated 
        fraction: minimum fraction of read pairs with correct orientation to call support for the assembly.

        """
        # initialising user parameters
        self.minmapq = minmapq
        self.maxinsert = maxinsert
        self.threshold = threshold
        self.step = step
        self.window = window
        self.fraction = fraction
        self.prior = prior
        self.sam = pysam.Samfile(bamfile, "rb" )
        self.fasta = pysam.Fastafile(fastafile)
        self.min_size = min_size

        # getting reads from bamfile
        self.all_reads = self.sam.fetch()
        self.references = self.sam.references
        self.lengths = self.sam.lengths
        
        # refdict: key=contig, val=contig length
        # read_stock: key=contig, val=aligned reads
        self.refdict = {}
        self.read_stock = {}
        for k,v in zip(self.references,self.lengths):
            self.refdict[k]=v
            self.read_stock[k] = self.get_reads(k, 0, v) 
        self.sizes = self.get_read_size_distribution()
        self.isize_median, self.isize_MAD = MAD(self.sizes)
        self.isize_mean, _ = meansd(self.sizes)
        self.isize_sd = 1.4826 * self.isize_MAD
        print self.isize_sd, self.isize_MAD, 1.4826 * self.isize_MAD

    def get_read_size_distribution(self):

        """

        Function to calculate global insert size distribution across the whole assembly
        Return a frequency table of insert sizes as a  dictionary with key = insert size, value = frequency

        """

        frq = collections.defaultdict(int) # dictionary of insert sizes
        found = {}
        for read in self.all_reads:
            # accept read based on mapq, contig alignemnt and insert size
            if (read.mapq > self.minmapq) and (read.rnext == read.tid) and (abs(read.tlen) < self.maxinsert):
                if read.qname in found and found[read.qname][0]==read.tid:
                    mate = found[read.qname]
                    isize = abs(max( mate[1]+mate[2]-read.pos,read.pos+read.rlen-mate[1] ))
                    frq[isize] += 1
                else:
                    found[read.qname] = (read.tid,read.pos,read.rlen)
        return frq

    def get_reads(self, ref, start, end):

        """

        Function to fetch reads aligned to a specific part of the assembled genome and return a list of aligned reads, where each list entry is a tuple:
        (read start position, read end position, read name, strand alignment) and strand alignment is a boolean indicating whether the two reads of a read pair align correctly to opposite strands.
        Reads are fetched that align to contig "ref" between positions "start" and "end".

        Arguments:
        ref: the name of the contig from which aligned reads are to be fetched.
        start: the position on the contig from which to start fetching aligned reads
        end: the position on the contig from which to end fetching aligned reads

        """

        # fetch all reads within a region
        # insert size: gap between end of one mate and start of next
        reads = self.sam.fetch(ref, start, end)
        read_stock = []
        found = {}
        for read in reads:
            if (read.rnext == read.tid):
                if read.qname in found and found[read.qname][0]==read.tid: # if mate maps to same contig
                    mate = found[read.qname] # fetch mate

                    # correctly ordering mates
                    if mate[1] > read.pos:
                        start_pos = read.pos + read.rlen
                        end_pos = mate[1]
                    else:
                        start_pos = mate[1] + mate[2]
                        end_pos = read.pos
                    # add mates to list of mates on that contig
                    # include strand orientation info
                    correct_strands = ((read.is_reverse) and not (read.mate_is_reverse)) or ((read.mate_is_reverse) and not (read.is_reverse))
                    read_stock.append((start_pos, end_pos, read.qname, correct_strands))
                else:
                    found[read.qname] = (read.tid,read.pos,read.rlen) # haven't reached mate yet
        return read_stock

    #@profile
    def make_tree(self, ref):

        """

        Function to construct an interval tree from reads aligning to a contig and return the interval tree.

        The interval tree stores nodes with properties start (start postition of interval), end (end position of interval) and other,
        which is a tuple of the mate pair name (string) and the strand alignment of the two paired reads (boolean).

        Arguments:
        ref: Reference ID of the contig for which the interval tree is to be constructed

        """

        bridges = self.read_stock[ref] 
        # insert first interval into tree
        s1, e1, name, correct_strands = bridges[0]
        tree = IntervalNode(s1, e1, other=(name, correct_strands))

        # insert the rest of the intervals
        for (start, end, name, correct_strands) in bridges[1:]:
            tree = tree.insert(start, end, other=(name, correct_strands))
        return tree

    def get_read_mappings(self, ref):

        """

        Function to calculate the fraction of reads pairs within a contig that align correctly to opposite strands.

        Return five arrays: the positions at which strand alignment was evaluated, the fraction correctly aligned, the fraction incorrectly aligned to the same strand, the unmapped 
        fraction and the fraction that have some other alignment issue.

        Arguments:
        ref: the reference id of the contig to be evaulated

        """
        dump_val = self.step
        positions = []
        same_strand = 0
        opp_strand = 0
        unmapped = 0
        other = 0

        # arrays of read mapping behaviour
        good_ratio = []
        unmapped_ratio = []
        bad_ratio = []
        other_ratio = []
        mini_pos = []

        reads = self.sam.fetch(reference = ref)
        # note that iterating in this manner works because the bam file is sorted.
        # create arrays containing fraction of correctly / incorrectly alinged reads
        for i, r in enumerate(reads):    
            mini_pos.append(r.pos)
            if r.mate_is_unmapped:
                unmapped += 1
            elif ((r.is_reverse) and not (r.mate_is_reverse)) or ((r.mate_is_reverse) and not (r.is_reverse)):
                same_strand += 1
            elif((r.is_reverse) and (r.mate_is_reverse)) or (not (r.mate_is_reverse) and not (r.is_reverse)):
                opp_strand += 1
            else:
                other += 1

            if (i+1) % dump_val == 0:
                total = same_strand + opp_strand + unmapped + other
                good_ratio.append(float(same_strand) / total)
                bad_ratio.append(float(opp_strand) / total)
                unmapped_ratio.append(float(unmapped) / total)
                other_ratio.append(float(other) / total)

                same_strand = 0
                opp_strand = 0
                unmapped = 0
                other = 0
                positions.append(np.mean(mini_pos))
                mini_pos = []

        return np.array(positions), np.array(good_ratio), np.array(bad_ratio), np.array(unmapped_ratio), np.array(other_ratio)

    def get_mapping_anomalies(self):

        """

        Function to determine the frequency of strand mapping anomalies across the entire genome assembly.

        Calls get_read_mappings for each contig larger than the aligned_assembly.min_size and returns:
        1) a dictionary with keys = contig reference IDs; values = list of positions and strand alignment ratios as described in get_read_mappings   
        2) a dictionary of anomalies wiht keys = contig reference IDs, values = [list of positions for which the ratio of correctly aligned strands < 0.75 (currently hard-coded), corresponding ratio of correctly aligned strands]

        """

        mapping_ratios = {} # key=contig, val=list of arrays of mapping behaviours
        anomalies = {}
        for w, (ref, length) in enumerate(self.refdict.iteritems()):
            if length > self.min_size: # consider only big contigs
                positions, good_ratio, bad_ratio, unmapped_ratio, other_ratio = self.get_read_mappings(ref)
                map_criterion = good_ratio < self.fraction
                pos_anomalies = positions[map_criterion]
                map_anomalies = good_ratio[map_criterion]
                mapping_ratios[ref] = [positions, good_ratio, bad_ratio, unmapped_ratio, other_ratio]
                anomalies[ref] = [pos_anomalies, map_anomalies]
        return mapping_ratios, anomalies

    def get_size_anomalies(self):

        """

        Function to determine the frequency of insert size anomalies across the entire genome assembly.

        Calls probability_of_readlength for each contig larger than the aligned_assembly.min_size and returns:
        1) a dictionary with keys = contig reference IDs; values = array of Zscores as described in probability_of_readlength   
        2) a dictionary of anomalies wiht keys = contig reference IDs, values = [list of positions for which  abs(z-score) > 2 (currently hard-coded), corresponding z-score value]

        """

        anomalies = {}
        zscores = {}
        all_probabilities = []
        stock_probabilities = {}

        for w, (ref, length) in enumerate(self.refdict.iteritems()):
            if length > self.min_size:
                tree = self.make_tree(ref) # build tree from all reads aligning to a contig
                positions = np.arange(self.step, length - self.window, self.step)
                probabilities = []
                for pos in positions:
                    if pos % 10000 == 0:
                        print pos
                    bridges = np.array(find_intersections(tree, pos-self.window, pos+self.window)) # fetch reads in windows across contig
                    bridge_lengths, strand_alignment = get_insertlengths(bridges) # get insert sizes and mapping behaviour
                    prob_lengths = probability_of_readlength(bridge_lengths, self.isize_mean, self.isize_sd, self.prior, length) # get prob. insert sizes from null
                    condition = strand_alignment == 1
                    D = np.sum(prob_lengths[condition]) # D is total assembly support
                    probabilities.append(D)
                    all_probabilities.append(D)
                stock_probabilities[ref] = [positions, np.array(probabilities)]
        p_mean = np.mean(np.array(all_probabilities)) # get contig mean and variance
        p_std = np.std(np.array(all_probabilities))

        for ref, [positions, probs] in stock_probabilities.iteritems():
            zscore = (probs - p_mean) / p_std # calculate position z score from contig mean, std
            # anomalies have Zscore < Threshold.
            # Note: threshold should be negative
            z_criterion = (zscore < self.threshold)
            z_anomalies = zscore[z_criterion]
            print ref, z_anomalies
            pos_anomalies = positions[z_criterion]
            zscores[ref] = [positions, zscore]
            anomalies[ref] = [pos_anomalies, z_anomalies] # list of anomaly locations and socres
        return zscores, anomalies, tree

    def get_anomalies(self, outfile, trim, img_name=None):

        """

        Function to determine the frequency of anomalous mate pair behaviour across the entire genome assembly and return a dictionary where:
        key = contig reference IDs, 
        value = list of postions within that contig where an assembly error is identified and the contig should be broken.

        Calls get_size_anomalies and get_mapping_anomalies for each contig larger than the aligned_assembly.min_size; makes a .csv file listing for each contig the positions of identified misassemblies and their corresponding anomalous scores.

        Arguments:
        outfile: name of file (including filepath) to store the list of contig misassemblies.

        Keyword Arguments:
        img_name: name of file (including filepath, not including filetype) to store plots of alignment quality

        """

        print "Anomaly detection"
        # get anomaly positions
        zscores, size_anomalies, tree = self.get_size_anomalies()
        map_ratios, map_anomalies = self.get_mapping_anomalies()

        break_points = {}

        # # make a wiggle file
        # print "Writing wiggle file"
        # wig_file = "%s.wig" % ("/media/rmurphy/sandbox/bash_scripts/test_TB")
        # with open(wig_file, "w") as wig:
        #     #wig.write("track type=wiggle_0 graphType=line color=0,0,255 altColor=255,0,0 name='Zscore' graphType=heatmap midRange=35:65 midColor=0,255,0\n")
        #     for ref, [positions, probs] in zscores.iteritems():
        #         print ref
        #         wig.write("fixedStep chrom=%s start=%s step=%s span=%s\n" % (ref, 1, self.step, 1)) 
        #         #print zscores[ref]
        #         #for vals in zscores[ref]:
        #         #    positions = vals[0]
        #         #    probs = vals[1]
        #         for prob in probs:
        #             wig.write("%s\n" % (prob))
        #         wig.write("\n")

        for w, (ref, [positions, probs]) in enumerate(zscores.iteritems()):
            # write all Z scores to a csv file
            print "Writing Z scores to file (%s)" % ref
            for pos, prob in zip(positions, probs):
                outfile.write("%s %s %s\n" %(ref, pos, prob))
            z_pos, z_anomalies = size_anomalies[ref] 
            print "z_anomalies:", z_anomalies
            print ref, z_pos
            map_positions, good_ratio, bad_ratio, unmapped_ratio, other_ratio = map_ratios[ref] 
            pos_anomalies, map_anom = map_anomalies[ref]

            anomaly_positions = sorted(z_pos.tolist())
            print ref, anomaly_positions

            # make list of positions to break this contig
            break_points[ref] = []
            
            if len(anomaly_positions) != 0:
                current = []
                for p in range(len(anomaly_positions)):
                    anom_pos = anomaly_positions[p]
                    if current == []:
                        current.append(anom_pos)
                    else:
                        if anom_pos - current[-1] <= trim:
                        # if anomalies are well separated flush current values to break_point
                            current.append(anom_pos)
                        else:
                            break_points[ref].append(np.mean(current))
                            current = [anom_pos]
                if current != []:
                    break_points[ref].append(np.mean(current))
                print ref, break_points[ref]

            if img_name != None:
                # plot zscores and anomalies
                fig, ax1 = plt.subplots()
                plt.subplots_adjust(bottom=0.15)
                ax1.set_xlim([0, max(map_positions)])
                ax1.set_xlabel("Position",size=24)
                ax1.set_ylabel("Assembly Support", size=24)
                plt.tick_params(axis='both', which='major', labelsize=20)
                lns1 = ax1.plot(positions, probs, c="k", label="Support")
                #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                anomalies_to_plot = sorted(break_points[ref])
                anomalies_y = [-10] * len(anomalies_to_plot)

                ax1.scatter(anomalies_to_plot, anomalies_y, c="r", marker = "o")
                print "anomalies", anomalies_to_plot

                # if name given, save image as .pdf and .png
                name = img_name + "_%s.pdf" % ref
                plt.savefig(name)
                name = img_name + "_%s.png" % ref
                plt.savefig(name)
                plt.cla()
        
        return break_points

    def breakContigs_double(self,outfile, breakpoints, trim):

        """

        Function to break a contigs at positions identified as assembly errors and write a new fasta file containing all contigs (both altered and unaltered).

        Makes a two-point break at the identified misassembly position, splitting at 5 Kb upstream and downstream of the misassembly and (currently) excluding the misassembled region.

        Arguments:
        outfile: name of the new fasta file (including filepath)
        breakpoints: dictionary of misassemblies. key = contig reference ID, value = list of misassembly positions within the contig
        trim: distance, in bases, to trim from each each edge of a breakpoint to remove misassembly (integer)

        """
        
        for k, v in breakpoints.iteritems():
            # print breakpoints
            print k, v
        newcontigs = []
        for contig, length in self.refdict.iteritems():
            #dna = self.fasta[contig] # sequence of contig
            dna = self.fasta.fetch(reference=contig) # sequence of contig
            if contig in breakpoints:
                splits = breakpoints[contig]
                splits.sort()
                prev = 0
                for s in splits: # iterate through breakpoints
                    print s
                    print "breaking contig"
                    if (s - prev > trim) and ((length - s) > trim): 
                        print "breaking here:", s
                        newcontigs.append(dna[int(prev):int(s-trim)]) # trim and append section before break
                        prev = s + trim # trim other end of break
                    else:
                        print "Too small!"
                newcontigs.append(dna[int(prev):])
            else:
                newcontigs.append(dna)

        # write new contigs to file
        newcontigs.sort(lambda x,y: cmp(len(x), len(y)),reverse=True)
        print "Writing new fasta..."
        for count, dna in enumerate(newcontigs):
            name = ">CONTIG_%d_length_%d"%(count,len(dna))
            print name
            outfile.write(name)
            outfile.write("\n")
            outfile.write(dna)
            outfile.write("\n")    

def main():
    # read command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Routine to identify and correct large-scale misassemblies in de novo assemblies')
    parser.add_argument('bam', metavar='bam', type=str, help='bam')
    parser.add_argument('fasta', metavar='fasta', type=str, help='scaffold fasta')
    parser.add_argument('outfile', metavar='outfile', type=str, help='Output file name')
    parser.add_argument('newfasta', metavar='newfasta', type=str, help='Fasta file for new contigs, including filepath')
    parser.add_argument('-min_size', metavar='min_size', type=int, default=10000, help='Minimum contig size to analyse')
    parser.add_argument('-img_name', metavar ='img_name', type=str, default=None, help='Name under which to save (optional) graphs of alignment quality. Default value: None (no graphs produced)')
    parser.add_argument('-trim', metavar ='trim', type=int, default=4000, help='Number of bases to trim from each side of an identified misassembly. Default value: 5000')
    parser.add_argument('-T', metavar ='T', type=float, default= -4.0, help='Threshold in Z score below which a misassembly is called. Default value: -4.0')
    parser.add_argument('-step_size', metavar ='step_size', type=int, default=1000, help='Step-size in bases to traverse contigs. Default value: 1000')
    parser.add_argument('-window', metavar ='window', type=int, default=200, help='Window size across which bridging mate pairs are evaluated. Default value: 200')
    parser.add_argument('-minmapq', metavar ='minmapq', type=int, default=40, help='Minimum MapQ value, above which a read pair is included in calculating population statistics. Default value: 40')
    parser.add_argument('-maxinsert', metavar ='maxinsert', type=int, default=30000, help='Maximum insert size, below which a read pair is included in calculating population statistics. Default value: 30000')
    parser.add_argument('-fraction', metavar ='fraction', type=int, default=0.75, help='Minimum fraction of read pairs with correct orientation to call support for the assembly. Default value: 0.75')  
    parser.add_argument('-prior', metavar ='prior', type=float, default=0.01, help='Prior probablility that the insert size is anomalous. Default value: 0.01')  


    args = parser.parse_args()

    # make assembly object
    f = aligned_assembly(args.bam, args.fasta, args.min_size, args.T, args.step_size, args.window, args.minmapq, args.maxinsert, args.fraction, args.prior)
    print "This is ok"
    
    # find anomalies
    with open(args.outfile, "w") as of:
        bps = f.get_anomalies(of, args.trim, args.img_name)

    # break contig at identified anomalies
    with open(args.newfasta, "w") as outfasta:
        f.breakContigs_double(outfasta, bps, args.trim)

if __name__ == "__main__":
    main()
