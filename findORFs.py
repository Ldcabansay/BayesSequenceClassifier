#!/usr/bin/env python3

import sequenceAnalysis 
import sys
import operator


class findORFs:
    
        
    def __init__ (self):
        '''constructor: creates orfFinder object from orfFinder class in sequenceAnalysis.py'''
        self.orfs = sequenceAnalysis.orfFinder()

    def findSeqOrfs(self, header, thisSeq, minOrf, maxOrf, startCodons, stopCodons, longestOrf):
        '''Method: passes args to self.orfs.findOrf method in order to find orfs for each header, sequence
        from FASTAfile'''
        self.orfs.findOrfs(header, thisSeq, minOrf, maxOrf, startCodons, stopCodons, longestOrf)
        

    def sortORF(self):
        '''Method: calls genomeORF dictionary from self.orfs. GenomeORFs contains all orfs found 
            in fasta file, each header is used as a key for each sequence. Using header key, each 
            orf in sequence is sorted in order of longest to shortest orf. The sorted orfs as a list are
            added to allorfs dict with key header.
        '''
        genomeOrfs = self.orfs.genomeORFs
        allorfs = {}
        for thisheader in genomeOrfs:
            orflist = []
            for frame in genomeOrfs[thisheader]:
                #orflist = []
                for orf in genomeOrfs[thisheader][frame]:
                    orflist.append(orf)
                    #allorfs[thisheader].append(orf)
            if orflist:
                sortedorfs = (sorted(orflist, key = lambda x: x[-1], reverse=True))
                allorfs[thisheader] = sortedorfs
        #print allorfs
        return allorfs            
       
    def writeOrf(self, header, file):
        """calls dictionary of sorted orfs from sortORFs() and writes to file the header and
            orfs found for each sequence in FASTA file.
        """
        thisheader = header
        allorfs = self.sortORF()
        
        if allorfs[thisheader]:
        #only write the sequences with orfs, 
        #if no orfs found and dictionary is empty it won't write to file
            seqOrfs = allorfs[thisheader]
            file.write(header)
            #print (seqOrfs)
            for i in seqOrfs:
                printorf = i
                #if i[0]=='-3':
                
                file.write("{x[0]:s} {x[1]:>5d}...{x[2]:>5d} {x[3]:>5d}\n".format(x = printorf))

    def printOrf(self, header, file):
        """calls dictionary of sorted orfs from sortORFs() and writes to file the header and
            orfs found for each sequence in FASTA file.
        """
        thisheader = header
        allorfs = self.sortORF()
        
        if allorfs[thisheader]:
        #only write the sequences with orfs, 
        #if no orfs found and dictionary is empty it won't write to file
            seqOrfs = allorfs[thisheader]
            #file.write(header)
            print header
            #print (seqOrfs)
            for i in seqOrfs:
                printorf = i
                #if i[0]=='-3':
                
                print ("{x[0]:s} {x[1]:>5d}...{x[2]:>5d} {x[3]:>5d}\n".format(x = printorf))
                        

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog -User can provide defined options for findORFs.py program to specify minimum orf size, reporting only longest orfs, or to use only specific start/stops codons. See usage for how to set options on commandline', 
                                             epilog = 'input filename and desired output filename must be specified according to usage syntax', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s inFileName [options] -option1[default]'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', default='output', nargs='?', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='prints only the longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (0,100,200,300,500,1000), default=100, action = 'store', help='specifies minimum number of bases for ORF gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', help='specified start Codon(s) used to define where ORF begins') #allows multiple start list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='specified stop Codon(s) used to define where ORF ends') #allows multiple stop list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def main(inCL=None): 
   
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)
    print (myCommandLine.args)
    # myCommandLine.args.inFile has the input file name
    fasta_file = myCommandLine.args.inFile
    sqReader = sequenceAnalysis.FastAreader(fasta_file)
    myGenome = findORFs()
    # myCommandLine.args.minGene is the minimum Gene length to include
    minOrf=myCommandLine.args.minGene
    maxOrf = 0
    # myCommandLine.args.longestGene is True if only the longest Gene is desired
    longestOrf=myCommandLine.args.longestGene
    # myCommandLine.args.start is a list of start codons
    startCodons=myCommandLine.args.start
    # myCommandLine.args.start is a list of stop codons
    stopCodons=myCommandLine.args.stop
    # myCommandLine.args.outFile has the output file name
    headers = []
    #print 'input'
    for header, sequence in sqReader.readFasta():
            #iterate through sequences from fastafile in sqReader generator
        myGenome.findSeqOrfs(header, sequence, minOrf, startCodons, stopCodons, longestOrf)
        headers.append(repr(header))
    if myCommandLine.args.outFile=='output':
        filename= 'results_'+repr(myCommandLine.args.inFile)+'.txt'
    else:
        filename = myCommandLine.args.outFile
    file = open(filename, "w")   
    for header in headers:
        myGenome.writeOrf(header, file)
    file.close
    
if __name__ == "__main__":
    main()  # delete the list when you want to run normally

#'tass2.fa','tass2ORFdata-ATG-100.txt','--longestGene'