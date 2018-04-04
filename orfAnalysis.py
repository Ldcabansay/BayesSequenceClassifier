#!/usr/bin/env python3
import pandas as pd
from genomeAnalyzer import genomeAnalyzer as ga
import sequenceAnalysis as sa
import readwriteExcel
import sys
import operator
import numpy as np
from findORFs import findORFs

'''
orfAnalysis.py is a program developed to use FASTA files to 
generate an excel dataset of ORF relative codon frequenies

Program Usage: See class CommandLine() for how to call program

Directory Requirements:
    orfAnalysis must be placed in directory containing the fasta files.
    fasta file names must be in format: (#sample)_speciesName.fna
    Example: 3_Brucella.fna
        ^This is the fasta file name for the third sample of the Brucella genome
    
    Currently the main function of the program must be altered to 
    provide species names as a list. Ex: species = ['Escherichia','Brucella', 'MycobacteriumPhage']

Dependencies: 
    orfAnalysis requires the following programs to be within the same directory as orfAnalysis:
        sequenceAnalysis.py, genomeAnalyzer.py, findORFs.py readwriteExcel.py

'''

class orfAnalysis:
    '''
        Creates an orfAnalysis object that can be used to parse the ORFs of DNA sequences and computes the
        the relative codon frequencies for each ORF parsed. ORF min and max can be 
        specified by the user 

        instantiation: 
            myorfAnalysis = orfAnalysis()

        Usage:
            orfAnalysis(fasta_file, sampleCount, species, minOrf, maxOrf)
    '''

    def __init__ (self, file, sampleCount, species, minOrf, maxOrf):
        '''
        constructor: Initializes lists and dictionaries used by the program. 
                    Constructor creates a findORFs() object to construct a 
                    dictionary containing the orfparameters a the fasta file
                    ORF parameters: orf frame, begin pos, end pos, orf length, orf sequence
                    
        '''
        
        self.species = species
        self.count = sampleCount        
        self.orfData = {}
        self.orfs = findORFs()
        self.headers = []
        #setting maxOrf to zero will instead not set a maximum Orf size
        startCodons=['ATG']
        stopCodons=['TAG','TGA','TAA']
        longestOrf=True


        sqReader = sa.FastAreader(file)
        for header, sequence in sqReader.readFasta():
            if maxOrf==0:
                maxOrf = len(sequence)
            self.headers.append(repr(header))
            self.orfs.findSeqOrfs(header, sequence, minOrf, maxOrf, startCodons, stopCodons, longestOrf)

    def orfHeaders(self):
        #returns the fastafile headers for each species
        return self.headers

    def analyseOrf(self):
        '''Method: Sorts ORFs found and passes sequence of each orf to the geOrfData function.
            getORfData returns the relative codon frequencies of each ORF. The method constructs 
            a dataframe and adds the relative codon frequencies of each ORF to that dataframe and 
            returns dataframe at the end of function call. 
        '''
        orfDF = pd.DataFrame()
        genomeORFs = self.orfs.sortORF()
        orfSequences = []
        orfSeqLengths = []
        orfCount = 0
        for orf in genomeORFs[self.headers[0]]:
            orfCount +=1
            #gen1orf = ['Genome1: ']
            orfSeq = orf[-1]
            
            orfSequences.append(orfSeq)
            self.getOrfData(orfSeq)
            orfSeqLengths.append(float(len(orfSeq)))
            self.orfData['species'] = self.species
            #cookie=np.array(self.orfData)
            df = pd.DataFrame(self.orfData, index=[0])
            orfDF = pd.concat([orfDF, df])
            self.orfData = {}

        #print orfDF
        #print orfDF.columns.values
        return orfDF

        
    def getOrfData(self, ORFsequence):
        '''Method: passes ORF sequence to sequence analysis object to get relative codon frequencies
            Creates new orfData dictionary entry with relative codon frequency value
        
        '''
        orfInfo = ga()
        orfInfo.addSequence(ORFsequence)
        orfGC = orfInfo.gcContent()
        orfcodonComp = orfInfo.codonAnalysis()
        self.orfData['GC'] = float(orfGC)
        for codon in range(len(orfcodonComp)):
            #avgRelCodonfreq = float()/numOrfs
            self.orfData[orfcodonComp[codon][0]] = orfcodonComp[codon][2]
    '''
    def orfDataFrame(self): 
        #prints results to commandline    
        df = self.orfDataFrame_all
        if self.count<2:
            print df
        else:
            #print data[:,1]
            print 'NONE'
    '''

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
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog -User can provide defined options for findORFs.py program to specify minimum orf size, reporting only longest orfs, or to use only specific start/stops codons. See usage for how to set options on commandline', 
                                             epilog = 'input filename and desired output filename must be specified according to usage syntax', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s outfile.xlsx [options] -option1[default]'
                                             )
        
        self.parser.add_argument('outFile', action = 'store', default='output', nargs='?', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='prints only the longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, default=100, action = 'store', help='specifies minimum number of bases for ORF gene length')
        self.parser.add_argument('-xG', '--maxGene', type=int, default=0, action = 'store', help='specifies max number of bases for ORF gene length')
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

    data = pd.DataFrame()
    headers = []
    species = ['Escherichia', 'MycobacteriumPhage']
    #species = ['Escherichia','Brucella', 'MycobacteriumPhage']
    #species = ['Escherichia','Brucella']
    #species = ['MycobacteriumPhage','Brucella']
    sampleCount = 0
    #species = 'Escherichia'

    minOrf=myCommandLine.args.minGene
    maxOrf=myCommandLine.args.maxGene


    for species in species:
        for i in range(0,24):
            count = i+1
            sampleCount +=1
            print 'sampleCount', sampleCount
            fasta_file = repr(count)+'_'+species+'.fna'
            mygenomes = orfAnalysis(fasta_file, sampleCount, species, minOrf, maxOrf)
            orfDF = mygenomes.analyseOrf()
            data = pd.concat([data, orfDF],axis=0)
            headers.append(mygenomes.orfHeaders())
    #print writedata.shape
    rwExcel = readwriteExcel
    sheetname = repr(len(species))+'EschBruc2'
    
    if myCommandLine.args.outFile != 'output':
        excelfile = myCommandLine.args.outFile
    else:
        excelfile=sheetname+'.xlsx'

    rwExcel.writeExcelData(data,excelfile,sheetname,startrow=1,startcol=1)

if __name__ == "__main__":
    main()