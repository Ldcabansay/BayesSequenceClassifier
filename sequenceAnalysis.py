#!/usr/bin/env python3
# Name: Louise Cabansay (lcabansa)
# Group Members: List full names (CATS usernames) or None
import sys
import operator

class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O ProteinParam.aa2abs280[]
 
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
 
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
 
    def __init__ (self):
	    #'''contructor: saves dictionaries for amino acid composition (aaComp), nucleotide composition (nucComp), and codon composition (codonComp) '''
        self.aaComp = {
            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0,
            '-': 0}
        self.nucComp = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0}
        
        self.codonComp = {
           
                'UUU': 0, 'UCU': 0, 'UAU': 0, 'UGU': 0, # UxU
                'UUC': 0, 'UCC': 0, 'UAC': 0, 'UGC': 0, # UxC
                'UUA': 0, 'UCA': 0, 'UAA': 0, 'UGA': 0, # UxA
                'UUG': 0, 'UCG': 0, 'UAG': 0, 'UGG': 0, # UxG
                # C
                'CUU': 0, 'CCU': 0, 'CAU': 0, 'CGU': 0, # CxU
                'CUC': 0, 'CCC': 0, 'CAC': 0, 'CGC': 0, # CxC
                'CUA': 0, 'CCA': 0, 'CAA': 0, 'CGA': 0, # CxA
                'CUG': 0, 'CCG': 0, 'CAG': 0, 'CGG': 0, # CxG
                # A
                'AUU': 0, 'ACU': 0, 'AAU': 0, 'AGU': 0, # AxU
                'AUC': 0, 'ACC': 0, 'AAC': 0, 'AGC': 0, # AxC
                'AUA': 0, 'ACA': 0, 'AAA': 0, 'AGA': 0, # AxA
                'AUG': 0, 'ACG': 0, 'AAG': 0, 'AGG': 0, # AxG
                # G
                'GUU': 0, 'GCU': 0, 'GAU': 0, 'GGU': 0, # GxU
                'GUC': 0, 'GCC': 0, 'GAC': 0, 'GGC': 0, # GxC
                'GUA': 0, 'GCA': 0, 'GAA': 0, 'GGA': 0, # GxA
                'GUG': 0, 'GCG': 0, 'GAG': 0, 'GGG': 0  # GxG
            } #rna comp dict
    
    def addSequence (self, thisSeq):
	'''method: iterates over new sequence appending counts of nucleotides, codons, and amino acids to relative dictionary'''
        for nuc in thisSeq.upper():
            if nuc in self.nucComp:
                self.nucComp[nuc] += 1.0
        rnaSeq = thisSeq.replace('T','U')
        #print rnaSeq
        for n in range(0,len(rnaSeq),3):
            codon = rnaSeq[n:n+3]
            #print codon
            if codon in self.codonComp:
                self.codonComp[codon] += 1.0
                aa = NucParams.rnaCodonTable[codon]
                if aa in self.aaComp:
                    self.aaComp[NucParams.rnaCodonTable[codon]] += 1.0
                
    def aaComposition(self):
        return self.aaComp
    def nucComposition(self):
        return self.nucComp
    def codonComposition(self):
        return self.codonComp
    def nucCount(self):
        return sum(self.nucComp.values())


class FastAreader :

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
        
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
            
            idx = 0
            for line in fileH:
                #commented out idx lines used to track processing of large files, prints progress of reading file every NNNN lines using idx%NNNN 
                idx += 1
                if idx % 10000== 0:
                    print("Step " + repr(idx))
                if line.startswith ('>'):
                    yield header,sequence
                    #print len(sequence)
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header,sequence
        #print 'length sequence', len(sequence)

class orfFinder :

    def __init__ (self):
        '''Constructor: initializes lists and dictionaries to be used by the orfFinder methods 
        '''
        self.orfComp = {'+1':[], '+2':[], '+3':[], '-1':[], '-2':[], '-3':[] }
        self.genomeORFs = {}
            #instantiates a dictionary of all orfs for genome (entire fasta file)
            #this dictionary will use header as key value each header corresponds to one sequence
        self.stopCodons = []
        self.startCodons = []

    def findOrfs(self, header, thisSequence, minOrf, maxOrf, startCodons, stopCodons, longestOrf):
        '''Method: Main director of function
        '''  
        for codon in startCodons:
            self.startCodons.append(codon)
        for codon in stopCodons: 
            self.stopCodons.append(codon)
        self.longestOrf = longestOrf
        if maxOrf>0:
            self.maxOrf = maxOrf
        else:
            self.maxOrf = len(thisSequence)

        # forward pass
        self.findORFsInOneDirection(thisSequence, minOrf, False)

        # reverse pass
        revseq = self.reverseComp(thisSequence)
        self.findORFsInOneDirection(revseq, minOrf, True)

        self.genomeORFs[repr(header)] = self.orfComp
        #print 'ORF COMP:    ', self.orfComp
        #for each seq (new header in fasta file, append orfComp to genome comp and then clear orfcomp dictionaries)
        self.orfComp = {'+1':[], '+2':[], '+3':[], '-1':[], '-2':[], '-3':[] }

    def findORFsInOneDirection(self, seq, minOrf, isReverse):
        '''Method: finds orfs for given sequence. Appends orfs found to orfcomp dictionary.
                 Handles both forward and reverse strands depending on boolean value of isReverse
            inputs: sequence: sequence as string
                   minOrf: smallest orf length considered
                   isReverse: boolean value (T/F)
        '''
        orflist = []
        lenSeq = len(seq)
#        print(seq, isReverse)
        for frame in range(0,3):
            # loop over the list of start positions returned by findStartPositions:
            lastStop = -1  # it will keep track of the last stop seen
            startPositions = self.findStartPositions(seq, frame)
            for startPos in startPositions:
                #loop over start positions to find the stop position for each start
                if self.longestOrf and startPos < lastStop:
                    continue
                isDanglingStop = (startPos == 0 and frame > 0)
                #handles cases with dangling stop
                if isDanglingStop:
                    stopPos = self.findNextStopCodon(seq, frame)
                else:
                    stopPos = self.findNextStopCodon(seq, startPos)
                lastStop = stopPos
                orflength = stopPos - startPos + 1
#                print(frame+1, startPos+1, stopPos+1, isDanglingStop, lastStop)

                if orflength <= self.maxOrf:
                    if orflength >= minOrf:
                        #append orf to OrfComp dictionary in forward strand coordinates
                        orfsequence = seq[startPos+1:stopPos+1]
                        if isReverse:
                            strFrame = "-" + str(frame+1)
                            beginCoord = lenSeq - stopPos
                            endCoord = lenSeq - startPos
                            self.orfComp[strFrame].append((strFrame, beginCoord, endCoord, orflength, orfsequence))
                        else:
                            strFrame = "+" + str(frame+1)
                            beginCoord = startPos + 1
                            endCoord = stopPos + 1
                            self.orfComp[strFrame].append((strFrame, beginCoord, endCoord, orflength, orfsequence))
            
    def findStartPositions(self, seq, frame):
        '''Method: Called by findORFsInOneDirection() to find start positions for each frame
            inputs: sequence, frame
            output: list of start positions for given frame
        '''
        # store start positions found
        startlist = []
        # loop over all nucleotides in sequence except the last two (where codons cannot start):
        for nuc in range(frame, len(seq)-2, 3):
            codon = seq[nuc:nuc+3]
            # check if the codon is the start codon:
            if codon in self.startCodons:
                #append the start position as first base of start codon
                startlist.append(nuc)  
            # handles dangling stop, where a stop codon is encountered before start
            if codon in self.stopCodons:
                if len(startlist)<1:
                    startlist.append(0)
        # return the list
        return startlist


    def findNextStopCodonHelper(self, seq, startPos):
        '''Method: Helps findNextStop() by iterating through sequence to find stop position for a given start, 
            ranging from the position of given start to len(seq)-2. Ignores nucs/codons before
            start position. 
        '''
        # iterate sequence from specified start (position of start codon/orf start)
        # from start to, but not including, len(seq), in jumps of 3    
        for nuc in range(startPos, len(seq)-2, 3):
            # check if the current position is the start codon we look for:
            if seq[nuc:nuc+3] in self.stopCodons:
                #return the position of the last base of the stop codon,
                # this position is the end of the orf
                return nuc+2
        return None


    def findNextStopCodon(self, seq, startPos):
        '''Method: Uses findNextStopCodonHelper() to find stop position for current given start.
            If no start is found, then there is a dangling start without a stop so append a stop at 
            len(seq) to signify a downstream stop for the orf. 
        '''
        #print 'start pos', startPos
        #print 'seq', seq
        #seq = seq.upper()
        # stoplist store temporary restults:
        stopPos = []
        # find the start position of the next stop codon:
        pos = self.findNextStopCodonHelper(seq, startPos)
    	# check that pos is not None, which would indicate that the codon was not found:
        if pos != None:
            return pos
        #if no stop is found for the start, add a stop to end of sequence
        else:
            return len(seq) - 1

    def reverseComp(self, thisSeq):
        '''Method: translates sequence into its reverse complement to produce bottom strand
        '''
        seq = thisSeq.upper()
        complementComp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join([complementComp[base] for base in reversed(seq)])

    #def findORFseq(self, thisSeq, beginPos, endPos):
    