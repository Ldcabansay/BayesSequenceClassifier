import os, struct
import matplotlib as plt
import numpy as np
import numpy.linalg as LA
import pandas as pd
import random
import operator
import readwriteExcel
import plotly
from plotly.graph_objs import *
import plotly.tools as tls


class ClassData:
    '''
        Defines an object to hold dataset, and also the objects trainingSet and testSet.
        trainingSet and testSet are from the randomly split dataset.
        User can define split as the percent of data to be trainingSet (as the decimal value of percentage)
        User can also define a specific randomSplit to duplicate results
 
        instantiation:     
            myData = ClassData()    

        usage:
            dataset = 'datafile' #map dataset to datafile
            split = 0.75 #percent of dataset to be training set (enter as decimal)
            randomSeed = 0 # for this method, if randomSeed = 0, a random value will be called to split the data
                    #set randomSeed to a specific integer in order to duplicate a specific random split result 
            [allDataset, TrainDataset, TestDataset] = myData.TestTrainDataSplit(dataset, split, randomSeed) 

    '''
    def __init__(self):
        '''Constructor: initializes arrays for the whole dataset, trainingSet, and testSet
        '''
        self.dataset = []
        self.trainingSet = []
        self.testSet = []
        
    def TestTrainDataSplit(self, dataset, split, randomSeed):
        '''Method: splits dataset into a randomized test set and training set
        '''
        #full dataset
        self.dataset = dataset
        if randomSeed==0:
            random.seed()
        else:
            random.seed(randomSeed)
        for x in range(len(dataset)):
                dataset[x] = dataset[x]
                if random.random() < split:
                    self.trainingSet.append(dataset[x])
                else:
                    self.testSet.append(dataset[x])
        return np.array(self.dataset), np.array(self.trainingSet), np.array(self.testSet)

    def dataLabels(self):
        '''Method: constructs arrays of the feature vector labels according to each of the datasets (all, train, test)
        '''
        TrainDataset = np.array(self.trainingSet)
        TestDataset = np.array(self.testSet)
        self.trainlabels = (TrainDataset[:,-1][...][None]).astype(str)
        #print 'train labels', self.trainlabels.shape
        self.testlabels = (TestDataset[:,-1][...][None]).astype(str)
        #print 'test labels', self.testlabels.shape
        self.labels = list(set(self.trainlabels[0,:]))
        self.classlabels = np.array(self.labels)
        print 'Classes: ', self.labels, '\n'
        
        return self.classlabels, self.trainlabels, self.testlabels

class PCA:
    '''
    PCA contains the methods to conduct a principle components analysis on a given data set. 

    instantiation:
        myPCA =  XUZCVPR()
    usage:
        XUZCVPR(self, dataset, Utrain, train)
    '''

    def XUZCVPR(self, dataset, Utrain, train):
        '''Method: Main PCA function that calculates principle components of the data set using a
            convariance matrix 
        '''
        #X is the matrix of original feature vectors (column w/ feature labels removed, but rows are preserved for class identification later)
        #future testing note: might be possible to preserve feature labels if dataframe was used instead of np.array, but np functions might not work
        self.X = np.array(dataset)
        if train==True:
            Uvector = np.mean(self.X,axis=0) #U is mean vector (the mean of each feature)
            self.U = np.array([Uvector])
        else: 
            self.U = Utrain 
        self.Z = self.X - self.U #Z is the difference vector after you subtract the mean (U) from each feature in all feature vectors
    
        self.meanZ = np.mean(self.Z,axis=0) # axis to calculate column means1
        self.meanZround = [round(x) for x in self.meanZ] #test to see if calculated mean correctly
        self.emptymeanZ=filter(lambda x:x != 0, self.meanZround) # all the column mean of z should be 0
    
        self.C = np.cov(self.Z.astype(float),rowvar=False) #covariance matrix
        self.Ctranspose = self.C.transpose()
        self.checkC = np.array_equal(self.C,self.Ctranspose)
    
        self.aEighV=LA.eigh(self.C)

        self.V = np.flipud(self.aEighV[1].T) #normalized eigen vectors
        self.Evals = np.flipud(self.aEighV[0]) #eigen values
        self.Vrows = self.V[0,:]
        self.checkVrows = (np.dot(self.C, self.Vrows))/(self.Evals[0]*self.Vrows)
        self.eig_pairs = [(np.abs(self.Evals[i]), self.V[:,i]) for i in range(len(self.Evals))]
        self.P=np.dot(self.Z,self.V.T) #Construct projection matrix for principal components
        self.R=np.dot(self.P,self.V) #R is analogous to Z vector (Xrec-U)
        self.Xrec = self.R+self.U #add U to R and see if Xrec is similar to original X matrix
    
        return np.array(self.X), np.array(self.U), np.array(self.Z), np.array(self.C), np.array(self.V), np.array(self.Evals), np.array(self.P), np.array(self.R), np.array(self.Xrec)
    
    def PCAcheck(self):
        '''Method: prints summary of PCA calculations to check if done correctly
        '''
        print 'X-shape: ' +repr(self.X.shape)
        print 'U-shape: ' +repr(self.U.shape)
        print 'Z-shape: ' +repr(self.Z.shape)
        print 'C-shape: ' +repr(self.C.shape)
        print 'V-shape: ' +repr(self.V.shape)
        print 'P-shape: ' +repr(self.P.shape)
        print 'R-shape: ' +repr(self.R.shape)
        print 'Xrec-shape: '+ repr(self.Xrec.shape)
        print 'meanZround: ' + repr(self.meanZround)
        print 'emptymeanZ: ' + repr(self.emptymeanZ)
        print 'check covariance matrix symmetry, C equals C.T : ' + repr(self.checkC)
        print 'Rows are eigenvectors if values are 1: ' + repr(self.checkVrows), '\n'
        print 'EigenVector, EigenValue pairs in descending order:'  
        print 'Note: Eigenvectors and values returned in order most to least importance', '\n', '\n'
        # Visually confirm that the list is correctly sorted by decreasing eigenvalues
        #for i in eig_pairs:
            #print(i[0])
    def dimensionReduction(self):
        # Make a list of (eigenvalue, eigenvector) tuples
        eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]

        # Sort the (eigenvalue, eigenvector) tuples from high to low
        eig_pairs.sort()
        eig_pairs.reverse()

        # Visually confirm that the list is correctly sorted by decreasing eigenvalues
        #print('Eigenvalues in descending order:')
        #for i in eig_pairs:
            #print(i[0])

    def visualizePCA(self):
        '''Method: Graph PCA results to show explained variance captured by n number of principle components
            requires plotly package
        '''
        tot = sum(self.Evals)
        var_exp = [(i / tot)*100 for i in sorted(self.Evals, reverse=True)]
        cum_var_exp = np.cumsum(var_exp)

        trace1 = Bar(
                x=['PC %s' %i for i in range(1,60)],
                y=var_exp,
                showlegend=False)

        trace2 = Scatter(
                x=['PC %s' %i for i in range(1,60)], 
                y=cum_var_exp,
                name='cumulative explained variance')

        data = Data([trace2])

        layout=Layout(
                yaxis=YAxis(title='Explained variance in percent'),
                title='Explained variance by different principal components')

        fig = Figure(data=data, layout=layout)
        py.iplot(fig)




class BayesClassifier: 
    '''A naive Bayesian classifier Object for n classes
        Classifier uses the PDF of Gaussian distribution in multiple dimensions

        initialization:
        myBayesClassifier = BayesClassifier()
        
        usage:
        ClassifyBayesD(TestDataset, TrainDataset, Classlabels, dimensions, allD)
        
        Usage arguments explained:

            TestDataset: Feature Vectors of unknown class membership
            TrainDataset: Feature Vectors of known class membership. Used to construct PDF
            and calculate probabilities

            ClassLabels: The labels for each class 

            Dimensions: number of feature vectors (dimensions) to use in training classifier (these are columns
            of the dataset)

            allD: boolean (True or False). If False (default) only the results of the specified dimension will be returned
            If true, all results up to specified dimensions will be returned by iterating through the classifier starting
            from 2 dimensions until specified dimensions (Recommended not to do this for very large number of dimensions)
            The print out of all dimensions lets the user see how classifier performs as dimensions are increased.
        
    '''
    def pdf(self, x,mu,sigma):
        '''Method: computes the probability density function for a Guassian distribution in multiple dimensions
        '''

        xf = x.astype(float) #make sure data is float
        muf = mu.astype(float) #mean vector

        d=np.alen(muf)#len first dimension of mean vector
        dfact1=(2*np.pi)**d
        dfact2=LA.det(sigma)
        fact=1/np.sqrt(dfact1*dfact2)
        xc=xf-muf
        isigma=LA.inv(sigma)
        #isigxc = np.dot(isigma,xc.T)
        #ex = np.dot(xc,isigxc)
        npdf = fact * np.exp(-0.5 * np.einsum('ij,jk,ik->i',xc,isigma,xc)) #note: getting einsum to work annoying
        return npdf  

    def BuildNDBayesianClassifier(self, Dataset, Classlabels, D):
        '''Constructor:  Initializes a ClassStat dictionary, with key class label to
         hold the distribution parameters for each class as values. 
        These values are stored in dictionaries as well so classStats is a dictionary of dictionaries.

        Imbedded dictionary for each class in the data is constructed containing key, values:
            Num: the total number of feature vectors in a class 
            Data: array of all the feature vectors belonging to the class
            Mean: mean vector of class' feature vectors 
            Cov: Covariance matrix of relating the feature vectors
            
            
        '''
        ClassStats = {} #dictionary to hold the distribution information for each class
        for n in range(len(Classlabels)):
            #initialize an entry in ClassStats dictionary for each class label
            ClassStats[Classlabels[n]]={}
            Class = Dataset[Dataset[:,-1] == Classlabels[n]]
            ClassData = (Class[:,:D]).astype(float)
            ClassStats[Classlabels[n]]['Num'] = len(Class)
            ClassStats[Classlabels[n]]['Data'] = ClassData
            ClassStats[Classlabels[n]]['Mean'] = np.mean(ClassData,axis=0)
            ClassStats[Classlabels[n]]['Cov'] = np.cov(ClassData, rowvar=False)
        #print ClassStats
        return ClassStats
    

    def ApplyNDBayesianClassifier(self, TestDataset, TrainDataset, Classlabels, D):
        '''Method: Uses ClassStat dictionary to get 'counts' of test data observations 
            by constructing a PDF of the training data and applying it on the test data.
            Helper function ResultLPBayesClassifier() to calculate probabilities and apply
            species label based on those probabilities. 
        '''
        ClassStats = self.BuildNDBayesianClassifier(TrainDataset, Classlabels, D)
        w=1; #width of the bin
        CountC_all = []
        for n in range(len(Classlabels)):
            NC = ClassStats[Classlabels[n]]['Num']
            UC = ClassStats[Classlabels[n]]['Mean']
            covC = ClassStats[Classlabels[n]]['Cov']
            testset = np.array((TestDataset[:,:D])).astype(float)
            countC = NC*w*self.pdf(testset, UC, covC) #gets number 
            #print countC
            CountC_all.append(countC)
        [resultlabel, resultprob]= self.ResultLPBayesClassifier(CountC_all, TestDataset, Classlabels)
        return np.array([resultlabel, resultprob])

    def ResultLPBayesClassifier(self, CountC_all, TestDataset, Classlabels):
        '''Method: Used by ApplyNDBayesianClassifier() to label results based on the highest probability 
            of value belonging to that class (the one with the max count)
        '''
        ClassCounts_all = np.array(CountC_all)
        #print ClassCounts_all.shape
        resultlabel = np.full(np.alen(TestDataset), "Indeterminate", dtype=object)
        resultprob = np.full(np.alen(TestDataset), 0 , dtype=float)
        for g in range(len(TestDataset)):
            CountXvalues = []
            for w in range(len(Classlabels)):
                count = ClassCounts_all[w][g]
                CountXvalues.append(count)
            max_value = max(CountXvalues)
            max_index = CountXvalues.index(max_value)
            label = Classlabels[max_index]
            resultlabel[g]=label
            #print sum(ClassCounts_all)
            resultprob = (ClassCounts_all[max_index][g]).astype('float')/sum(ClassCounts_all)
        return resultlabel, resultprob

    def ClassifyBayesD(self, TestDataset, TrainDataset, Classlabels, dimensions, allD):
        '''Method: This function runs the Bayes Classifier in D dimensions. This is used for
             running a classifier using PCA results, D representing number of principle components.
        '''
        DBayesPCAResults=[]
        #if user wants all PCA dimension results printed up to specified dimension(s)
        if allD == True: 
            for d in range(dimensions-1):
                D=d+2
                nBResults = self.ApplyNDBayesianClassifier(TestDataset, TrainDataset, Classlabels, D)
                DBayesPCAResults.append(nBResults)
            return np.array(DBayesPCAResults)
        #if user wants the PCA dimension results for only the specified dimension(s)
        else:
            D=dimensions
            nBResults = self.ApplyNDBayesianClassifier(TestDataset, TrainDataset, Classlabels, D)
            return np.array(nBResults)
        
    
class ClassifierPerformanceMetrics:
    '''ClassifierPerformanceMetrics gives the user a summery breakdown of how the classifier performed
        It initializes a ClassiferObject to perform classification within certain parameters
        The program prints the summary to std out.
        
        instantiation:  
            myClassifier = ClassifierPerformanceMetrics()
        
        Usage: 
            ResultArrayBayesPCAperformance = myClassifier.BayesPCAperformance(testSetBayes, trainSetBayes, classlabels, dimensions, allD)

        Usage arguments explained: (Note: same args as for BayesClassifier, 
                                    designed this way for future application
                                    to other classifiers)

            TestDataset: Feature Vectors of unknown class membership
            TrainDataset: Feature Vectors of known class membership. Used to construct PDF
            and calculate probabilities

            ClassLabels: The labels for each class 

            Dimensions: number of feature vectors (dimensions) to use in training classifier (these are columns
            of the dataset)

            allD: boolean (True or False). If False (default) only the results of the specified dimension will be returned
            If true, all results up to specified dimensions will be returned by iterating through the classifier starting
            from 2 dimensions until specified dimensions (Recommended not to do this for very large number of dimensions)
            The print out of all dimensions lets the user see how classifier performs as dimensions are increased.
    '''
    
    def __init__(self):
        '''Constructor: initialize BayesClassifier object
        '''
        self.BayesClassifier = BayesClassifier()

    def PerformanceMetrics(self, Resultlabels, Dataset, PositiveLabel, Classlabels):
        '''Method: Calculates performance metrics for Bayes Classifier. Each of the classlabels are
            is evaluated and given performance metrics. Currently the classifier is evaluated by Accuracy, Sensitivity, Specificity, and 
            Positive Predictive Value. Function also prints metrics in readable format.
        '''
        OutputCL = (Resultlabels).astype('str')#output labels
        #print OutputCL.shape
        GroundTruth = (Dataset[:,-1]).astype('str') #true labels
        Classcomps = OutputCL == GroundTruth #where the output label is the true label
        #print Classcomps
        
        TrueP=0
        FalseP=0
        TrueN=0
        FalseN=0
        Num=0

        for i in range(len(GroundTruth)):
            if Classcomps[i]== True:
                if OutputCL[i] == PositiveLabel:
                    TrueP+=1
                else:
                    TrueN+=1
            elif Classcomps[i]==False:
                if OutputCL[i] != PositiveLabel:
                    FalseN+=1
                else:
                    FalseP+=1
    
        Accuracy = float((TrueP+TrueN))/(TrueP+TrueN+FalseP+FalseN)
        Sensitivity = float((TrueP))/(TrueP+FalseN)
        Specificity = float((TrueN))/(FalseP+TrueN)
        #Calculate Positive Predictive Values for each class
        if (TrueP+FalseP)==0: #handles cases where no instances of a class label are positive
            PPV = 0
        else:
            PPV = float((TrueP))/(TrueP+FalseP)
        stringmeasures = ['TrueP', 'FalseP', 'TrueN', 'FalseN', 'Accuracy','Sensitivity', 'Specificity','PPV']
        measures_values = [TrueP, FalseP, TrueN, FalseN, Accuracy, Sensitivity, Specificity, PPV]
        print 'Classifier Performance:'
        print '     Positive Class Label: '+ repr(PositiveLabel)
        for i in range(len(stringmeasures)):
            print '         '+stringmeasures[i]+ ': '+repr(measures_values[i])
    
        return [TrueP, FalseP, TrueN, FalseN, Accuracy, Sensitivity, Specificity, PPV]


    def BayesPCAperformance(self, TestDataset, TrainDataset, Classlabels, dimensions, allD):
        '''Method: Function mediates between Bayes Classifier functions and performance metrics.
            this function is currently the one called by main() to perform classification. This 
            function passes the datasets to Bayes Classifier functions and then passes those results to 
            performance metrics for evaluation and printing output to stdout.
        '''
        BayesPCAperformance =[]
        DResults = self.BayesClassifier.ClassifyBayesD(TestDataset, TrainDataset, Classlabels, dimensions, allD)
        #DResults is an array of result labels and result probability of testset data
        PC=len(DResults)
        if allD==False:
            PC = 1
        for i in range(PC):
            if allD==False:
                print repr(dimensions)+' Principal Components '+ repr(len(Classlabels))+'-Class Bayes Classifier Performance:'
            else:
                print repr(i+2)+' Principal Components Bayes Classifier Performance:'
            labelstring = []
            reallabelstring = []
            #print DResults[0][0]
            OutputCL = ((DResults[0])).astype('str') #testset result labels
            realCL=(TrainDataset[:,-1]).astype('str')#training set result labels
            #iterates through results and formats output for printing
            for j in range(len(Classlabels)):
                labelnum = OutputCL[OutputCL == Classlabels[j]]#number of results from classlabel[j]
                lstring = repr(Classlabels[j])+ ': '+ repr(len(labelnum))
                labelstring.append(lstring)
                reallabelnum = realCL[realCL == Classlabels[j]]
                reallstring = repr(Classlabels[j])+ ': '+ repr(len(reallabelnum))
                reallabelstring.append(reallstring)
            print '  Training Output ->' + repr(reallabelstring) + ' Total: ' +repr(len(realCL))
            print '  Testing Output -->' + repr(labelstring)+ ' Total: ' +repr(len(OutputCL))
            totalPPV = []
            for cl in range(len(Classlabels)):
                PositiveLabel=Classlabels[cl]
                if allD==True:
                    nDBayesPerformance = self.PerformanceMetrics(DResults[i][0], TestDataset, PositiveLabel, Classlabels)
                    totalPPV.append(nDBayesPerformance[-1])
                else:
                    nDBayesPerformance = self.PerformanceMetrics(DResults[0], TestDataset, PositiveLabel, Classlabels)
                    totalPPV.append(nDBayesPerformance[-1])
                BayesPCAperformance.append(nDBayesPerformance)
        
            print 'Avg PPV: '+ repr(np.mean(totalPPV))   
            print '\n'
        return np.array(BayesPCAperformance)
        

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

    Usage:
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)
    print (myCommandLine.args)  
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog -User can provide defined options for bayesClassifier.py program to specify % of data to split into train and testing subsets, setting random seed, number of dimensions, and whether or not to print all dimensions results. See usage for how to set options on commandline', 
                                             epilog = 'input filefile must be in .xlsx specified on commandline according to usage syntax', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s inFileName [options] -option1[default]'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('-d', '--dimensions', type=int, default=3, action = 'store', help='specify the number of dimensions (PCA components) to use in Bayes Classifier, default 3')
        self.parser.add_argument('-pC', '--printPCAcheck', action = 'store', nargs='?', const=True, default=False, help='if True prints out XUZCVPR summary of PCA')        
        #will add these later:
        #self.parser.add_argument('-p', '--printPCA', action = 'store', nargs='?', const=True, default=False, help='if True prints out PCA X and Xrec results all possible principle components')
        #self.parser.add_argument('outFile', action = 'store', default='output', nargs='?', help='output file name') 

        self.parser.add_argument('-aD', '--allDimensions', action = 'store', nargs='?', const=True, default=False, help='if False(default) prints only the Bayes Classifier results of specified # of dimensions, if True prints Bayes Classifier result upto and including specified # of dimensions (ie just 5 or 1, 2, 3, 4, and 5')
        self.parser.add_argument('-rS', '--randomSeed', type=int, default=0, action = 'store', help='specify a randomSeed value for reproducibility of results')
        self.parser.add_argument('-dS', '--dataSplit', type=float, default=0.75, action = 'store', help='specified percentage of data (as float) to be randomly split as training set, the remaining will be used as test (ex: .80 for 80%)')
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

    rwExcel = readwriteExcel
    dataset = np.array(rwExcel.readExcel(myCommandLine.args.inFile))
    
    myData = ClassData()
    randomSeed = myCommandLine.args.randomSeed
    split = myCommandLine.args.dataSplit
    [allDataset, TrainDataset, TestDataset]= myData.TestTrainDataSplit(dataset, split, randomSeed)
    [classlabels, trainlabels, testlabels] = myData.dataLabels()

    #Run PCA on training data
    trainVectorsPCA = TrainDataset[:,:-1] #only the feature vectors, remove class labels
    myPCA = PCA()#initialize PCA object
    [X, U, Z, C, V, EighVals, P, R, Xrec] = myPCA.XUZCVPR(trainVectorsPCA, Utrain=1, train=True)
    if myCommandLine.args.printPCAcheck:
        myPCA.PCAcheck()

    #plotPCA results
    #PCA.visualizePCA()
    
    
    #Prep classifier testdata based on PCA results
    testZ = (TestDataset[:,:-1])-U
    testP = np.dot(testZ,V.T)
    
    trainSetBayes = np.concatenate((P, trainlabels.T), axis=1)
    testSetBayes = np.concatenate((testP, testlabels.T), axis=1)

    dimensions = myCommandLine.args.dimensions
    allD = myCommandLine.args.allDimensions
    myClassifier = ClassifierPerformanceMetrics()
    ResultArrayBayesPCAperformance = myClassifier.BayesPCAperformance(testSetBayes, trainSetBayes, classlabels, dimensions, allD)
    
    #uncomment to write results to file 

    #filename = 'Results'+repr(myCommandLine.args.inFile)
    #file = open(filename, "w")   
    #file.write(ResultArrayBayesPCAperformance)
    #file.close

   
if __name__ == "__main__":
        main() 