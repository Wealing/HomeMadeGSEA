#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:02:00 2017

@author: Daniel Hernandez

Homemade GSEA
"""
from __future__ import division
import csv, scipy, numpy, math, matplotlib

def rankgenes(my_data,group1,group2,GeneNames):
    # rank genes according to differential expression
    # uses simple wilcoxon test
    stat,pvalue=[],[],
    for row in my_data[1:,]:#ASSUMING FIRST ROW HEADER
        stat_tmp,pvalue_tmp=scipy.stats.ranksums(row[group1],row[group2])
        stat.append(stat_tmp)
        pvalue.append(pvalue_tmp)
    #order is the index vector of the genes
    order=sorted(range(len(stat)), key=lambda k: stat[k])
    #order = numpy.argsort(stat)#THERE ARE TRANSCRIPTS WITH EXACTLY SAME DISTRIBUTIONS

    orderedGenes = [ GeneNames[i] for i in order]
    return orderedGenes
    
def readArrayData(arrayData):
       #Input Array data
    with open(arrayData, mode='r') as infile:
        
        # get phenotype groups
        reader = csv.reader(infile, delimiter='\t')
        groups=next(reader)#ASSUMING FIRST ROW HEADER
        uniqueGroups=list(set(groups[1:]))#ASSUMING FIRST COLUMN LABEL
        
        #get names of Genes
        GeneNames = []
        for row in reader:
            GeneNames.append(row[0])
        print(arrayData)
    return groups, uniqueGroups,GeneNames
    
def calcES(pathwaysFilename,orderedGenes):    
    # Calculate Enrichment score
    with open(pathwaysFilename, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')

        # 
        pathwayName, pathwayDesc,measuredFromPathway =[], [],[]
        N=len(orderedGenes)
        runningSumPerPathwayAll,ESperPathAll=[],[]
        trailAll=[]
        for pathway in reader:
            trail = []
            
            pathwayName.append(pathway[0])
            pathwayDesc.append(pathway[1])
            
            #Ns=len(pathway[2:])
            Ns=len(set(orderedGenes).intersection(pathway[2:])) 
            measuredFromPathway.append(set(orderedGenes).intersection(pathway[2:]))
            runningSumPerPathway,ESperPath=0,0
            
            # running sum max value is the Enrichment Score
            if Ns==0: ESperPath, runningSumPerPathway=None,None
            else:
                for gene in orderedGenes[0:]:
                    #print(gene)
                    
                    
                    if gene in pathway[2:]:
                        runningSumPerPathway += math.sqrt( (N-Ns)/Ns )
                    else:
                        runningSumPerPathway += - math.sqrt( Ns/(N-Ns) )
                        
                    ESperPath=max(ESperPath,math.fabs(runningSumPerPathway))
                    trail.append(runningSumPerPathway)
            runningSumPerPathwayAll.append(runningSumPerPathway)
            trailAll.append(trail)
            ESperPathAll.append(ESperPath)
#            plt.plot(trail)
#            print(trail)
            
        return ESperPathAll, trailAll, pathwayName, pathwayDesc,measuredFromPathway
            
    
    
    
    
def secondAssignment(pathwaysFilename, arrayData,randomPermNum=50, outfile_name ='secondAssignment.tsv'):
    
    groups, uniqueGroups,GeneNames = readArrayData(arrayData)
    
    ###### run first time #######
    # first group indeces
    group1 = numpy.where(numpy.array(groups)==uniqueGroups[0])[0]
    # second group indeces
    group2 = numpy.where(numpy.array(groups)==uniqueGroups[1])[0]
    #save all data to numpy array #THIS WILL NOT SCALE
    my_data = numpy.genfromtxt(arrayData, delimiter='\t')
    orderedGenes = rankgenes(my_data,group1,group2,GeneNames)
    #takes 3 seconds, slowest step
    ESperPath, trailAll, pathwayName, pathwayDesc,measuredFromPathway=calcES(pathwaysFilename,orderedGenes)
    
    
    
    ###### run GSEA with 50 random permutations to generate null distribution #######
    ESperPathNull=[]
    for i in range(0,randomPermNum):# number of random permutations, the more, the more robust the findings
        permutedCategories = list(numpy.random.permutation(groups[1:]))
        # first group indeces
        group1 = numpy.where(numpy.array([groups[0]]+permutedCategories)==uniqueGroups[0])[0]
        # second group indeces
        group2 = numpy.where(numpy.array([groups[0]]+permutedCategories)==uniqueGroups[1])[0]
        orderedGenes = rankgenes(my_data,group1,group2,GeneNames)
        
        ESperPathNullTemp, _,_,_,_ = calcES(pathwaysFilename,orderedGenes)
        ESperPathNull.append(ESperPathNullTemp)
        print(str(i+1)+' of '+ str(randomPermNum) +' random permutations being added to null distribution' )


    
    ##############
    #compare single real ES with random distribution to determine over expression
    nullDistributionAll,statisticAll, pvalueAll = [],[],[]
    for i in range(len(pathwayName)):
        
        nullDistribution = []
        for j in range(len(ESperPathNull)):
            nullDistribution.append(ESperPathNull[j][i])
        nullDistributionAll.append(nullDistribution)
        
        if ESperPath[i] is not None:
            statistic, pvalue = scipy.stats.ttest_1samp(nullDistribution,ESperPath[i])
        else: statistic, pvalue = 0,1 # THESE WILL NEED TO BE REMOVED
        statisticAll.append(statistic)
        pvalueAll.append(pvalue)
    order = numpy.argsort(statisticAll)
    ESnormalizedPerPath2 = 1-(numpy.array(statisticAll)-min(statisticAll))/(max(statisticAll)-min(statisticAll))
    
    
    
    ####################3
    # output pathways ordered by p-value, actually sorting on the statistic though.
    print('\n\n text file '+ outfile_name+' has been created with results\n\n' )
    
    matplotlib.pyplot.plot(trailAll[order[0]])
#    matplotlib.pyplot.plot(trailAll[order[250]])
#    matplotlib.pyplot.plot(trailAll[order[500]])

    with open(outfile_name, "w") as text_file:
        text_file.write('index\tEnrichment\tp-value\tgene set name\tdifferentially expressed genes in gene set\n')
        for i in order:
            text_file.write(str(i+1) +'\t'+str(ESnormalizedPerPath2[i]) + '\t' + str(pvalueAll[i])+ '\t'+str(pathwayName[i])+  '\n')
           
        
groups= secondAssignment('pathways.txt','leukemia.txt',randomPermNum=50)