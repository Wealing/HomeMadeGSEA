#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:02:00 2017

@author: Daniel F Hernández: dan.fhg@gmail.com
"""
import csv, scipy.stats, numpy

def hypergeoPval(N,K,n,k):
    
    '''
    N is the population size,
    K is the number of success states in the population,
    n is the number of draws,
    k is the number of observed successes,
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
    '''
    
    #check if the inputs are correct
    if ~all(isinstance(x,int) for x in [N,K,n,k]):
        ValueError
    
    #print( 'p-value P(X= ' + str(k) + '): ' + str(scipy.stats.hypergeom.pmf(k,N,K,n)))
    #print( 'p-value P(X>= ' + str(k) + '): ' + str(scipy.stats.hypergeom.sf(k - 1,N,K,n)))
    
    return scipy.stats.hypergeom.sf(k - 1,N,K,n)
    
def firstAssignment(filename,enrichedGenes,outfile_name = 'firstAssignment.tsv'):
    
    # load the enriched genes
    with open(enrichedGenes, mode='r') as infile:
        reader = csv.reader(infile)
        enrichedGenes=next(reader)
        enrichedGenes = [x.strip(' ') for x in enrichedGenes]
                    
    with open(filename, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')

        # get the total number of genes (needed to calculate p-vavlue)
        allGenes = set()
        pathwayName, pathwayDesc = [],[]
        for row in reader:
            pathwayName.append(row[0])
            pathwayDesc.append(row[1])
            allGenes.update(row[2:]) #MUST MAKE SURE SAME GENES ARE WRITTEN SAME WAY
        allGenesNum = len(allGenes) 
        
        # calculate p value
        infile.seek(0)# restart beginning
        p,enrichedInPathway,pathwayLen,listEnrichedInPathway=[],[],[],[]
        for row in reader:

            # save number of enriched genes in metabolic pathways
            listEnrichedInPathway.append(set(enrichedGenes).intersection(row[2:]))
            enrichedInPathway.append(len(listEnrichedInPathway[-1] ))
            pathwayLen.append(len(row[2:]))
            p.append(hypergeoPval(allGenesNum,len(enrichedGenes),pathwayLen[-1],enrichedInPathway[-1]))
            
        
    # output pathways ordered by p-value
    vals = numpy.array(p)
    sort_index = numpy.argsort(vals)
    
    with open(outfile_name, "w") as text_file:
        text_file.write('index\tEnrichment\tp-value\tgene set name\tdifferentially expressed genes in gene set\n')
        for num,i in enumerate(sort_index):

            if num <     len(p)-sum(vals==1):
                text_file.write(str(i+1) + '\t' + str(enrichedInPathway[i])+' of '+str(pathwayLen[i]) + '\t' + str(p[i]) + '\t' + str(pathwayName[i]) +'\t'+ str(listEnrichedInPathway[i]) + '\n')
            else:
                #text_file.write(str(i+1) + '\t' + str(enrichedInPathway[i])+'/'+str(pathwayLen[i]) + '\t' + str(p[i]) + '\t' + str(pathwayName[i]) + '\n')
                text_file.write('no other gensets had differentially expressed genes')
                break
    return vals
        
firstAssignment('pathways.txt','enrichedGenes.csv')
