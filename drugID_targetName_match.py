# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 15:15:50 2020

@author: arn4va
"""


import pandas as pd
import sys

print ('Number of arguments:'+ str(len(sys.argv))+ 'arguments.')
print ('Argument List:'+ str(sys.argv))

dataCSV=sys.argv[1]
modelEXCEL=sys.argv[2]

def IDmatch(dataCSV,modelEXCEL):
##    
#    dataCSV='allPharm.csv'
#    modelEXCEL='fib617_references.xlsx'
    data=pd.read_csv(dataCSV)
    model= pd.read_excel(modelEXCEL, skiprows=[0])
    
    drugIDs=list()
    geneSymbols=list()
    geneNames=list()
    
    model_genes=model['gene name'].dropna() #remove NaN's
    #Get list of genes in model
    geneList=list()
    for gene in model_genes:
        genes=gene.strip().split('; ')
        for g in genes:
            geneList.append(g)
            
    #match gene symbols to gene names and drug IDs
    for i in geneList:
        pharm_select = data[data['Gene Name'].str.lower()==i.lower()]
        if not pharm_select.empty:
            IDs=pharm_select['Drug IDs'].values[0].split(';')
            for d in IDs:
                drugIDs.append(d)
                geneSymbols.append(i)
                geneNames.append(pharm_select['Name'].values[0])
                
    
    matchDF=pd.DataFrame({'Drug':drugIDs,'Target':geneNames})
    
    matchDF.to_csv('drugID_names_matched.csv',index = False)
    

IDmatch(dataCSV,modelEXCEL)