# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:23:42 2020

@author: arn4va
"""
#import nltk
#nltk.download('punkt')

import sys
import numpy as np
import pandas as pd
from urllib.request import Request, urlopen
#from nltk.tokenize import sent_tokenize, word_tokenize
import re

def striphtml(data):
    p = re.compile(r'<.*?>')
    return re.split(p,data)

#read in csv with drug IDs
#df=pd.read_csv('drugID_names_matched.csv')
df=pd.read_csv(sys.argv[1])

print(df)

def webScrapeDrugAction(df):
    #get drug IDs
    IDs=df.Drug.values
    
    actions=list()
    
    k=1
    url_test='https://www.drugbank.ca/drugs/DB12010'
    
    for y in range(0,len(df)):
        ID=IDs[y].strip()
        url=('https://www.drugbank.ca/drugs/'+ID)
        print(url)
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        webUrl=urlopen(req)
        page=webUrl.read().decode('utf-8') #decode converts bytes to string object
        
        stripped=striphtml(page)
        
        while('' in stripped) : 
            stripped.remove('') 
            
        #Get the name of the drug target
        target=df.Target[y]
        
        #Find the drug target dropdown box text
        
        indices = [i for i, x in enumerate(stripped) if x == target]
        
        #The word 'kind' comes after the target name in the targets box on the drugbank
        #page, so determine is this word follows the target name
        
        keep=list()
        
        for k in indices:
            print(stripped[k+1])
            if stripped[k+1]=='Kind':
                keep.append(k)
        if(len(keep)!=1):#insure at least and only 1 index was saved to keep
            print('WARNING: Drug ID '+ID+' Action Was Not Correctly Assigned, Inspect Webpage!' )
            actions.append('NoMatch')
        
        #establish search range following correct target name index
        if(len(keep)==1):
            
            ant_bool=False
            ag_bool=False
            for i in range (keep[0],keep[0]+11):
                word=stripped[i]
                
                if (word.lower() == 'inhibitor') or (word.lower()=='antagonist') or (word.lower()=='antagonists' or (word.lower() == 'inhibitors')):
                    ant_bool=True
                    break
                    
                else:
                    ant_bool=False
                    
            for i in range (keep[0],keep[0]+11):   
                word=stripped[i]
                
                if (word.lower() == 'agonist') or (word.lower() == 'agonists') or (word.lower() == 'inducer') : #add 'inducer'
                    ag_bool=True
                    break
                else:
                    ag_bool=False
                    
            if (ant_bool & ~ag_bool):
                actions.append('Antagonist')
                
            elif (~ant_bool & ag_bool):
                actions.append('Agonist')
            
            elif (~ant_bool & ~ag_bool):
                actions.append('None')
            
            else :
                actions.append('Conflicting')    
                
    print(actions)
    df['Actions_from_Python']=actions
    df.to_csv('DrugActionsTest_Output.csv')
    
    exitStatus='Finished!'
    return exitStatus, actions

webScrapeDrugAction(df)
print('Done!')
                
                
                
                
                
            