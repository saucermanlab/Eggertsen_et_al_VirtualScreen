import pandas as pd
#import os
#os.chdir("C:/Users/HP/Documents/RemoteSauce/DRUG_CODE")


import sys

print ('Number of arguments:'+ str(len(sys.argv))+ 'arguments.')
print ('Argument List:'+ str(sys.argv))

# read in matlab ODE file as text file

targetCSV=sys.argv[1]
actionCSV=sys.argv[2]
model=sys.argv[3]
idCSV=sys.argv[4]  

###actionCSV='"DrugActionsTest_Output.csv"';
###targetCSV='"allPharm.csv"';
###idCSV='"drugIDs.csv"';
###model='Cardiac_Hypertrophy_Model.xlsx';

# remove excess quotations [does this need it?]
targetCSV = targetCSV.strip('\"')
actionCSV = actionCSV.strip('\"')
model = model.strip('\"')
idCSV = idCSV.strip('\"')

def ID_pull(targetCSV, actionCSV, model, idCSV):

    # import drug target data
    targets = pd.read_csv(targetCSV)
    drugs = targets['Drug IDs']  # DrugBank ID
    geneName = targets['Gene Name'] 
    drugName = pd.read_csv(idCSV)
    ###drugClass = pd.read_excel("All Drugs Comp NonComp Classification.xlsx")
    drugAction = pd.read_csv(actionCSV)
    drugAction['Action']=drugAction['Actions_from_Python']
    drugAction_strip = pd.Series([0]*len(drugAction))
    for i in range(0,len(drugAction)):
        drugAction_strip[i] = drugAction['Drug'][i].lstrip()
    drugAction['Drug'] = drugAction_strip
    
    # import model
    nodes = pd.read_excel(model, skiprows=[0]);
    rnodeID = nodes['ID']
    rnodeGene = nodes['gene name']; 
    
    nodeGene = []; nodeID = []; 
    for i in range(0,len(rnodeID)):
        if (isinstance(rnodeGene[i],float)) == False:
            nodeID.append(rnodeID[i])
            nodeGene.append(rnodeGene[i])
    
    nodeGeneIDs = [None]*len(nodeGene); geneIDs = []; nodeIDs = [];
    # separate multiple genes and assign to each node for given model
    for i in range(0,len(nodeGene)):
        nodeGeneIDs[i] = nodeGene[i].replace(" ", "").split(";")
        for j in range(0,len(nodeGeneIDs[i])):
            geneIDs.append(nodeGeneIDs[i][j])
            nodeIDs.append(nodeID[i])
    
    # pull drugs associated with gene target
    dlist = []; glist = [];
    for i in geneIDs:
        if (i in list(geneName)) == True:
            dlist.append(drugs[list(geneName).index(i)])
            glist.append(i)
    drugs = dlist; geneName = glist;
    
    #organize to single gene and multiple effected nodes
    nodeEffect = [None]*len(geneName);
    for i in range(0,len(geneName)):
        ngroup = [];
        for j in range(0,len(geneIDs)):
            if geneName[i] == geneIDs[j]:
                ngroup.append(nodeIDs[j])
                nodeEffect[i] = ngroup;
    
    # separate multiple drugs for single gene
    drugIDs = [None]*len(drugs)
    for i in range(0,len(drugs)):
        drugIDs[i] = drugs[i].split("; ")
    
    nname = []; gname = []; dname = []; dclass = []; daction = []; did = [];
    for i in range(0,len(drugIDs)): 
        tset = drugIDs[i]
        for j in range(0,len(tset)):
            # match DrugBank ID to drug name for each associated gene target
            drug_id = drugIDs[i][j]
            if drugName['DrugBank ID'].str.contains(drug_id).sum() > 0:
                drug_name = drugName[drugName['DrugBank ID'] == drug_id]['Name'].iloc[0]
                gname.append(geneName[i])
                nname.append(nodeEffect[i])
                did.append(drug_id)
                dname.append(drug_name)
                # list action (agonist or antagonist) for each drug
                if (drug_id in list(drugAction['Drug'])) == True:
                    drug_action = drugAction[drugAction['Drug'] == drug_id]['Action'].iloc[0]
                    if (drug_action.lower() == 'agonist') or (drug_action.lower() == 'inducer'):
                        drug_action = 'Agonist'
                    elif (drug_action.lower() == 'antagonist') or (drug_action.lower() == 'inhibitor'):
                        drug_action = 'Antagonist'
                elif (drug_name in list(drugAction['Drug'])) == False:
                    drug_action = 'N/A'
                daction.append(drug_action)
    
    dclass = [None]*len(dname) # DELETEME after including this part of the code...
    
    # place Gene Name, Drug Name and Classification into the
    # Drug Target Matrix (DTM)
    dt = {'Drug ID': did,'Drug Name': dname,'Action': daction,'Gene Name': gname,'Node Name': nname,'Class': dclass}
    DTM = pd.DataFrame(dt, columns = ['Drug ID','Drug Name','Action','Gene Name','Node Name','Class'])
    
    # populate DrugsToSimulate matrix
    isagonist = []; agonisttarget = []; agonisttargetgeneid = []; agonisttargetindex = [];
    isantagonist = []; antagonisttarget = []; antagonisttargetgeneid = []; antagonisttargetindex = [];
    simdrugs = []; 
    for i in range(0,len(DTM)):
        if daction[i] == 'Agonist':
            isagonist.append('Yes'); isantagonist.append('No');
            agonisttarget.append(tuple(nname[i])); agonisttargetgeneid.append(gname[i]);
            if len(nname[i])<2:
                agonisttargetindex.append(str(list(nodes['ID']).index(nname[i][0])+1))
            else:
                nhold = [];
                for k in range(0,len(nname[i])):
                    nhold.append(str(list(nodes['ID']).index(nname[i][k])+1))
                agonisttargetindex.append(';'.join(nhold))
            antagonisttarget.append(''); antagonisttargetgeneid.append(''); antagonisttargetindex.append('');
        elif daction[i] == 'Antagonist':
            isagonist.append('No'); isantagonist.append('Yes');
            antagonisttarget.append(tuple(nname[i])); antagonisttargetgeneid.append(gname[i]);
            if len(nname[i])<2:
                antagonisttargetindex.append(str(list(nodes['ID']).index(nname[i][0])+1))
            else:
                nhold = [];
                for k in range(0,len(nname[i])):
                    nhold.append(str(list(nodes['ID']).index(nname[i][k])+1))
                antagonisttargetindex.append(';'.join(nhold))
            agonisttarget.append(''); agonisttargetgeneid.append(''); agonisttargetindex.append('');
        else:
            isagonist.append(''); isantagonist.append(''); agonisttarget.append(''); agonisttargetgeneid.append(''); 
            agonisttargetindex.append(''); antagonisttarget.append(''); antagonisttargetgeneid.append('');
            antagonisttargetindex.append(''); 
        # similar drugs algorithm
        simset = []; checkset = [DTM['Drug Name'][i]];
        for j in range(0,len(DTM)):
            if DTM['Drug Name'][j] != DTM['Drug Name'][i]:
                if DTM['Gene Name'][j] == DTM['Gene Name'][i]:
                    if DTM['Action'][j] == DTM['Action'][i]:
                        simset.append(DTM['Drug Name'][j])
                        checkset.append(DTM['Drug Name'][j])
        simdrugs.append(';'.join(simset))
        
        
    ds = {'Drug':tuple(dname),'IsAgonist':tuple(isagonist),'AgonistTargetGeneID':tuple(agonisttargetgeneid),
          'AgonistTarget':tuple(agonisttarget),'AgonistTargetIndex':tuple(agonisttargetindex),'IsAntagonist': tuple(isantagonist),
          'AntagonistTargetGeneID':tuple(antagonisttargetgeneid),'AntagonistTarget':tuple(antagonisttarget),'AntagonistTargetIndex':tuple(antagonisttargetindex),
          'DrugAction':tuple(['Competitive']*len(DTM)),'SimilarDrugs':tuple(simdrugs)}
    DSM = pd.DataFrame(ds)
    DSM = DSM.drop_duplicates(); DSM = DSM.reset_index(drop=True);
    
    #drop empty rows
    dropdex = [];
    for i in range(0,len(DSM)):
        if DSM['IsAgonist'][i] == '':
            dropdex.append(i);
        DSM['AgonistTarget'][i]=";".join(DSM['AgonistTarget'][i])
        DSM['AntagonistTarget'][i]=";".join(DSM['AntagonistTarget'][i])
    DSM = DSM.drop(dropdex).reset_index(drop=True);
    DSM = DSM.sort_values(by='Drug').reset_index(drop=True)
    foo = lambda a: ";".join(a)
    DSM = DSM.groupby(by='Drug').agg({'Drug':'first',
                                      'IsAgonist':'first',
                                      'AgonistTargetGeneID':foo,
                                      'AgonistTarget':foo,
                                      'AgonistTargetIndex':foo,
                                      'IsAntagonist':'first',
                                      'AntagonistTargetGeneID':foo,
                                      'AntagonistTarget':foo,
                                      'AntagonistTargetIndex':foo,
                                      'DrugAction':'first',
                                      'SimilarDrugs':'first'})  ###check this line!!!!!###
    
    for i in range(0,len(DSM)):
        DSM['AgonistTargetIndex'][i]=";".join(list(set(DSM['AgonistTargetIndex'][i].split(';'))))
        DSM['AntagonistTargetIndex'][i]=";".join(list(set(DSM['AntagonistTargetIndex'][i].split(';'))))
        DSM['AgonistTarget'][i]=";".join(list(set(DSM['AgonistTarget'][i].split(';'))))
        DSM['AntagonistTarget'][i]=";".join(list(set(DSM['AntagonistTarget'][i].split(';'))))
        
    # write to CSV file
    DSM.to_csv('DrugsToSimulate.csv',index = False,header = True)
    
ID_pull(targetCSV,actionCSV,model,idCSV)