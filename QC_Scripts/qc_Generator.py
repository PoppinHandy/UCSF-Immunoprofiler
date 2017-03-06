# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:04:40 2017

@author: andyp
"""

import pandas as pd
import os
import math
import numpy as np

# Query the db for specific rows and cols that the user inputs and report any percentage between stains greater than 10% for QC purposes
def make_QC_Sheet(filename, df):
    """Takes in a file to be written and generates a standard QC file sheet."""
    row = ["CD45+", "CD3+ all", "CD4+", "Q3: CD8a+ , CD4-", "NK", "B-cells", "Lymphocytes+", "HLADR+, Lymphocytes-", "Calculations"]
    dfDict = dict()
    dfDict["%NK/CD45"] = list()
    cd = dict() # Calculations dicitonary for later
    
    # For each name in the list
    for r in row:
        # NK Cells
        if (r == "NK"):
            nk_Stain2 = df.loc[:, "Cell Type"] == "HLADR-, CD3-, CD56+ (NK)"
            df_nk2_Query = df[nk_Stain2]
            dfDict["%NK/CD45"].append(df_nk2_Query)

            nk_Stain4 = df.loc[:, "Cell Type"] == "CD16+ NK cells"
            df_nk4_Query = df[nk_Stain4]
            dfDict["%NK/CD45"].append(df_nk4_Query)
                         
        # B-cells
        elif (r == "B-cells"):
            nk_B = df.loc[:, "Cell Type"] == "B-cells"
            df_B_Query = df[nk_B]
            dfDict["%B/CD45"] = df_B_Query

        
        elif (r == "CD4+"):
            cd3 = df.loc[:, "Parent"] == "CD3+ all"
            cd4 = df["Cell Type"].str.contains("CD4+", regex=False)
            df_Query = df[cd3 & cd4]
            dfDict["%CD4/CD3"] = calculateCD8AndCD4(df_Query)
            
        elif ("CD8a+" in r):
            cd3 = df.loc[:, "Parent"] == "CD3+ all"
            cd8 = df.loc[:, "Cell Type"].str.contains("CD8a+", regex=False)
            df_Query = df[cd3 & cd8]
            dfDict["%CD8/CD3"] = calculateCD8AndCD4(df_Query)
                   
        elif(r == "Calculations"):
            stain2 = df.loc[:, "Stain Name"] == "Stain 2"
            df_Calculations = df[stain2]
            cd = calculations(df_Calculations)
                              
        elif ("Lymphocytes+" in r):
            parent = df.loc[:, "Parent"] == "CD45+"
            hasLymph = df.loc[:, "Cell Type"] == r
            df_Query = df[hasLymph & parent]
            dfDict["%Lymphocytes/CD45"] = df_Query

        else:
            hasCell = df.loc[:, "Cell Type"] == r
            df_Query = df[hasCell]
            dfDict[r] = df_Query
            
    df_Samp = set(df["Sample Name"].values)                
    sampleSet = set()
    # Building a set of names for each individual sample
    for sn in df_Samp:
        name = sn.split("_")[0] + "_" + sn.split("_")[1]
        sampleSet.add(name)
        
    sampleList = sorted(sampleSet)
    for s in sampleList:
        indFile = "QC/QC_" + s + ".tsv"
        iF = open(indFile, 'w')
        p = s + "\n"
        filename.write(p)
        iF.write(p)
        dfQC = formatQCSheet(s, dfDict)
        dfQC_formatted = pd.pivot_table(dfQC, index="Cell Type", columns="Stain Name", values="Percentage")
        dfQC_formatted = dfQC_formatted.reindex(["CD45+", "CD3+ all", "CD4", "CD8", "%NK/CD45", "%B/CD45", "%Lymphocytes/CD45", "HLADR+, Lymphocytes-"])
        for calcs in cd[0].keys():
            if calcs == s:
                dfQC_formatted.loc["%Lymphocytes/CD45", "Stain 2"] = cd[0][calcs]
        for calcs in cd[1].keys():
            if calcs == s:
                dfQC_formatted.loc["HLADR+, Lymphocytes-", "Stain 2"] = cd[1][calcs]
        dfQC_formatted.to_csv(filename, sep="\t")
        dfQC_formatted.to_csv(iF, sep="\t")
        cellTypeList = ["CD45+", "CD3+ all", "CD4", "CD8", "%NK/CD45", "%B/CD45", "%Lymphocytes/CD45", "HLADR+, Lymphocytes-"]
        compare_Percentages(filename, dfQC, cellTypeList)
        filename.write("\n\n")
        compare_Percentages(iF, dfQC, cellTypeList)
        
        includeAnnotations(filename, df, s)
        includeAnnotations(iF, df, s)
        iF.close()
        
def calculateCD8AndCD4(dfCD):
    """Returns a dataframe with true percentages of CD8s and CD4s."""
    df_Samp = set(dfCD["Sample Name"].values)                
    sampleSet = set()
    finalDF = pd.DataFrame(columns=["Sample Name", "Stain Name", "Cell Type", "Percentage"])
    
    # Building a set of names for each individual sample
    for sn in df_Samp:
        name = sn.split("_")[0] + "_" + sn.split("_")[1]
        sampleSet.add(name)
        
    for sampName in sampleSet:
        onlySample = dfCD[dfCD["Sample Name"].str.contains(sampName)]
        stains = set(onlySample["Stain Name"].values)
        for s in stains:
            onlyStain = onlySample[onlySample.loc[:, "Stain Name"] == s]
            if len(onlyStain) > 1:
                stainName = s
                ct = dfCD["Cell Type"].values[0]
                if("CD8a+" in ct):
                    ct = "CD8"
                elif("CD4+" in ct):
                    ct = "CD4"
                totalCount = np.sum(onlyStain["Count"].values)
                if (totalCount > 0 and onlyStain["Count"].values[0] > 0):
                    totalParent = (onlyStain["Count"].values[0])/((onlyStain["Percentage of Parent"].values[0])/100)
                else:
                    totalParent = (onlyStain["Count"].values[1])/((onlyStain["Percentage of Parent"].values[1])/100)
                percentage = (totalCount/totalParent)*100
                finalDF.loc[finalDF.shape[0]] = [sampName, stainName, ct, percentage]
                
            else:
                stainName = s
                ct = dfCD["Cell Type"].values[0]
                if("CD8a+" in ct):
                    ct = "CD8"
                elif("CD4+" in ct):
                    ct = "CD4"
                p = dfCD["Percentage of Parent"].values[0]
                finalDF.loc[finalDF.shape[0]] = [sampName, stainName, ct, p]
    return finalDF
    
def formatQCSheet(sampleName, dfDict):
    """Formats all the population information from flojo into QC standard format."""
    QC = pd.DataFrame(columns=["Sample Name", "Stain Name", "Cell Type", "Percentage"])
    for dfName in dfDict.keys():
        if dfName == "%NK/CD45":
            for d in dfDict[dfName]:
                hasSample = d[d["Sample Name"].str.contains(sampleName, regex=False)]
                ct = dfName
                stainList = list(hasSample["Stain Name"].values)
                percentList = list(hasSample["Percentage of Grandparent"].values)
                for s in range(len(stainList)):
                    QC.loc[QC.shape[0]] = [sampleName, stainList[s], ct, percentList[s]]
        elif dfName == "%B/CD45":
            df = dfDict[dfName]
            hasSample = df[df["Sample Name"].str.contains(sampleName, regex=False)]
            ct = dfName
            stainList = list(hasSample["Stain Name"].values)
            percentList = list(hasSample["Percentage of Grandparent"].values)
            for s in range(len(stainList)):
                QC.loc[QC.shape[0]] = [sampleName, stainList[s], ct, percentList[s]]
                
        elif dfName == "%CD8/CD3" or dfName == "%CD4/CD3":
            df = dfDict[dfName]
            hasSample = df[df["Sample Name"].str.contains(sampleName, regex=False)]
            QC = QC.append(hasSample, ignore_index=True)
        else:
            df = dfDict[dfName]
            hasSample = df[df["Sample Name"].str.contains(sampleName, regex=False)]
            ct = dfName
            stainList = list(hasSample["Stain Name"].values)
            percentList = list(hasSample["Percentage of Parent"].values)
            for s in range(len(stainList)):
                QC.loc[QC.shape[0]] = [sampleName, stainList[s], ct, percentList[s]]
    return QC
    
def compare_Percentages(filename, df_Query, cellTypeList):
    """Takes in a list of samples and compares their percentages to see if any have a percent difference > 10%. Writes a note below row if any conflicting stains have > 10%"""
    notes = list()
    for ct in cellTypeList:
        df_Match_Stain = df_Query[df_Query.loc[:, "Cell Type"] == ct]
        if len(df_Match_Stain) > 1:
            df_Match_Stain = df_Match_Stain.sort_values(["Stain Name"])
            stain_Comparison = list(df_Match_Stain.loc[:, "Stain Name"])
            count_Comparison = list(df_Match_Stain.loc[:, "Percentage"])   # Inputting the counts into a list for comparison
            
            # Tricky subtraction, subtracts all the values in the list from each other, if greater than 10% than report
            for count in range(len(count_Comparison)):
                for count_forwards in range(count, len(count_Comparison)):
                    diff = math.fabs(float(count_Comparison[count]) - float(count_Comparison[count_forwards]))
                    #print(str(float(count_Comparison[count])) + "-" + str(float(count_Comparison[count_forwards])))
                    if (diff > 10):
                        notes_Stmt = ct + " " + str(stain_Comparison[count]) + " ," + str(stain_Comparison[count_forwards]) + "\tPercent Difference: " + ("%.3f" % diff)
                        notes.append(notes_Stmt)
    if(len(notes) > 0):
        filename.write("\n")
        for n in notes:
            filename.write(n)
            filename.write("\n")

def includeAnnotations(filename, df, sampleName):
    dfSample = df[df["Sample Name"].str.contains(sampleName)]
    dfAnnotation = dfSample[dfSample.loc[:, "Annotation"].notnull()]
    annotateTextList = list()
    if len(dfAnnotation) > 0:
        stainOnly = sorted(set(dfAnnotation["Stain Name"].values))
        for s in stainOnly:
            stainOnlyDF = dfAnnotation[dfAnnotation.loc[:, "Stain Name"] == s]
            cellTypeList = list(stainOnlyDF["Cell Type"].values)
            for c in cellTypeList:
                annotateText = c + " in " + s + " has an N/A."
                annotateTextList.append(annotateText)
                
    if(len(annotateTextList) > 0):
        for n in annotateTextList:
            filename.write(n)
            filename.write("\n")  
    filename.write("\n")
            
# Function made purely for calculating for Stain 2
def calculations(df_C):
    """Calculates lymphocyte and myeloid percents."""
    end_Set = set()
    df_Samples = set(df_C["Sample Name"].values)
    lymphocytes_Dict = dict()
    myeloids_Dict = dict()
    calculations_List = list()
    
    # Building a set of names for each individual sample
    for sn in df_Samples:
        name = sn.split("_")[0] + "_" + sn.split("_")[1]
        end_Set.add(name)
    
    for sample in end_Set:
        df_Match_Name = df_C[df_C["Sample Name"].str.contains(sample)]
        cd45_List = list(df_Match_Name[df_Match_Name["Cell Type"] == "CD45+"]["Count"].values)
        cd45 = float(cd45_List[0])
            
        # For lymphocytes
        b_cell = df_Match_Name["Cell Type"] == "B-cells"
        NK = df_Match_Name["Cell Type"] == "HLADR-, CD3-, CD56+ (NK)"
        CD3 = df_Match_Name["Cell Type"] == "CD3+ all"
        
        if (cd45 != 0):             
            BNKCD3 = pd.to_numeric(df_Match_Name[b_cell | NK | CD3]["Count"]).sum()
            lymphocytes_Dict[sample] = float(BNKCD3/cd45) * 100
        else:
            lymphocytes_Dict[sample] = "N/A"
            
        # For myeloids
        HLADR_List = list(df_Match_Name[df_Match_Name["Cell Type"] == "CD3-, HLADR+"]["Percentage of Parent"].values)
        b_cell_List = list(df_Match_Name[b_cell]["Percentage of Parent"].values)
        hladr = float(HLADR_List[0])/100
        b_cell_Percent = float(b_cell_List[0])/100            
        myeloids_Dict[sample] = float((1 - b_cell_Percent) * hladr) * 100
        
    calculations_List.append(lymphocytes_Dict)
    calculations_List.append(myeloids_Dict)
    return calculations_List
      
def deleteCols(df, colList):
    """Deletes all columns in a list from the dataframe."""
    for c in colList:
        del df[c]
    return df
    
if __name__ == '__main__':

    # Formatting db
    df_all_concat = pd.read_csv("populationDB.csv")
    df = deleteCols(df_all_concat, ["ID", "Parent ID"])
    
    # Makes a folder containing the query outputs
    pop_file = "QC/QC_all_wsps.tsv"
    os.makedirs(os.path.dirname(pop_file), exist_ok=True)
    f = open(pop_file, "w")
    make_QC_Sheet(f, df)
    f.close()
    
    
            