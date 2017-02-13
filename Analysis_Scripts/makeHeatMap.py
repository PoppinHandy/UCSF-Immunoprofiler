# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 12:03:24 2017

@author: andyp
"""

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import glob
from scipy.stats.stats import pearsonr
import numpy as np
import glob
import os

def readCSV(fileLocation):
    pearsonDF = pd.read_csv(fileLocation, index_col='Unnamed: 0', sep=",")
    return pearsonDF
    
def convertToHeatmap(df, name, isTumor, pdfFolder):
    plt.ioff()
    if len(df) > 2:
        if(isTumor):
            os.makedirs(os.path.dirname(pdfFolder), exist_ok=True)
            fig, ax = plt.subplots(figsize=(10,10))
            sns.heatmap(df, cmap="RdBu_r")
            ax.set_title(name, fontsize=15)
            pdfName = pdfFolder + name + ".pdf"
            fig.savefig(pdfName, bbox_inches='tight')
            plt.close()
            
        else:
            os.makedirs(os.path.dirname(pdfFolder), exist_ok=True)
            fig, ax = plt.subplots(figsize=(10,10))
            sns.heatmap(df, cmap="RdBu_r")
            ax.set_title(name, fontsize=15)
            pdfName = pdfFolder + name + ".pdf"
            fig.savefig(pdfName, bbox_inches='tight')
            plt.close()
    
def convertToPearson(df, cellType, stainName, mfiFolder):
    """Takes in a dataframe and converts the values to pearson correlation values."""
    colLen = list()
    for c in df.columns:
        colLen.append(len(df[c]))
    maxValue = np.max(colLen)
    for c in df.columns:
        if (len(df[c]) != maxValue):
            del df[c]
            
    col = list(df.columns.values)
    pearsonDict = dict()
    pValueDict = dict()
    for c in range(len(col)):
        name = col[c]
        if name not in pearsonDict.keys():
            pearsonDict[name] = dict()
            pValueDict[name] = dict()
        for c2 in range(len(col)):
            if (c2 != c):
                name2 = col[c2]
                x = df[name].values
                y = df[name2].values
            
                # Filtering out any blanks since pearsonr dislikes nan
                nas = np.logical_or(np.isnan(x), np.isnan(y))
                pearsonDict[name][name2], pValueDict[name][name2] = pearsonr(x[~nas], y[~nas])
    pearsonDF = pd.DataFrame(pearsonDict, index=pearsonDict.keys(), columns=pearsonDict.keys())
    #pearsonDF = pearsonDF[list(pearsonDF.columns.values)]
    pearsonDF = pearsonDF.reindex(index=list(pearsonDF.columns.values))
    pearsonDF = pearsonDF.dropna(axis=1, how='all')
    pearsonDF = pearsonDF.dropna(axis=0, how='all')
    f = mfiFolder + "Correlation of Percent MFI in " + cellType + " " + stainName  + ".csv"
    file = open(f, 'w')
    pearsonDF.to_csv(file)
    file.close()
    return pearsonDF
    
def readinMFIDF(location):
    mfiDFList = list()
    for file in glob.glob(location):
        mfiDFList.append(pd.read_csv(file, sep="\t", header=0))
    mfiDF = pd.concat(mfiDFList)
    return mfiDF
   
def createAllMFIs(mfiDF, mfiTumorFolder, mfiNormalFolder, pearsonTumorFolder, pearsonNormalFolder, heatmapTumorFolder, heatmapNormalFolder):
    """Creates raw MFI Tables, separated by Tumor and Normal samples, of all populations."""
    stains = set(mfiDF["Stain Name"].values)
    for st in stains:
        onlyStain = mfiDF[mfiDF.loc[:, "Stain Name"] == st]
        subpops = set(onlyStain["Parent"].values)
        for sp in subpops:
            onlySubpops = onlyStain[onlyStain.loc[:, "Parent"] == sp] 
            if len(onlySubpops) > 1:
                onlySubpops = deleteCols(onlySubpops, ["Parent", "Grandparent", "Stain Name"])
                tumorDF = onlySubpops[onlySubpops["Sample Name"].str.contains("_T")]
                normalDF = onlySubpops[onlySubpops["Sample Name"].str.contains("_N")]
                tumorDF = tumorDF.dropna(axis=1, how='all')
                normalDF = normalDF.dropna(axis=1, how='all')
                
                if (len(tumorDF) > 0):
                    finalTumorTable = pd.pivot_table(tumorDF, index="Sample Name")
                    outputMFITable(finalTumorTable, sp, st, True, mfiTumorFolder)
                    dfPearson = convertToPearson(finalTumorTable, sp, st, pearsonTumorFolder)
                    heatMapName = "Correlation of Percent MFI in " + sp + " " + st + " Tumor"
                    convertToHeatmap(dfPearson, heatMapName, True, heatmapTumorFolder)
                    
                if(len(normalDF) > 2):
                    #finalNormalTable = normalDF.reindex(index=normalDF["Sample Name"].values)
                    finalNormalTable = pd.pivot_table(normalDF, index="Sample Name")
                    outputMFITable(finalNormalTable, sp, st, False, mfiNormalFolder)
                    dfPearson = convertToPearson(finalNormalTable, sp, st, pearsonNormalFolder)
                    heatMapName = "Correlation of Percent MFI in " + sp + " " + st + " Normal"
                    convertToHeatmap(dfPearson, heatMapName, False, heatmapNormalFolder)

def createMFIByCT(mfiDF, cellType, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder):
    """Creates raw MFI Tables, separted by Tumor and Normal, filtered out by Cell Types"""
    mfiDF = mfiDF[mfiDF.loc[:, "Parent"] == cellType]
    stains = set(mfiDF["Stain Name"].values)
    for st in stains:
        onlyStain = mfiDF[mfiDF.loc[:, "Stain Name"] == st]
        if len(onlyStain) > 1:
            onlyStain = deleteCols(onlyStain, ["Parent", "Grandparent", "Stain Name"])
            tumorDF = onlyStain[onlyStain["Sample Name"].str.contains("_T")]
            normalDF = onlyStain[onlyStain["Sample Name"].str.contains("_N")]
            tumorDF = tumorDF.dropna(axis=1, how='all')
            normalDF = normalDF.dropna(axis=1, how='all')
            if (len(tumorDF) > 0):
                finalTumorTable = pd.pivot_table(tumorDF, index="Sample Name")
                outputMFITable(finalTumorTable, cellType, st, True, mfiTumorFolder)
                dfPearson = convertToPearson(finalTumorTable, cellType, st, pearsonTumorFolder)
                heatMapName = "Correlation of Percent MFI in " + cellType + " " + st + " Tumor"
                convertToHeatmap(dfPearson, heatMapName, True, heatmapTumorFolder)
                
            if(len(normalDF) > 2):
                finalNormalTable = pd.pivot_table(normalDF, index="Sample Name")
                outputMFITable(finalNormalTable, cellType, st, False, mfiNormalFolder)
                dfPearson = convertToPearson(finalNormalTable, cellType, st, pearsonNormalFolder)
                heatMapName = "Correlation of Percent MFI in " + cellType + " " + st + " Normal"
                convertToHeatmap(dfPearson, heatMapName, False, heatmapNormalFolder)
    
def deleteCols(df, colList):
    """Deletes all columns in a list from the dataframe."""
    for c in colList:
        del df[c]
    return df
    
def outputMFITable(df, parentName, stainName, isTumor, mfiFolder):
    if (isTumor):
        folder = mfiFolder
        fileHeader = folder + "Tumor " + parentName + " " + stainName + " Function" + ".csv"
        f = open(fileHeader, 'w')
        df.to_csv(f)
        f.close()
    else:
        folder = mfiFolder
        fileHeader = folder + "Normal " + parentName + " " + stainName + " Function" + ".csv"
        f = open(fileHeader, 'w')
        df.to_csv(f)
        f.close()
    
def setUpFolders(folderList):
    for f in folderList:
        os.makedirs(os.path.dirname(f), exist_ok=True)
        
if __name__ == "__main__":
    
    
    mfiOption = input("Build tables from MFIs? [y/n] ")
    popOption = input("Select population? [y/n] (if n, then program does all populations) ")
    if (mfiOption == "y"):
        if (popOption == "y"):
            ct = input("What populations are you looking for? (list names separated by semicolons (;), no spaces in between semicolons) ")
            cellTypeList = ct.split(";")
            mfiLocation = "MFI/*.tsv"
            try:
                mfiDF = readinMFIDF(mfiLocation)
            except FileNotFoundError:
                print("File of MFIs not found, make sure they are in a folder labeled MFI in the same location as this script.")
            else:
                majorFolder = "Analysis Output For Queries/"
                mfiTumorFolder = majorFolder + "MFI Tables/Tumor/"
                mfiNormalFolder = majorFolder + "MFI Tables/Normal/"
                pearsonTumorFolder = majorFolder + "Pearson Tumor/"
                pearsonNormalFolder = majorFolder + "Pearson Normal/"
                heatmapTumorFolder = majorFolder + "Heat Maps/Tumor/"
                heatmapNormalFolder = majorFolder + "Heat Maps/Normal/"
                setUpFolders([majorFolder, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder])
                
                for cellType in cellTypeList:
                    # Pull out by stain, then by subpopulation
                    createMFIByCT(mfiDF, cellType, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder)
                    
        elif (popOption == "n"):
            mfiLocation = "MFI/*.tsv"
            try:
                mfiDF = readinMFIDF(mfiLocation)
            except FileNotFoundError:
                print("File of MFIs not found, make sure they are in a folder labeled MFI in the same location as this script.")
            else:
                majorFolder = "Analysis Output/"
                mfiTumorFolder = majorFolder + "MFI Tables/Tumor/"
                mfiNormalFolder = majorFolder + "MFI Tables/Normal/"
                pearsonTumorFolder = majorFolder + "Pearson Tumor/"
                pearsonNormalFolder = majorFolder + "Pearson Normal/"
                heatmapTumorFolder = majorFolder + "Heat Maps/Tumor/"
                heatmapNormalFolder = majorFolder + "Heat Maps/Normal/"
                setUpFolders([majorFolder, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder])
                
                # Pull out by stain, then by subpopulation
                createAllMFIs(mfiDF, mfiTumorFolder, mfiNormalFolder, pearsonTumorFolder, pearsonNormalFolder, heatmapTumorFolder, heatmapNormalFolder)
                #files = glob.glob("Analysis Outputs/MFI Tables/Tumor/*.csv")

        else:
            print("You put in an incorrect response! Please enter y or n next time.")
            raise
            
    elif (mfiOption == "n"):
        print("Will be implemented at a later point.")
#==============================================================================
#                 for f in files:
#                     name = f.split('.')[0]
#                     df = readCSV(f)
#                     fPDF = "Correlation of Percent MFI in " + name
#                     dfPearson = convertToPearson(df, f)
#                     convertToHeatmap(dfPearson, fPDF)
#==============================================================================