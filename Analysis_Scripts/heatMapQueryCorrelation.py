# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 15:03:16 2017

@author: andyp
"""

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import glob
from scipy.stats.stats import pearsonr
import numpy as np
import os

def deleteCols(df, colList):
    """Deletes all columns in a list from the dataframe."""
    for c in colList:
        del df[c]
    return df
        
def convertToHeatmap(df, name, pdfFolder):
    """Converts a pearson correlation chart into a heatmap."""
    plt.ioff()
    if len(df) > 2:
        os.makedirs(os.path.dirname(pdfFolder), exist_ok=True)
        fig, ax = plt.subplots(figsize=(len(df), len(df)))
        cbar = fig.add_axes([.905, .6, .02, .1])
        sns.heatmap(df, ax = ax, cmap="RdBu_r", cbar_ax=cbar, cbar=True)
        ax.set_title(name, fontsize=20)
        pdfName = pdfFolder + name + ".pdf"
        plt.yticks(rotation=0, fontsize=15)
        plt.xticks(rotation=90, fontsize=15)
        fig.savefig(pdfName, bbox_inches='tight')
        plt.close()
    
    
def mfiToHeatmap(mfiDF, cellTypeList, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder):
    """Pipeline that creates MFI Tables, Pearsons and heatmap from inputted data, separted by Tumor and Normal, filtered out by Cell Types"""
    stains = set(mfiDF["Stain Name"].values)
    for st in stains:
        onlyStain = mfiDF[mfiDF.loc[:, "Stain Name"] == st]
        
        # Weird naming issue
        onlyStain['HLADR HI'].fillna(onlyStain['LIVEDEAD_HLADR HI'], inplace=True)
        del onlyStain['LIVEDEAD_HLADR HI']
        
        tumorList = list()
        normalList = list()
        
        for cellType in cellTypeList:
            onlyCell = onlyStain[onlyStain.loc[:, "Parent"] == cellType]
            if len(onlyCell) > 1:
                onlyCell = deleteCols(onlyCell, ["Grandparent", "Stain Name"])
                tumorDF = onlyCell[onlyCell["Sample Name"].str.contains("_T")]
                normalDF = onlyCell[onlyCell["Sample Name"].str.contains("_N")]
                tumorDF = tumorDF.dropna(axis=1, how='all')
                normalDF = normalDF.dropna(axis=1, how='all')
                
                # Have to make sure more than one valid entry for each df or else program will complain once it tries to calculate pearson's
                if (len(tumorDF) > 2):
                    finalTumorTable = pd.pivot_table(tumorDF, index=["Sample Name"], columns=["Parent"])
                    finalTumorTable = finalTumorTable.swaplevel(axis=1)
                    tumorList.append(finalTumorTable)
                    
                if(len(normalDF) > 2):
                    finalNormalTable = pd.pivot_table(normalDF, index=["Sample Name"], columns=["Parent"])
                    finalNormalTable = finalNormalTable.swaplevel(axis=1)
                    normalList.append(finalNormalTable)
                    
        if (len(tumorList) > 0):
            if(len(tumorList) == 1):
                tumor = tumorList[0]
                f = mfiTumorFolder + st + ".csv"
                fi = open(f, 'w')
                tumor.to_csv(fi)
                fi.close()
                pearsonDF = convertToPearson(tumor, st, pearsonTumorFolder)
                convertToHeatmap(pearsonDF, st, heatmapTumorFolder)
                
            else:
                tumor = pd.concat(tumorList, axis=1)
                f = mfiTumorFolder + st + ".csv"
                fi = open(f, 'w')
                tumor.to_csv(fi)
                fi.close()
                pearsonDF = convertToPearson(tumor, st, pearsonTumorFolder)
                convertToHeatmap(pearsonDF, st, heatmapTumorFolder)
                
        if(len(normalList) > 0):
            if(len(normalList) == 1):
                normal = normalList[0]
                f = mfiNormalFolder + st + ".csv"
                fi = open(f, 'w')
                normal.to_csv(fi)
                fi.close()
                pearsonDF = convertToPearson(normal, st, pearsonNormalFolder)
                convertToHeatmap(pearsonDF, st, heatmapNormalFolder)
                
            else:
                normal = pd.concat(normalList, axis=1)
                f = mfiNormalFolder + st + ".csv"
                fi = open(f, 'w')
                normal.to_csv(fi)
                fi.close()
                pearsonDF = convertToPearson(normal, st, pearsonNormalFolder)
                convertToHeatmap(pearsonDF, st, heatmapNormalFolder)

def convertToPearson(df, stainName, mfiFolder):
    """Takes in a dataframe and converts the values to pearson correlation values."""
    colLen = list()
    for c in df.columns:
        colLen.append(len(df[c]))
    maxValue = np.max(colLen)
    for c in df.columns:
        if (len(df[c]) != maxValue):
            del df[c]
            
    col = list(df.columns.values) # format is (cell type, mfi gate name)
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
                
    # dictionary keys are tuples read in from multiindex
    mdx = pd.MultiIndex.from_tuples(pearsonDict.keys())
    pearsonDF = pd.DataFrame(pearsonDict, index=mdx, columns=mdx)
    pearsonDF = pearsonDF.sort_index(axis=1)
    pearsonDF = pearsonDF.sort_index(axis=0)
    #pearsonDF = pearsonDF[list(pearsonDF.columns.values)]
    pearsonDF = pearsonDF.reindex(index=list(pearsonDF.columns.values))
    pearsonDF = pearsonDF.dropna(axis=1, how='all')
    pearsonDF = pearsonDF.dropna(axis=0, how='all')
    f = mfiFolder + "Correlation of Percent MFI in " + stainName  + ".csv"
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
    
def setUpFolders(folderList):
    for f in folderList:
        os.makedirs(os.path.dirname(f), exist_ok=True)
        
if __name__ == "__main__":
    popOption = input("Select population? [y/n] (if n, then program does all populations in each stain) ")
    if (popOption == "y"):
        ct = input("What populations are you looking for? (list names separated by semicolons (;), no spaces in between semicolons) ")
        cellTypeList = ct.split(";")
        #cellTypeList = ["CD45+", "CD3+ all"]
        mfiLocation = "MFI/*.tsv"
        try:
            mfiDF = readinMFIDF(mfiLocation)
        except FileNotFoundError:
            print("File of MFIs not found, make sure they are in a folder labeled MFI in the same location as this script.")
        else:
            majorFolder = "Heat Map Query Output/"
            mfiTumorFolder = majorFolder + "MFI Tables/Tumor/"
            mfiNormalFolder = majorFolder + "MFI Tables/Normal/"
            pearsonTumorFolder = majorFolder + "Pearson Tumor/"
            pearsonNormalFolder = majorFolder + "Pearson Normal/"
            heatmapTumorFolder = majorFolder + "Heat Maps/Tumor/"
            heatmapNormalFolder = majorFolder + "Heat Maps/Normal/"
            setUpFolders([majorFolder, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder])
            mfiToHeatmap(mfiDF, cellTypeList, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder)
            
    elif(popOption == "n"):
        mfiLocation = "MFI/*.tsv"
        try:
            mfiDF = readinMFIDF(mfiLocation)
        except FileNotFoundError:
            print("File of MFIs not found, make sure they are in a folder labeled MFI in the same location as this script, and that each file is separated by stain and sample.")
        else:
            cellTypeList = set(mfiDF["Parent"].values)         
            majorFolder = "Heat Map Query Output/"
            mfiTumorFolder = majorFolder + "MFI Tables/Tumor/"
            mfiNormalFolder = majorFolder + "MFI Tables/Normal/"
            pearsonTumorFolder = majorFolder + "Pearson Tumor/"
            pearsonNormalFolder = majorFolder + "Pearson Normal/"
            heatmapTumorFolder = majorFolder + "Heat Maps/Tumor/"
            heatmapNormalFolder = majorFolder + "Heat Maps/Normal/"
            setUpFolders([majorFolder, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder])
            mfiToHeatmap(mfiDF, cellTypeList, mfiTumorFolder, mfiNormalFolder, pearsonNormalFolder, pearsonTumorFolder, heatmapTumorFolder, heatmapNormalFolder)
    else:
        print("You put in an incorrect response! Please enter y or n next time.")
        raise