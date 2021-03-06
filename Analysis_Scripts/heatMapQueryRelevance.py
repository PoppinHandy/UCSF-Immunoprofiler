# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:48:17 2017

@author: andyp
"""
import os
import seaborn as sns
import pandas as pd
import glob
import matplotlib.pyplot as plt

def readinMFIDF(location):
    mfiDFList = list()
    for file in glob.glob(location):
        mfiDFList.append(pd.read_csv(file, sep="\t", header=0))
    mfiDF = pd.concat(mfiDFList)
    return mfiDF
    
def setUpFolders(folderList):
    for f in folderList:
        os.makedirs(os.path.dirname(f), exist_ok=True)
 
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
        fig, ax = plt.subplots(figsize=(df.shape[1], df.shape[0]))        
        cbar = fig.add_axes([.905, .6, .02, .1])
        sns.heatmap(df, ax = ax, cmap="Blues", cbar_ax=cbar, cbar=True)
        ax.set_title(name, fontsize=20)
        pdfName = pdfFolder + name + ".pdf"
        plt.yticks(rotation=0, fontsize=15)
        plt.xticks(rotation=90, fontsize=15)
        ax.set_xlabel("") # Adds a label called Parent-None without this
        fig.savefig(pdfName, bbox_inches='tight')
        plt.close()
        
def mfiToHeatmap(mfiDF, cellTypeList, mfiTumorFolder, mfiNormalFolder, heatmapTumorFolder, heatmapNormalFolder):
    """Pipeline that creates MFI Tables, Pearsons and heatmap from inputted data, separted by Tumor and Normal, filtered out by Cell Types"""
    stains = set(mfiDF["Stain Name"].values)
    
    # Removing stains in sample names to avoid redundancy
    df_Samp = list(mfiDF["Sample Name"].values)                
    sampleSet = list()
    for sn in df_Samp:
        name = sn.split("_")[0] + "_" + sn.split("_")[1]
        sampleSet.append(name)
    mfiDF.loc[:, "Sample Name"] = sampleSet

    for st in stains:
        onlyStain = mfiDF[mfiDF.loc[:, "Stain Name"] == st]
        
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
                convertToHeatmap(tumor, st, heatmapTumorFolder)
                
            else:
                tumor = pd.concat(tumorList, axis=1)
                f = mfiTumorFolder + st + ".csv"
                fi = open(f, 'w')
                tumor.to_csv(fi)
                fi.close()
                convertToHeatmap(tumor, st, heatmapTumorFolder)
                
        if(len(normalList) > 0):
            if(len(normalList) == 1):
                normal = normalList[0]
                f = mfiNormalFolder + st + ".csv"
                fi = open(f, 'w')
                normal.to_csv(fi)
                fi.close()
                convertToHeatmap(normal, st, heatmapNormalFolder)
                
            else:
                normal = pd.concat(normalList, axis=1)
                f = mfiNormalFolder + st + ".csv"
                fi = open(f, 'w')
                normal.to_csv(fi)
                fi.close()
                convertToHeatmap(normal, st, heatmapNormalFolder)
                
if __name__ == "__main__":
    popOption = input("Make sure the files are from the Separate By Stain Folder generated by parsePopulations script. Select population? [y/n] (if n, then program does all populations in each stain) ")
    if (popOption == "y"):
        ct = input("What populations are you looking for? (list names separated by semicolons (;), no spaces in between semicolons) ")
        cellTypeList = ct.split(";")
        mfiLocation = "MFI/*.tsv"
        try:
            mfiDF = readinMFIDF(mfiLocation)
        except FileNotFoundError:
            print("File of MFIs not found, make sure they are in a folder labeled MFI in the same location as this script.")
        else:
            majorFolder = "% Composition vs Sample Heat Map Output/"
            mfiTumorFolder = majorFolder + "MFI Tables/Tumor/"
            mfiNormalFolder = majorFolder + "MFI Tables/Normal/"
            heatmapTumorFolder = majorFolder + "Heat Maps/Tumor/"
            heatmapNormalFolder = majorFolder + "Heat Maps/Normal/"
            setUpFolders([majorFolder, mfiTumorFolder, mfiNormalFolder, heatmapTumorFolder, heatmapNormalFolder])
            mfiToHeatmap(mfiDF, cellTypeList, mfiTumorFolder, mfiNormalFolder, heatmapTumorFolder, heatmapNormalFolder)
            
    elif(popOption == "n"):
        mfiLocation = "MFI/*.tsv"
        try:
            mfiDF = readinMFIDF(mfiLocation)
            if 'HLADR HI' in mfiDF.columns:
                if 'LIVEDEAD_HLADR HI' in mfiDF.columns:
                    mfiDF['HLADR HI'].fillna(mfiDF['LIVEDEAD_HLADR HI'], inplace=True)
                    del mfiDF['LIVEDEAD_HLADR HI']
            
        except FileNotFoundError:
            print("File of MFIs not found, make sure they are in a folder labeled MFI in the same location as this script, and that each file is separated by stain and sample.")
        else:
            cellTypeList = set(mfiDF["Parent"].values)         
            majorFolder = "% Composition vs Sample Heat Map Output/"
            mfiTumorFolder = majorFolder + "MFI Tables/Tumor/"
            mfiNormalFolder = majorFolder + "MFI Tables/Normal/"
            heatmapTumorFolder = majorFolder + "Heat Maps/Tumor/"
            heatmapNormalFolder = majorFolder + "Heat Maps/Normal/"
            setUpFolders([majorFolder, mfiTumorFolder, mfiNormalFolder, heatmapTumorFolder, heatmapNormalFolder])
            mfiToHeatmap(mfiDF, cellTypeList, mfiTumorFolder, mfiNormalFolder, heatmapTumorFolder, heatmapNormalFolder)
    else:
        print("You put in an incorrect response! Please enter y or n next time.")
        raise