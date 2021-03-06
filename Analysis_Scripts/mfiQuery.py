# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 12:52:46 2017

@author: andyp
"""

import pandas as pd
import glob
import os

def readinMFIDF(location):
    mfiDFList = list()
    for file in glob.glob(location):
        mfiDFList.append(pd.read_csv(file, sep="\t", header=0))
    mfiDF = pd.concat(mfiDFList)
    return mfiDF

def queryMFI(df, query, majorFolder, filename):
    f = open(filename, 'w')
    for q in query:
        f.write(q)
        f.write("\n")
        q = q.upper() + " HI"
        query = df[df[q].notnull()]
        query = query.dropna(axis=1, how='all')
        query = pd.pivot_table(query, index="Sample Name")
        query.to_csv(f)
        
        # To make individual files per query
        file2 = majorFolder + q + ".csv"
        f2 = open(file2, 'w')
        query.to_csv(f2)
        f2.close()
        
        f.write("\n")
    f.close()
        
if __name__ == "__main__":
    ct = input("Make sure the files are from the Separate By Stain Folder generated by parsePopulations script. What MFI gates are you looking for? (list names separated by semicolons (;), no spaces in between semicolons, DO NOT INCLUDE HI) ")
    cellTypeList = ct.split(";")
    fn = input("What do you want the file (with all the queries together) to be named? ")
    mfiLocation = "MFI/*.tsv"
    try:
        mfiDF = readinMFIDF(mfiLocation)
    except FileNotFoundError:
        print("File of MFIs not found, make sure they are in a folder labeled MFI in the same location as this script.")
    else:
        majorFolder = "MFI Query Output/"
        os.makedirs(os.path.dirname(majorFolder), exist_ok=True)
        
        # Weird naming issue
        if 'HLADR HI' in mfiDF.columns:
            if 'LIVEDEAD_HLADR HI' in mfiDF.columns:
                mfiDF['HLADR HI'].fillna(mfiDF['LIVEDEAD_HLADR HI'], inplace=True)
                del mfiDF['LIVEDEAD_HLADR HI']
        fileLoc = majorFolder + fn + ".csv"
        queryMFI(mfiDF, cellTypeList, majorFolder, fileLoc)
        