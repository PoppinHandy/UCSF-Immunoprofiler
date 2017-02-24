# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 15:49:38 2017

@author: andyp
"""

import xml.etree.ElementTree as ET
import os
import glob
import math
import pandas as pd

def generatePopulationsSummary(workspace_path, common_count_limit, rare_count_limit, include_mfis):
    
    """Takes in the path of the workspace and produces a pandas database containing information about each population pulled from the wsp xml. Outputs the database."""
    # Importing the workspace in
    tree = ET.parse(workspace_path)
    root = tree.getroot()
    
    # Creating the headers for tsv
    stain_Table = pd.DataFrame(columns=["Stain Name", "Sample Name", "Cell Type", "Count", "ID", "Parent ID"])
    stain_Table["Count"] = stain_Table["Count"].astype('float64')
    
    # Finding the tube name
    for tube_info in root.findall(".//Sample/SampleNode"):
        tube_Name = tube_info.attrib.get("name")
        # Finds the XML nodes with population information and prints it to file
        for pop in tube_info.findall(".//Population"):
            ids = 0
            parent_id = 0
            stain_group = pop.attrib.get("owningGroup")
            name = pop.attrib.get("name")
            count = int(pop.attrib.get("count"))
            # if false, will not read in mfi gate entries into the dataframe
            if (include_mfis == False):
                if (not " hi" in name):            
                    if (stain_group == ""):
                        if "treg" in tube_Name:
                            stain_group = "Stain 1"
                        elif "nktb" in tube_Name:
                            stain_group = "Stain 2"
                        elif "sort" in tube_Name:
                            stain_group = "Stain 3"
                        elif "dc" in tube_Name:
                            stain_group = "Stain 4"
                        elif "innate" in tube_Name:
                            stain_group = "Stain 5"
                        else:
                            stain_group = "N/A"
                    # Building the percentage composition table
                    for i in pop.findall("Gate"):
                        s = list(i.attrib.values())
                        
                        # xml_Order and checkID is because sometimes, the xml parser switches the order of id and parent id for some weird reason so this is a failsafe to adapt to that weird bug
                        xml_Order = list(i.attrib.keys())
                        checkID_Order = xml_Order[0].split("}")[1]
                        if (len(s) == 1):
                            ids = s[0]
                        else:
                            if (checkID_Order == "id"):
                                ids = s[0]
                                parent_id = s[1]
                            else:
                                ids = s[1]
                                parent_id = s[0]
                                    # Adding data to be displayed to the final table
                    stain_Table.loc[stain_Table.shape[0]] = [stain_group, tube_Name, name, count, ids, parent_id]
                
    stain_Table = stain_Table.sort_values(["Sample Name", "Stain Name"]) # Ordering the table to prevent future life stresses
    
    # Partitioning percentages
    stain_Table["Parent"] = stain_Table.apply(lambda row: defineParents(stain_Table, row, 1), axis=1)
    stain_Table["Percentage of Parent"] = definePercentageComposition(stain_Table, stain_Table["Parent ID"], stain_Table["Count"], 1)
    stain_Table["Grandparent"] = stain_Table.apply(lambda row: defineParents(stain_Table, row, 2), axis=1)
    stain_Table["Percentage of Grandparent"] = definePercentageComposition(stain_Table, stain_Table["Parent ID"], stain_Table["Count"], 2)
    #print(stain_Table)
    
    # Cleaning up and formatting, setting as category takes up less memory and faster
#==============================================================================
#     stain_Table["Stain Name"] = stain_Table["Stain Name"].astype('category')
#     stain_Table["Parent"] = stain_Table["Parent"].astype('category')
#     stain_Table["Percentage of Parent"] = stain_Table["Percentage of Parent"].astype("float64")
#     stain_Table["Percentage of Grandparent"] = stain_Table["Percentage of Grandparent"].astype("float64")
#     stain_Table["Cell Type"] = stain_Table["Cell Type"].astype('category')
#==============================================================================
    stain_Table = stain_Table[["Sample Name", "Cell Type", "Stain Name", "Parent", "Count", "Percentage of Parent", "ID", "Parent ID", "Grandparent", "Percentage of Grandparent"]] 
    return stain_Table
    
def defineParents(df, row_count, generations):
    """Finds the identity of the parent using parent ids. Applies to all rows of df"""
    
    # Keeps track of IDs as the family finder goes up the tree
    current_ID = row_count["Parent ID"]
    parent_Name = ""
    if row_count["Parent ID"] == 0:
        return "N/A"
    else:
        for i in range(generations):
            if current_ID == 0:
                return "N/A"
            parent_match = df[df.loc[: ,"ID"] == current_ID]
            parent_index = parent_match.index.tolist() # Keeps track of the index for the match (should only be 1)
            current_ID = parent_match.at[parent_index[0], "Parent ID"] # Updates current ID to equal the next parent ID to move up the ancestry tree
            parent_Name = parent_match.at[parent_index[0], "Cell Type"]
        return parent_Name
        
def definePercentageComposition(df, parentID, count, generations):
    """Finds the percentage of the parent count."""
    
    parent_Count = list()   # Keeps parent counts to be divided by count array at the end
    for pid in parentID:
         # Keeps track of IDs as the family finder goes up the tree
        current_ID = pid
        if pid == 0:
            parent_Count.append(0)
        else:
            pc = 0
            for i in range(generations):
                if current_ID == 0:
                    pc = 0
                    break
                parent_match = df[df.loc[: ,"ID"] == current_ID]
                parent_index = parent_match.index.tolist() # Keeps track of the index for the match (should only be 1)
                current_ID = parent_match.at[parent_index[0], "Parent ID"] # Updates current ID to equal the next parent ID to move up the ancestry tree
                pc = parent_match.at[parent_index[0], "Count"] # Updates the parent's count for final calculations
            parent_Count.append(pc)
    percent_Comp = (count/parent_Count) * 100
    return percent_Comp  
    
def filterOutMFIs(df):
    """Takes out any cell types with " hi" in it to separate MFIs from population summary."""
    noMFIs = df[~df["Cell Type"].str.contains(" hi")]
    return noMFIs

def deleteCols(df, colList):
    """Deletes all columns in a list from the dataframe."""
    for c in colList:
        del df[c]
    return df
    
def generateMFISummary(wspPath, mfiFile, mfiPercentFile, mfiCountFile):
    """Generates a summary of mfi values and mfi percentages of cell populations."""                
    # Getting xml parsing ready
    tree = ET.parse(wspPath)
    root = tree.getroot()
    
    # Prepping dataframd and dictionary for info
    mfiTable = pd.DataFrame(columns=["Sample Name", "Stain Name", "Cell Type", "ID", "Parent ID", "Count", "MFI Value"])
    mfiPopulationTable = pd.DataFrame(columns=["Sample Name", "Stain Name", "Cell Type", "MFI ID", "MFI Value"])
    
    # Finding the tube name
    for tube_info in root.findall(".//Sample/SampleNode"):
        tube_Name = tube_info.attrib.get("name")
        
        # Finds the XML nodes with population information
        for pop in tube_info.findall(".//Population"):
            ids = 0
            parent_id = 0
            stain_group = pop.attrib.get("owningGroup")
            cellType = pop.attrib.get("name")
            count = int(pop.attrib.get("count"))
            mfiValue = math.inf
            mfiID = "N/A"
            
            if (stain_group == ""):
                if "treg" in tube_Name:
                    stain_group = "Stain 1"
                elif "nktb" in tube_Name:
                    stain_group = "Stain 2"
                elif "sort" in tube_Name:
                    stain_group = "Stain 3"
                elif "dc" in tube_Name:
                    stain_group = "Stain 4"
                elif "innate" in tube_Name:
                    stain_group = "Stain 5"
                else:
                    stain_group = "N/A"

            # Building the percentage composition table
            for i in pop.findall("Gate"):
                s = list(i.attrib.values())
                
                # xml_Order and checkID is because sometimes, the xml parser switches the order of id and parent id for some weird reason so this is a failsafe to adapt to that weird bug
                xml_Order = list(i.attrib.keys())
                checkID_Order = xml_Order[0].split("}")[1]
                if (len(s) == 1):
                    ids = s[0]
                else:
                    if (checkID_Order == "id"):
                        ids = s[0]
                        parent_id = s[1]
                    else:
                        ids = s[1]
                        parent_id = s[0]
                        # Only mfi gate values read in 
                        
            for mfi in pop.findall("./Subpopulations/Statistic"):
                if (" hi" in cellType):
                    if (mfi.attrib.get('value') == '\ufffd'):
                        mfiValue = math.inf
                    else:
                        mfiValue = float(mfi.attrib.get('value'))
                else:
                    mfiID = mfi.attrib.get('id')
                    if (mfi.attrib.get('value') == '\ufffd'):
                        mfiValue = math.inf
                    else:
                        mfiValue = float(mfi.attrib.get('value'))
            
                        
            # Adding data to be displayed to the final table: ["Sample Name", "Stain Name", "Cell Type", "ID", "Parent ID", "Counts", "MFI Values"]
            mfiTable.loc[mfiTable.shape[0]] = [tube_Name, stain_group, cellType, ids, parent_id, count, mfiValue]
            mfiPopulationTable[mfiPopulationTable.shape[0]] = [tube_Name, stain_group, cellType, mfiID, mfiValue]
            
    # Partitioning percentages
    mfiTable["Parent"] = mfiTable.apply(lambda row: defineParents(mfiTable, row, 1), axis=1)
    mfiTable["Percentage of Parent"] = definePercentageComposition(mfiTable, mfiTable["Parent ID"], mfiTable["Count"], 1)
    mfiTable["Grandparent"] = mfiTable.apply(lambda row: defineParents(mfiTable, row, 2), axis=1)
    
    # Outputting mfi files to both percent and value
    outputMFI(mfiFile, mfiTable, "MFI Value")
    outputMFI(mfiPercentFile, mfiTable, "Percentage of Parent")
    outputMFI(mfiCountFile, mfiTable, "Count")
    
def outputMFI(mfiFolder, mfiTable, col_values):
    """Format mfi table into a readable chart for Prism."""
    mfiTable = mfiTable[mfiTable["Cell Type"].str.contains(" hi")]
    if col_values == "MFI Value":
        mfiTable = deleteCols(mfiTable, ["Percentage of Parent", "Parent ID", "ID", "Count"])
    if col_values == "Percentage of Parent":
        mfiTable = deleteCols(mfiTable, ["MFI Value", "Parent ID", "ID", "Count"])
    if col_values == "Count":
        mfiTable = deleteCols(mfiTable, ["MFI Value", "Parent ID", "ID", "Percentage of Parent"])
        
    # Stain for separation by stain and sample
    stains = ["Stain 1", "Stain 2", "Stain 3", "Stain 4", "Stain 5"]
    df_Samples = set(mfiTable["Sample Name"].values)
    sample_Set = set()
    
    # Building a set of names for each individual sample
    for sn in df_Samples:
        name = sn.split("_")[0] + "_" + sn.split("_")[1]
        sample_Set.add(name)
    
    # Making a new file for each sample 
    for sn in sample_Set:
        mfiFile = mfiFolder + sn + "_" + col_values + ".tsv" 
        f = open(mfiFile, 'w')
        mfiSeparateBySample = mfiTable[mfiTable["Sample Name"].str.contains(sn)]
       # Sorting each category by stain
        for st in stains:
            stainExist = mfiSeparateBySample.loc[:, "Stain Name"] == st
            mfiSeparateByStain = mfiSeparateBySample[stainExist]
            
            if len(mfiSeparateByStain) > 0:
                
                # Adding this to read in files for the heatmap script
                if col_values == "Percentage of Parent":
                    dbStain = "WSP_Outputs/MFIs/MFI ByStain/" 
                    db =  dbStain + sn + "_" + col_values + "_" + st + ".tsv"
                    os.makedirs(os.path.dirname(dbStain), exist_ok=True)
                    f2 = open(db, 'w')    
                    
                    mfiSeparateByStain = mfiSeparateByStain.sort_values(["Parent", "Grandparent"])
                    mfiSeparateByStain = pd.pivot_table(mfiSeparateByStain, index=["Sample Name", "Stain Name", "Parent", "Grandparent"], columns="Cell Type", values=col_values)
                    mfiSeparateByStain.to_csv(f, sep="\t")
                    f.write("\n")

                    mfiSeparateByStain.columns=[c.upper() for c in mfiSeparateByStain.columns]
                    mfiSeparateByStain.to_csv(f2, sep="\t")
                    f2.write("\n")
                    f2.close()
                else:
                    mfiSeparateByStain = mfiSeparateByStain.sort_values(["Parent", "Grandparent"])
                    mfiSeparateByStain = pd.pivot_table(mfiSeparateByStain, index=["Sample Name", "Stain Name", "Parent", "Grandparent"], columns="Cell Type", values=col_values)
                    mfiSeparateByStain.to_csv(f, sep="\t")
                    f.write("\n")
        f.close()
        
        
def mfiSummaryRunner(wspLocation):
    """File setup and method calling for making the mfi summary files"""
    
    mfiFiles = "WSP_Outputs/MFIs/MFI Values/"
    mfiPercentFiles = "WSP_Outputs/MFIs/MFI Percents/"
    mfiCountFiles = "WSP_Outputs/MFIs/MFI Counts/"
    mfiAll = "WSP_Outputs/MFIs/allMFis.tsv"
    
    # Creating new folder
    os.makedirs(os.path.dirname(mfiFiles), exist_ok=True)
    os.makedirs(os.path.dirname(mfiPercentFiles), exist_ok=True)
    os.makedirs(os.path.dirname(mfiCountFiles), exist_ok=True)
    os.makedirs(os.path.dirname(mfiAll), exist_ok=True)
    
    for f in glob.glob(wspLocation):
        generateMFISummary(f, mfiFiles, mfiPercentFiles, mfiCountFiles)
    
    allMFIFiles = list()
    allMFIFiles.append(glob.glob("WSP_Outputs/MFIs/MFI Values/*.tsv"))
    allMFIFiles.append(glob.glob("WSP_Outputs/MFIs/MFI Percents/*.tsv"))
    allMFIFiles.append(glob.glob("WSP_Outputs/MFIs/MFI Counts/*.tsv"))
    
    with open(mfiAll, 'w') as file:
        for outer in allMFIFiles:
            for f in outer:
                splitName = f.split("\\")[1]
                file.write("Showing " + splitName + "\n")
                for line in open(f, 'r'):
                    file.write(line)
        file.write("\n")
    
        
def populationSummaryRunner(wspLocation):
    """File setup and method calling for making the population summary file. Also returns a list of dataframes for the whole file"""
    
    fileName = "WSP_Outputs/wsp_PopulationSummary.tsv"
    
    # Makes a folder containing the wsp outputs
    os.makedirs(os.path.dirname(fileName), exist_ok=True)
        
    # Setting count limits
    common_Count_Limit = math.inf
    rare_Count_Limit = math.inf
    
    # Open file for each workspace
    f_pop = open(fileName, 'w')
    for f in glob.glob(wspLocation):
        
        # Getting the name of the workspace
        wsp_name = os.path.split(f)[1].split(".")[0]
        
        # Annoying file naming system
        f = f.replace('\\','/')
        
        # Writing header
        f_pop.write(str(wsp_name))
        f_pop.write("\n")
        
        # Getting the dataframe for the wsp
        st = generatePopulationsSummary(f, common_Count_Limit, rare_Count_Limit, False)
        formatted_st = deleteCols(st, ["ID", "Parent ID", "Grandparent", "Percentage of Grandparent"])
        #f_pop.write("Population Common Count Limit: " + str(common_Count_Limit) + "\n")
        #f_pop.write("Population Rare Count Limit: " + str(rare_Count_Limit) + "\n")      
        f_pop.write(formatted_st.to_csv(sep="\t", index=False))
        f_pop.write("\n") # End of one wsp section
    f_pop.close()
    
# For importing to other modules
def createGlobalDB(wspLocation):
    dfList = list()
    for f in glob.glob(wspLocation):
        st = generatePopulationsSummary(f, math.inf, math.inf, False)
        dfList.append(st) 
    df_all_concat = pd.concat(dfList)
    df_all_concat = filterOutMFIs(df_all_concat)
    f = open("populationDB.csv", 'w')
    df_all_concat.to_csv(f)   
    
if __name__ == "__main__":    
    wspLocation = "*.wsp"
    populationSummaryRunner(wspLocation)
    mfiSummaryRunner(wspLocation)
    createGlobalDB(wspLocation)
                