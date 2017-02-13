# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:02:52 2017

@author: andyp
"""

import parse_Populations
import os
import pandas as pd

# Query the db for specific rows and cols that the user inputs and report any percentage between stains greater than 10% for QC purposes
def locate_Population(filename, col_name, row_name, df):
    """Takes in a file to be written, column name, and list of row names to query all the wsp dataframes. Writes answer to query to file in tsv format as well as noting any populations that have a percent difference greater than 10% between stains."""
    
    row_name = row_name.split(";")
    for c in col_name:
        # For each name in the list
        for r in row_name:
            #notes = list() # For percent comparisons between stains
            header = "\n" + c + ": " + r + "\n" 
            filename.write(header)
            
            # This is where analyzing starts, if you change the dataframe it refers to, it will affect everything from here on
            df_QueryCount = df[df.loc[:, c] == r]
            df_QueryPercent = df[df.loc[:, c] == r]
            filename.write("Percentage of Parent\n")
            df_Percent = formatTable(df_QueryPercent, "Percentage of Parent")
            df_Percent.to_csv(filename, sep="\t")
            filename.write("\n")
            filename.write("Counts\n")
            df_Count = formatTable(df_QueryCount, "Count")
            df_Count.to_csv(filename, sep="\t")
#==============================================================================
#             df_samples = set(df_Query["Sample Name"].values)
#             end_Set = set()
#==============================================================================
            
#==============================================================================
#             # Building a set of names for each individual workspace
#             for sn in df_samples:
#                 name = sn.split("_")[0] + "_" + sn.split("_")[1]
#                 end_Set.add(name)
#==============================================================================
                

    #==============================================================================
    #         # Looking for instances of each sample name in the dataframe
    #         for i in end_Set: 
    #             df_Match_Name = df_Query[df_Query["Sample Name"].str.contains(i)]   # Extremely convenient substring match method for pandas
    #             if(len(df_Match_Name) > 1):
    #                 count_Comparison = list(df_Match_Name.loc[:, "Percentage of Parent"])   # Inputting the counts into a list for comparison
    #                 
    #                 # Tricky subtraction, subtracts all the values in the list from each other, if greater than 10% than report
    #                 for count in range(len(count_Comparison)):
    #                     for count_forwards in range(count, len(count_Comparison)):
    #                         diff = math.fabs(float(count_Comparison[count]) - float(count_Comparison[count_forwards]))
    #                         #print(str(float(count_Comparison[count])) + "-" + str(float(count_Comparison[count_forwards])))
    #                         if (diff > 10):
    #                             notes_Stmt = i + ": Stain " + str(count + 1) + " and Stain " + str(count_forwards + 1) + " are greater than 10%! Their difference is: " + ("%.3f" % diff)
    #                             notes.append(notes_Stmt)
    #==============================================================================
        
#==============================================================================
#         for n in notes:
#             filename.write("\n")
#             filename.write(n)
#         filename.write("\n")
#==============================================================================

def formatTable(df, value):
    dfList = list()
    if value == "Percentage of Parent":
        dfFormat = deleteCols(df, ["Percentage of Grandparent", "Count", "ID", "Parent ID"])
    elif value == "Count":
        dfFormat = deleteCols(df, ["Percentage of Grandparent", "Percentage of Parent", "ID", "Parent ID"])
    df_Sample = set(dfFormat["Sample Name"].values)
    sample_Set = set()
    for sn in df_Sample:
        name = sn.split("_")[0] + "_" + sn.split("_")[1]
        sample_Set.add(name)
    sample_Set = list(sample_Set)
    sample_Set.sort()
    for sample in sample_Set:
        dfFormatted = dfFormat[dfFormat["Sample Name"].str.contains(sample)]
        dfFormatted.loc[:,"Sample Name"] = sample
        dfFormatted = pd.pivot_table(dfFormatted, index=["Sample Name", "Parent", "Grandparent"], columns=["Cell Type", "Stain Name"], values=value)
        dfList.append(dfFormatted)
    dfFinal = pd.concat(dfList)
    return dfFinal
    
def deleteCols(df, colList):
    """Deletes all columns in a list from the dataframe."""
    for c in colList:
        del df[c]
    return df
        
        
if __name__ == '__main__':
    
    colQuery = list()
    df = parse_Populations.getGlobalDB()
    
    oneFile = input("All in one file? [y/n] ") 
    cn = input("What column are you looking for? Enter one of these choices [Sample Name, Cell Type, Stain Name, Parent, Count, Percentage of Parent] exactly as printed (enter 0 to quit): ")
    cn = cn.strip("\n").title()
    if (cn == "0"):
        raise SystemExit
    # If one file is no, then keep looping
    if oneFile == "n":
        while(cn != str(0)):
            rn = input("What specific row in that column are you looking for? ")
            rn = rn.strip("\n")
            colQuery.append(cn)
            # Added a loop to continue using user input
            while(rn != "0"):
                pop_file = "Query_Outputs/" + cn.replace(" ", "_") + "_and_" + rn.replace(" ", "_") + "_all_wsps.tsv"
                
                # Makes a folder containing the query outputs
                os.makedirs(os.path.dirname(pop_file), exist_ok=True)    
            
                f = open(pop_file, "w")
                locate_Population(f, colQuery, rn, df)
                f.close()
                
                rn = input("What specific row in that column are you looking for? (Enter 0 to exit) ")
                rn = rn.strip("\n")
                
                if (rn == 0):
                    break
            cn = input("What column are you looking for? Enter one of these choices [Sample Name, Cell Type, Stain Name, Parent, Count, Percentage of Parent] exactly as printed (enter 0 to quit): ")
            colQuery = list()
            rowQuery = list()
    
    # If not one file, have user do all of it in one list
    else:
        while (cn != str(0)):
            colQuery.append(cn)
            rn = input("What specific row[s] in that column are you looking for? (list names separated by semicolons (;), no spaces in between semicolons) ")
            rn = rn.strip("\n")
            cn = input("What column are you looking for? Enter one of these choices [Sample Name, Cell Type, Stain Name, Parent, Count, Percentage of Parent] exactly as printed (enter 0 to quit): ")
        pop_file = "Query_Outputs/all_Queries.tsv"
        # Makes a folder containing the query outputs
        os.makedirs(os.path.dirname(pop_file), exist_ok=True)   
        f = open(pop_file, 'w')
        locate_Population(f, colQuery, rn, df)
        f.close()