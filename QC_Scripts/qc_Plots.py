# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 11:52:50 2017

@author: andyp
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def makeGraph(cell_Type, df, dataType):
    """Makes a pyplot comparing stains across all samples of a certain cell type and saves all graphs to a PDF in QC Folder along with corresponding error notes.
    """
    stains = ["Stain 1", "Stain 2", "Stain 3", "Stain 4"]
    cell_Type_Print = cell_Type.replace(" ", "_")
    cell_Type_Print = cell_Type_Print.replace(":", "")
    error_Notes = 'QC_Graphs/GraphNotes_' + cell_Type_Print + "_" + dataType + '.txt'
    source = 'QC_Graphs/StainGraphs_' + cell_Type_Print + "_" + dataType + '.pdf'
    e = open(error_Notes, 'w')
    
    plt.ion()
    with PdfPages(source) as pdf:
        for sc in range(len(stains)):
            sc2 = sc + 1
            if (sc2 < len(stains)):
                for sc3 in range(sc2, len(stains)):
                    fig, ax = plt.subplots(figsize=(10,10))
                    s1 = ax.scatter([], [], s=60)
                    
                    if(dataType == "Percentage of Parent"):
                        ax.set_xlim(-5, 105)
                        ax.set_ylim(-5, 105)
                    else:
                        ax.set_xlim(0, 100000)
                        ax.set_ylim(0, 100000)
                        
                    # Locating the right dataframe for the stains using cell type
                    ct= df.loc[:, "Cell Type"] == cell_Type
                    cell_DF = df[ct]
                    graph_Stain_Points(s1, ax, cell_DF, stains[sc], stains[sc3], dataType, e)
                    fig.canvas.draw()
                    x = [ax.get_xlim()[0], ax.get_xlim()[1]]
                    y = [ax.get_ylim()[0], ax.get_ylim()[1]]
                    
                    # Inputting R-squared stats
                    ax_Data = s1.get_offsets()     
                    
                    if (len(ax_Data) > 0):
                        r_sq = calculate_Linear_RSq(ax_Data)
                        r_sq_str = "R-Squared: " + ("%.3f" % r_sq)
                        ax.text(0, y[1], r_sq_str, verticalalignment='top', horizontalalignment = "left", fontsize=15, color='red')    
                        
                        # Inputting line of slope 1 for reference
                        s_line = ax.plot(x, y, "r--")
                        
                        # Inputting labels for clarity
                        ax.set_title(cell_Type + " " + dataType)
                        ax.set_xlabel(stains[sc]); ax.set_ylabel(stains[sc3])
                        pdf.savefig()
                        plt.close()
                    else:
                        plt.close()
    e.close()

# Function made for returning points for plotting stains
def graph_Stain_Points(plotName, ax, df, x, y, dataType, errorFile):
    """Takes in a Figure object, an Axes object, a dataframe, the name of the stains, and an error file handler. Parses out each x and y point by stain and plots them on the graph"""
    df_Samples = set(df["Sample Name"].values)  # Getting the samples ordered to plot correctly by sample
    sample_Set = set()
    
    for sn in df_Samples:
        name = sn.split("_")[0] + "_" + sn.split("_")[1]
        sample_Set.add(name)
    
    # Iterate through each sample and see if stain points can be plotted 
    for sample in sample_Set:
        df_Match_Sample = df[df["Sample Name"].str.contains(sample)]
        x_Match = df_Match_Sample[df_Match_Sample.loc[:, "Stain Name"] == x]
        y_Match = df_Match_Sample[df_Match_Sample.loc[:, "Stain Name"] == y]
        
        if (len(x_Match) == 0):
            no_Match = sample + " has no " + x + "\n"
            errorFile.write(no_Match)
        elif (len(y_Match) == 0):
            no_Match = sample + " has no " + y + "\n"
            errorFile.write(no_Match)
        else:
            x1 = list(x_Match.loc[:, dataType].values)
            x1_Value = x1[0]  
            y1 = list(y_Match.loc[:, dataType].values)
            y1_Value = y1[0]
            point = [x1_Value, y1_Value]
            arr = plotName.get_offsets()    # get_offsets returns the x, y values for scatter cause PathCollections object
            arr = np.append(arr, point)
            plotName.set_offsets(arr)
            ax.annotate(sample, point)    

def calculate_Linear_RSq(data):
    """Assumes best fit line is y = x and returns r_squared value"""
    xdata = []
    ydata = []
    for points in range(len(data)):
        xdata.append(data[points][0])
        ydata.append(data[points][1])
        
    # Build set of predicted y values
    predicted_y = []
    for x in range(len(xdata)):
        predicted_y.append(xdata[x])
    sum_of_squared = 0
    sum_of_means = 0
    actual_mean = np.mean(ydata)
    for y in range(len(ydata)):
        sum_of_squared += np.power((ydata[y] - predicted_y[y]), 2)
        sum_of_means += np.power((ydata[y] - actual_mean), 2)
    #print("Sum of means: " + str(sum_of_means) + " Sum of squared: " + str(sum_of_squared))
    r_sq = (1 - (sum_of_squared/sum_of_means))
    return r_sq 
            
if __name__ == '__main__':

    graph_file = "QC_Graphs/"
    
    # Makes a folder containing the query outputs
    os.makedirs(os.path.dirname(graph_file), exist_ok=True)
    df = pd.read_csv("populationDB.csv")
    
    graph_List = ["CD45+", "CD3+ all", "Q1: CD8a- , CD4+"]
    for graphs in graph_List:     
        makeGraph(graphs, df, "Percentage of Parent")
    makeGraph("CD45+", df, "Count")