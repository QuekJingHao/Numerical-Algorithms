"""
Purpose:    This Python script interfaces with the C++ core numerical program for QM3 project
            This script should be used in a Windos environment

Module:     PC4230 Quantum Mechanics III

Author:     A0183722E Quek Jing Hao

Date:       May 2022
"""
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as mpp
import math
import glob
import matplotlib
import os
import os.path
from matplotlib.pyplot import figure
import shutil
from shutil import copyfile
#----------------------------------------------------------------------------------------------------------------------------
#                                       *** Function declaratons ***
#----------------------------------------------------------------------------------------------------------------------------                                
def Plotting_Final_Results(X_Array, Y_Array, title, Marker, Label, x_label, y_label, save_name):
    figure(figsize=(15, 8), dpi=80)
    mpp.title(title, fontsize = 14)
    for k in range(len(Y_Array)):
        mpp.plot(X_Array[k], Y_Array[k], Marker[k], label = Label[k])
    mpp.xlabel(x_label, fontsize = 14); mpp.ylabel(y_label, fontsize = 14)
    mpp.legend(prop={'size': 12}, loc = "upper right"); mpp.grid(True)
    mpp.savefig(save_name)
    mpp.show()


def Plot_Single(X_Data, Y_Data, Title, x_label, y_label, save_name):
    figure(figsize=(15, 8), dpi=80)
    mpp.title(Title, fontsize = 14)
    mpp.plot(X_Data, Y_Data, '-b')
    mpp.xlabel(x_label, fontsize = 14); mpp.ylabel(y_label, fontsize = 14)
    mpp.grid(True)
    mpp.show()
    mpp.savefig(save_name)

#----------------------------------------------------------------------------------------------------------------------------
#                                *** Move over csv files from Main folder  ***
#----------------------------------------------------------------------------------------------------------------------------                                
Main_path = "C:/Users/Sam Quek/OneDrive/Documents/C++/Projects/QM3 Project Full C++ Implementation/Main/"
Py_path = "C:/Users/Sam Quek/OneDrive/Documents/C++/Projects/QM3 Project Full C++ Implementation/Py/"

filename = ["Psi Initial.csv", "Psi Evolved.csv", "Simulated Transition Probability.csv",
            "Theoretical Transition Probability.csv"]


for i in range(len(filename)):
    shutil.copy(os.path.join(Main_path, filename[i]), Py_path)


#----------------------------------------------------------------------------------------------------------------------------
#                               *** Read data from csv file exported by C++ ***
#----------------------------------------------------------------------------------------------------------------------------                                
csv_files = glob.glob(Py_path + "/*.csv")

df_list = []
for files in csv_files:
    df_list.append(pd.read_csv(files))

# Extract the data from the dataframe list
Qn1_Psi_Evolved = df_list[0]
Qn1_Psi_Initial = df_list[1]

X = np.array(Qn1_Psi_Initial['X'])
Psi_Initial = np.array(Qn1_Psi_Initial['Y'])
Psi_Evolved = np.array(Qn1_Psi_Evolved['Y'])


Qn34_Theoretical_Probability = df_list[2]
Qn34_Transition_Probability = df_list[4]

Time = np.array(Qn34_Theoretical_Probability['X'])
Theoretical =  np.array(Qn34_Theoretical_Probability['Y'])
Simulated = np.array(Qn34_Transition_Probability['Y'])

#----------------------------------------------------------------------------------------------------------------------------
#                                           *** Question 1 ***
#----------------------------------------------------------------------------------------------------------------------------                                
Psi_Data = [Psi_Initial, Psi_Evolved]
Markers = ['-k', '-r']
Labels = ['Initial Eigenstate', 'Final Eigenstate']

Plotting_Final_Results([X, X], Psi_Data, r'Evolved and Initial Eigenstates Plotted Against Position $x$', 
                       Markers, Labels, r'$X$', 'Probability', '1.png')

#----------------------------------------------------------------------------------------------------------------------------
#                                          *** Question 3, 4 ***
#----------------------------------------------------------------------------------------------------------------------------                                
Probabilities = [Simulated, Theoretical]
Markers_34 = ['-b', '-k']
Labels_34 = ['From Simulation', 'From Theoretical Expression']


Plotting_Final_Results([Time, Time], Probabilities, 'First-order Transition Probability Against Time (No RWA)',
                        Markers_34, Labels_34, r'Time, $t$', 'Transition Probability', '2.png')


#Plot_Single(Time, Theoretical, 'First-order Transition Probability Against Time (No RWA)', r'Time, $t$', 'Transition Probability', '3.png')