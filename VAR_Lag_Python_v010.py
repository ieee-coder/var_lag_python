#******************************************************************************
#   Name of the program: OptLag-HJC
#   Title of the program:OptLag-HJC: Python Module to Estimate and Determine 
#   the Optimal Lag Order in a VAR Model by Minimizing an Information Criterion.
#   Version: 1.0
# 
#   This program is released under GNU General Public License, version 3.
# 
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#   This software has been developed by Dr. Alan Mustafa under supervision of 
#   Prof. Abdulnasser Hatemi-J (Hatemi-J, 2012).
#Contacts:
#    - Prof. Abdulnasser Hatemi-J: AHatemi@uaeu.ac.ae
#    - Dr. Alan Mustafa: Alan.Mustafa@ieee.org
# 
#    Date: October 2020
#
#    © 2020 Prof. Abdulnasser Hatemi-J and Dr. Alan Mustafa
# 
#******************************************************************************
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
#import scipy
#from scipy.stats import chi2
#######################################
#from tkinter import *
import tkinter as tk
#from tkinter import Menu
from tkinter.filedialog import askopenfilename
#from tkinter import messagebox
#from glob import glob
from datetime import datetime
import textwrap
###############################################################################
#                         Start of GUI
###############################################################################


#========================= functions =========================================
class create_window_menu_UI(tk.Frame):

    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.master = master
        #============= Application Title --------------------------------------
        #self.lblTitle = tk.Label(self, text="Granger Causality Test", font=("Helvetica", 16))
        self.lblTitle = tk.Label(self, text="OptLag-HJC: Determining The optimal lag order of a VAR model", font=("Helvetica", 16))
        self.lblTitle.grid(row=0, column=0, columnspan=10, sticky="NWNESWSE")        
        
        #============= Drawing Horizontal Line --------------------------------
        hr = tk.Frame(self,height=2,width=650,bg="green")
        hr.grid(row=1, column=0, columnspan=10)
        
        #============= Dataset File Selection ---------------------------------
        self.lblDatasetFile = tk.Label(self, text="Dataset file:", font=("Helvetica", 12))
        self.lblDatasetFile.grid(row=2, column=0, columnspan=3, sticky="E")
        
        var_DatasetFile = tk.StringVar()
        self.tbx_DatasetFile = tk.Entry(self, textvariable=var_DatasetFile, font=("Helvetica", 12), state="disabled", justify="right")
        self.tbx_DatasetFile.grid(row=2, column=3, columnspan=3, sticky="W")
        
        self.btnSelectDataFile = tk.Button(self, text="Select the Dataset File", command=lambda: var_DatasetFile.set(os.path.split(askopenfilename())[1]), font=("Helvetica", 12))
        self.btnSelectDataFile.grid(row=2, column=6, columnspan=4, sticky="W")

        #============= Maximum Lag Number -------------------------------------
        self.lblMaxLag = tk.Label(self, text="Maximum Lag Number:", font=("Helvetica", 12))
        self.lblMaxLag.grid(row=3, column=0, columnspan=3, sticky="E")        
        
        var_maxLag = tk.StringVar()
        self.tbx_maxLag = tk.Entry(self, textvariable=var_maxLag, font=("Helvetica", 12))
        var_maxLag.set(5)
        self.tbx_maxLag.grid(row=3, column=3, columnspan=4, sticky="W")

        #============= Calculate the Lag Order --------------------------------
        self.btnCalcOptLag = tk.Button(self, text="Determine Optimal Lag", command=lambda: calc_optlag(int(self.tbx_maxLag.get()), var_DatasetFile,self), font=("Helvetica", 12))
        self.btnCalcOptLag.grid(row=3, column=6, columnspan=4, sticky="W")
       
        #============= Adding a Blank lable -----------------------------------
        #self.lblBlank = tk.Label(self, text=" ", font=("Helvetica", 8))
        #self.lblBlank.grid(row=3, column=9, columnspan=1, sticky="W")
 
       #============= Adding a Blank Row --------------------------------------
        self.lblBlank = tk.Label(self, text=" ", font=("Helvetica", 8))
        self.lblBlank.grid(row=4, column=0, columnspan=1, sticky="W")
        
       #============= Adding References ---------------------------------------
        #txtRef = "- Hatemi-J, A. (2003) A new method to determine optimal lag order in stable and unstable VAR models, Applied Economics Letters, Vol 10, 135-137. \n - Hatemi-J, A. (2008) Forecasting properties of a new method to determine optimal lag order in stable and unstable VAR models, Applied Economics Letters, Vol 15, 239-243."
        #txtRef = '''- Hatemi-J, A. (2003) A new method to determine optimal lag order in stable and unstable 
        #VAR models, Applied Economics Letters, Vol 10, 135-137. \n- Hatemi-J, A. (2008) Forecasting properties 
        #of a new method to determine optimal lag order in stable and unstable VAR models, Applied Economics Letters, 
        #Vol 15, 239-243.'''

        #self.lblRef = tk.Label(self, text=txtRef, font=("Helvetica", 8), wraplength=600, justify="left")
        #self.lblRef.grid(row=8, column=0, columnspan=9, sticky="W")
        
        #============= Close Button -------------------------------------------
        self.btnExit = tk.Button(self, text="Close", command=self.master.destroy, font=("Helvetica", 12))
        self.btnExit.grid(row=5, column=7, sticky="E")
       
        #============= Printing Credits ---------------------------------------
        #self.btnCredit = tk.Button(self, text="Credit", command=cmdMsgBoxCredit, font=("Helvetica", 12))
        #self.btnCredit.grid(row=6, column=3, sticky="W")
       
        #============= BlankSpace              ---------------------------------
        #self.lblblnkSpace = tk.Label(self, text="", font=("Helvetica", 8))
        #self.lblblnkSpace.grid(row=6, column=0, columnspan=10, sticky="E")

        #============= Drawing Horizontal Line --------------------------------
        hr = tk.Frame(self,height=2,width=650,bg="green")
        hr.grid(row=7, column=0, columnspan=10)
        
        #============= Output Message         ---------------------------------
        self.lblOutputMsg = tk.Label(self, text="", font=("Helvetica", 12), anchor="w")
        self.lblOutputMsg.grid(row=8, column=0, columnspan=9, sticky="w")
                 
        #============= Output EndNote         ---------------------------------
        #txtEndNote = '''Note:\n\nA copy of the same report and a chart has been added to the same folder as the program reside in.'''

        #self.lblEndNote = tk.Label(self, text=lblEndNote, font=("Helvetica", 8), wraplength=600, justify="left")
        self.lblEndNote = tk.Label(self, text="", font=("Helvetica", 8), wraplength=600, justify="left")
        self.lblEndNote.grid(row=9, column=0, columnspan=9, sticky="W")
        

###############################################################################
#                           End of GUI                                        #
###############################################################################

###############################################################################
#             Start of Calculations: Asymmetric Causality Test                #
###############################################################################

#def cmdMsgBoxCredit():
#    messagebox.showinfo('Credit', 'Some credits!')

def y_calc(dataset, hjc_p_value):
    return dataset.iloc[hjc_p_value:, :].values;
    
def calc_z_p_all(dataset,p_value):
    all_ones = np.ones(dataset.shape[0]-p_value)
    all_ones = all_ones.reshape((len(all_ones), 1))
    
    Z_All_ = all_ones

    for p in range(1,p_value+1):
        Z_All_ = np.concatenate((Z_All_, dataset.iloc[p_value - p:-p, :].values), axis=1)
    return Z_All_;

def z_all_calc(dataset,hjc_p_value):
    return calc_z_p_all(dataset,hjc_p_value);
    
def inv_z_z_T_calc(z_all, hjc_p_value, p_value):
    return np.linalg.inv(np.matmul(z_all.transpose(),z_all));
    
def B_hat_calc(y,z_all,inv_z_z_T, p, p_value):
    return np.matmul(np.matmul(y.transpose(),z_all),inv_z_z_T);
    

def get_hjc(dataset, hjc_p_value, p_value):

    z_all = z_all_calc(dataset,hjc_p_value)
    
    y = y_calc(dataset, hjc_p_value)
    
    inv_z_z_T = inv_z_z_T_calc(z_all, hjc_p_value, p_value)
    
    B = B_hat_calc(y,z_all,inv_z_z_T, hjc_p_value, p_value)
    D = y.transpose() - np.matmul(B, z_all.transpose())
    
    omega = (np.matmul(D, D.transpose()))/(T - k * hjc_p_value - 1)
    omega_all = np.linalg.det(omega)
    
    
    hjc = np.log(omega_all) + (hjc_p_value * (((k**2 * np.log(T)) + (2 * k**2 * np.log(np.log(T)))))/(2 * T))
    
    return hjc, B, inv_z_z_T, omega

################################################
def calc_optlag(p_value,dataset_file,self):
    dataset = pd.read_csv(dataset_file.get())
    
    #=============== creating Output Report File ==============================
    now = datetime.now()
#    print("now =", now)
#    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
#    print("date and time =", dt_string)    
    dt_string2 = now.strftime("%Y%m%d-%H%M%S")
    f_name = dataset_file.get()
    
    fn = f_name.split(".")
#    print("File name : ", fn[0])
#    print("File extension : ", fn[1])
#    fileName = 'optimal_lag_order_rprt_for_' + dataset_file.get() + '.txt'
    fileName = 'optimal_lag_order_rprt_for_' + fn[0] + "-" + dt_string2 + '.txt'
    output_rprt_file = open(fileName,'w')
    output_rprt_file.write('###############################################################################\n')
    output_rprt_file.write('#                                                                             #\n')
    output_rprt_file.write('#                           OPTIMAL LAG ORDER REPORT                          #\n')
    output_rprt_file.write('#                             ' + fn[0] + '                  \n')
    output_rprt_file.write('#                                                                             #\n')
    output_rprt_file.write('###############################################################################\n')
    output_rprt_file.close
    
    output_rprt_file = open(fileName,'a')
    

    
    
    global B
    B = []
    
    global B_hat_all
    B_hat_all = []
    
    global inv_z_z_T
    inv_z_z_T = []
    
    global inv_z_z_T_all
    inv_z_z_T_all = []
    
    global omega
    omega = []
    
    global O_all
    O_all = []
    
    global T
    global k
    T = dataset.shape[0]
    k = dataset.shape[1]
    
    
    hjc = []
    lag_values = []
#    B_flatten = []
#    B_flatten_all = []
#    B_elements = []
    
    B_elements_appened = []
    
    B_elements_appened.append([])
    
    for hjc_p_value in range(1,p_value+1):
    
        hjc.append(get_hjc(dataset, hjc_p_value, p_value)[0])
    
        B = get_hjc(dataset, hjc_p_value, p_value)[1]
        B_hat_all.append(B)
        
        #inv_z_z_T
        inv_z_z_T = get_hjc(dataset, hjc_p_value, p_value)[2]
        inv_z_z_T_all.append(inv_z_z_T)
        
        omega = get_hjc(dataset, hjc_p_value, p_value)[3]
        O_all.append(omega)
        
    
        lag_values.append(hjc_p_value)
    
    global lag_min
    lag_min = np.argmin(hjc) + 1
    
    hjc_min = min(hjc)
    hjc_max = max(hjc)
    
    annotation_coord_x = lag_min + 0.5
    annotation_coord_y = (hjc_min + hjc_max)/2
    
#    print("=========================================================================================")
    print("")
    print("The Optimal Lag Order = ", lag_min, "with the value of ", hjc_min)
    print("")
    print("=========================================================================================")
    print('Lag values: ', lag_values)
    print('HJC: ', hjc)

    txtRef = '''- Hatemi-J, A. (2003) A new method to determine optimal lag order in stable and unstable VAR models, Applied Economics Letters, Vol 10, 135-137.\n\n
    - Hatemi-J, A. (2008) Forecasting properties of a new method to determine optimal lag order in stable and unstable VAR models, Applied Economics Letters, Vol 15, 239-243.'''
#    print(txtRef)
    print("=========================================================================================")
    
    
    output_rprt_file.write('\n\n')
#    output_rprt_file.write('=============================================================\n\n')
    output_rprt_file.write("The Optimal Lag Order is [%d] with the value of [%.5f].\n\n" % (lag_min, hjc_min))
    output_rprt_file.write('\n')
    output_rprt_file.write('===============================================================================\n')
    output_rprt_file.write('\n')
    output_rprt_file.write('References:\n')
    output_rprt_file.write(textwrap.fill(txtRef,75))
    

    self.lblOutputMsg["text"] = "The Optimal Lag Order is: %d\n with the value of: %.5f\n\n" % (lag_min, hjc_min)
    txtEndNote = '''Note:\nA copy of the same report and a chart has been added to the same folder as the program resides in.'''
    self.lblEndNote["text"] = "%s\n" % (txtEndNote)
    
    
    plt.plot(lag_values, hjc, color = 'blue')
    plt.axis([lag_values[0]-1, lag_values[-1]+1, hjc_min-0.05, hjc_max+0.05])
    plt.suptitle('The Curve of HJC Values for Maximum Lag of %i - %s' %(p_value,fn[0]), fontsize=12, fontweight='bold')
    plt.title('The optimal lag order is [%d] with the value of [%.5f]' %(lag_min,hjc_min), fontsize = 12)

    plt.xlabel('Lag Number')
    plt.ylabel('HJC Value')
    plt.legend(['Information Values'])
    plt.annotate('Optimal lag = %i' % lag_min, xy = (lag_min,hjc_min),
                 xycoords='data',
                 xytext=(annotation_coord_x,annotation_coord_y),
                 textcoords='data', 
                 arrowprops=dict(facecolor='orange', arrowstyle='simple', connectionstyle='arc3, rad=0.3'),
                 size=15, va="center", ha="center",
                 )
    
    HJC_eq="HJC=  ln⁡(|Ω̂p|)+p((k^2  ln⁡T+2k^2  ln⁡(ln⁡T ))/2T),p=0…lag_order_value"
    plt.figtext(0.5, -0.11, HJC_eq, wrap=True, horizontalalignment='center', fontsize=12)

    #plt.savefig('MinLag.pdf')

#    plt.savefig('Optimal_Lag_Order.png')
    plt.savefig('Optimal_Lag_Order_chart_for_' + fn[0] + "-" + dt_string2 + '.png')
    plt.show()
    

    output_rprt_file.close()
    return;

###############################################################################
#              End of Calculations: Calculating Optimal Lag Order             #
###############################################################################

#In our main function, create the GUI and pass it to our App class
def main():
#    window= Tk()
    window= tk.Tk()
    #window.title("Causality Tests 0.1 (CT 0.1)")
    window.title("Optimal lag order (OptLag-HJC 1.0)")
    window.geometry('650x400')

    create_window_menu_UI(window).grid(row=0, column=1, columnspan=1)
    window.mainloop()

#Run the main function
if __name__ == "__main__":
    main()

