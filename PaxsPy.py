# -*- coding: utf-8 -*-
'''
This script reads ICP-MS output and calculates Pa activities, which can then be plugged into the 'master' spreadsheet
Yuxin Zhou
yzhou@ldeo.columbia.edu
'''
import numpy as np
import numpy.ma as ma
from scipy import stats # for linear regression
from Tkinter import Tk
from tkFileDialog import askopenfilenames, asksaveasfilename
import sys
import pandas as pd

spike_answer = str(raw_input("Are you using 2006-2 UTh spike and 2017-1a Pa spike? If not, click no and search \'MixedPa' in script and change its values. [y] or n:") or 'y')
if spike_answer == 'n':
    sys.exit()
figure_answer = str(raw_input("Do you want to inspect ICPMS raw output in figures?[y] or n:") or 'y')

## check OS
#if platform.system() == 'Windows':
#    spike_answer = ctypes.windll.user32.MessageBoxA(0, "Are you using 2006-2 UTh spike and 2017-1a Pa spike? If not, click no and search \'MixedPa' in script and change its values", "Spike?", 4)
#    if spike_answer == 7:
#        sys.exit()
#elif platform.system() == 'Darwin':
#    window = Tk()
#    window.wm_withdraw()
#    tkMessageBox.showinfo(title="Spike?", message="Are you using 2006-2 UTh spike and 2017-1a Pa spike? If not, click no and search \'MixedPa' in script and change its values")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
file_names = askopenfilenames(title="Select all the ICPMS output files and a \'sample_info' excel file") # show an "Open" dialog box and return the path to the selected file

def return_five_point_avg(file_name):
    # start reading from row 12, which are name/unit/blank
#    txt_handle = np.genfromtxt(file_name, delimiter='\t', skip_header=12)
    txt_handle = np.genfromtxt(file_name, delimiter='\t')
    if np.all(np.isnan(txt_handle[0])):# if first line is empty
        txt_handle = txt_handle[12:]
    else:
        txt_handle = txt_handle[7:]
    if figure_answer == 'y':
        txt_handle_r = np.transpose(txt_handle)
        data_df = pd.DataFrame(data=txt_handle_r[1:-1,:])
        data_df.plot(sharex=True,title=file_name)
        plt.ion()
        plt.show()
    # get rid of first column (mass) and last column (nan)
    txt_handle = txt_handle[:,1:-1]
    # If not blank, check and remove outliers
    if 'Blank' not in file_name and 'blank' not in file_name:
        txt_handle = reject_outliers(txt_handle)
    # average accros five points
    five_point_avg = ma.mean(txt_handle.reshape(len(txt_handle)/5, 5, -1),axis=1)
    # A second check for outliers after the five point average, except when the file is Blank
    if 'Blank' not in file_name and 'blank' not in file_name:
        print file_name + " # outliers: " + str(ma.count_masked(five_point_avg))
        return reject_outliers(five_point_avg)
    else:
        return five_point_avg

def reject_outliers(data, m = 2.):
    # from https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
#    # median approach
#    d = np.abs(data - np.median(data, axis=1)[:,None])
#    mdev = np.median(d,axis=1)
#    s = d/mdev[:,None]# if mdev!=0 else 0.
    # avg approach
    d = np.abs(data - np.mean(data,axis=1)[:,None])
    s = d/np.mean(data,axis=1)[:,None]
    return ma.array(data,mask=np.greater(s,m))

#%% process blanks and stds. Calculate tailcrxn slope and intercept
names = [name for name in file_names if 'Blank' in name or 'blank' in name or
     'Th_std' in name]
if not names:
    raise RuntimeError('No blank or std files found!')

# set up lists for tail corrections
# the three lists are for 231, 232, 233
Th_std_tailCrxn = [[],[],[]]
blank_Th_tailCrxn = [[],[],[]]

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    two_hundred_run_avg = np.mean(five_point_avg, axis=1)
    if 'Blank' in file_name or 'blank' in file_name:
        blank_Th_tailCrxn[0].append(two_hundred_run_avg[0]) # Th231
        blank_Th_tailCrxn[1].append(two_hundred_run_avg[1]) # Th232
        blank_Th_tailCrxn[2].append(two_hundred_run_avg[2]) # Th233
    if 'Th_std' in file_name:
        Th_std_tailCrxn[0].append(two_hundred_run_avg[0]) # Th231
        Th_std_tailCrxn[1].append(two_hundred_run_avg[1]) # Th232
        Th_std_tailCrxn[2].append(two_hundred_run_avg[2]) # Th233
        
# now do TailCrxn before processing UTh data file
# two arrays to store intercepts and slopes in the sequence of
# 236 vs 238, 234 vs 238, 229 vs 232, 230 vs 232, 234 vs 232
intercepts_tailCrxn = np.zeros(2)
slopes_tailCrxn = np.zeros(2)
correlations_tailCrxn = np.zeros(2)

Th231_TailCrxn = np.concatenate((Th_std_tailCrxn[0], blank_Th_tailCrxn[0]))
Th232_TailCrxn = np.concatenate((Th_std_tailCrxn[1], blank_Th_tailCrxn[1]))
slopes_tailCrxn[0], intercepts_tailCrxn[0], correlations_tailCrxn[0] = stats.linregress(Th232_TailCrxn, Th231_TailCrxn)[:3]
Th233_TailCrxn = np.concatenate((Th_std_tailCrxn[2], blank_Th_tailCrxn[2]))
slopes_tailCrxn[1], intercepts_tailCrxn[1], correlations_tailCrxn[1] = stats.linregress(Th232_TailCrxn, Th233_TailCrxn)[:3]

#%% SRM
names = [name for name in file_names if 'SRM' in name]
if not names:
    raise RuntimeError('No SRM files found!')

# set up lists to store the 3 SRM_a
SRM_238235_avg = []
SRM_238235_std = []
SRM_238235_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_238235_avg = np.mean(five_point_avg[2
                                                    ,:]/five_point_avg[1,:])
    SRM_238235_avg.append(two_hundred_run_238235_avg)
    two_hundred_run_238235_std = np.std(five_point_avg[2
                                                   ,:]/five_point_avg[1,:])/np.sqrt(five_point_avg.shape[1])
    SRM_238235_std.append(two_hundred_run_238235_std)
    two_hundred_run_238235_RSD = two_hundred_run_238235_std/two_hundred_run_238235_avg
    SRM_238235_RSD.append(two_hundred_run_238235_RSD)
    
#%% MixPa
names = [name for name in file_names if 'zmix' in name]
if not names:
    raise RuntimeError('No zmix files found!')

# set up lists to store the 3 MixPa
MixPa_233231_avg = []
MixPa_233231_std = []
MixPa_233231_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)

    # tail correction
    five_point_avg[0,:] -= slopes_tailCrxn[0] * five_point_avg[1,:] + intercepts_tailCrxn[0]
    five_point_avg[2,:] -= slopes_tailCrxn[1] * five_point_avg[1,:] + intercepts_tailCrxn[1]
    
    two_hundred_run_233231_avg = np.mean(five_point_avg[2
                                                    ,:]/five_point_avg[0,:])
    MixPa_233231_avg.append(two_hundred_run_233231_avg)
    two_hundred_run_233231_std = np.std(five_point_avg[2
                                                   ,:]/five_point_avg[0,:])/np.sqrt(five_point_avg.shape[1])
    MixPa_233231_std.append(two_hundred_run_233231_std)
    two_hundred_run_233231_RSD = two_hundred_run_233231_std/two_hundred_run_233231_avg
    MixPa_233231_RSD.append(two_hundred_run_233231_RSD)
    
#%% MixPa
names = [name for name in file_names if 'nochem_mixPa' in name]

if len(names) >1:
    raise RuntimeError('Multiple nochem_mixPa!')
if not names:
    nochem_mixPa_flag = False
else:
    nochem_mixPa_flag = True
    
if nochem_mixPa_flag:
    file_name = names[0]
    five_point_avg = return_five_point_avg(file_name)
    
    # tail correction
    five_point_avg[0,:] -= slopes_tailCrxn[0] * five_point_avg[1,:] + intercepts_tailCrxn[0]
    five_point_avg[2,:] -= slopes_tailCrxn[1] * five_point_avg[1,:] + intercepts_tailCrxn[1]
    
    nochem_mixPa_233231_avg = np.mean(five_point_avg[2
                                                    ,:]/five_point_avg[0,:])
    nochem_mixPa_233231_std = np.std(five_point_avg[2
                                                   ,:]/five_point_avg[0,:])/np.sqrt(five_point_avg.shape[1])
    nochem_mixPa_233231_RSD = two_hundred_run_233231_std/two_hundred_run_233231_avg

#%% sample results

# if this is UTh data file
names = [name for name in file_names if '_Pa.txt' in name]
if not names:
    raise RuntimeError('No Pa files found!')
names.sort()

# set up the 2d array as in master spreadsheet
# Columns: 231/233_avg	231/233_RSD
# Rows: Pa1-num_samples
num_samples = len(names)
master = np.zeros((num_samples,2))

for i, file_name in enumerate(names):
    five_point_avg = return_five_point_avg(file_name)
    
    # tail correction
    five_point_avg[0,:] -= slopes_tailCrxn[0] * five_point_avg[1,:] + intercepts_tailCrxn[0]
    five_point_avg[2,:] -= slopes_tailCrxn[1] * five_point_avg[1,:] + intercepts_tailCrxn[1]
    
    master[i,0] = np.mean(five_point_avg[0,:]/five_point_avg[2,:])
    two_hundred_run_231233_std = np.std(five_point_avg[0
                                                   ,:]/five_point_avg[2,:])/np.sqrt(five_point_avg.shape[1])
    master[i,1] = two_hundred_run_231233_std/master[i,0]
    
#%% ez reduction
    
# sample info. Exclude $ in file name in case that file is open
names = [name for name in file_names if 'info' in name and '$' not in name]
if len(names)>1:
    raise RuntimeError('More than one sample info file')
if not names:
    raise RuntimeError('Sample info file cannot be found. The file must have \'info\' in file name')
sample_info_type = ''
if names[0][-3:] == 'txt':
    sample_info_type = 'txt'
    try:
        sample_info = np.genfromtxt(names[0], delimiter='\t',dtype=None,skip_header=1)
    except ValueError:
        raise ValueError('In reading file ' + names[0] + ', value error!')
elif names[0][-4:] == 'xlsx':
    sample_info_type = 'xlsx'
    try:
        sample_info = pd.read_excel(names[0],header=None,skiprows=1)
    except ValueError:
        raise ValueError('In reading file ' + names[0] + ', value error!')
else:
    raise ValueError(names[0] + ' is not either a txt or excel file and cannot be processed')

## MixSpikeMassBias
# mass bias
SRM = np.mean(SRM_238235_avg)
SRM_RSD = np.sqrt((np.sum((np.array(SRM_238235_avg) * np.array(SRM_238235_RSD))**2)))/3/SRM
accepted_238235 = 137.55
accepted_238235_RSD = 0.50*0.01
mass_bias_per_amu = (SRM/accepted_238235-1)/3
mass_bias_per_amu_RSD = np.sqrt((SRM_RSD**2+accepted_238235_RSD**2))

# MixedPa
avg_233Pa_231Pa = np.mean(MixPa_233231_avg)
avg_233Pa_231Pa_RSD = np.sqrt(np.sum((np.array(MixPa_233231_avg) * np.array(MixPa_233231_RSD))**2))/3/avg_233Pa_231Pa
Pa231_cncn_pg_g = 38
Pa231_soln_mass_in_spike_g = 0.2088
Pa233_soln_mass_in_spike_g = 5.3203
Pa233_mass_in_spike_pg_g = avg_233Pa_231Pa*Pa231_cncn_pg_g*Pa231_soln_mass_in_spike_g/Pa233_soln_mass_in_spike_g
decay_days = int(input('Enter the number of days from Pa233 spike preparation to Pa clean up column: '))

## write avg_233Pa_231Pa to txt file for Marty
#with open("avg_233_231Pa.txt", "w") as text_file:
#    text_file.write("Purchase Amount: {}".format(avg_233Pa_231Pa))
#%%
##UnSpike
Pa231_Pa233_mass_bias_crtd = (1+(-2*mass_bias_per_amu))*master[:,0]
Pa231_Pa233_mass_bias_crtd_RSD = np.sqrt(master[:,1]**2+(2*mass_bias_per_amu_RSD)**2)
if sample_info_type == 'txt':
    Pa233_pg = Pa233_mass_in_spike_pg_g * (sample_info['f4']/1000)
elif sample_info_type == 'xlsx':
    Pa233_pg = Pa233_mass_in_spike_pg_g * (sample_info[4]/1000)
Pa231_pg = Pa233_pg * Pa231_Pa233_mass_bias_crtd * 231 / 233
Pa231_pg_RSD = np.sqrt(Pa231_Pa233_mass_bias_crtd_RSD**2 + avg_233Pa_231Pa_RSD**2)
if sample_info_type == 'txt':
    if not (sample_info['f0']=='BLANK').any():
        raise RuntimeError('Cannot determine from sample name in sample info which sample is blank. Name it BLANK')
    blank_index = np.argwhere(sample_info['f0']=='BLANK')
    Pa231_pg_minus_blank = Pa231_pg - Pa231_pg[blank_index]
elif sample_info_type == 'xlsx':
    if not (sample_info[0]=='BLANK').any():
        raise RuntimeError('Cannot determine from sample name in sample info which sample is blank. Name it BLANK')
    blank_index = sample_info[0][sample_info[0]=='BLANK'].index[0]
    Pa231_pg_minus_blank = Pa231_pg - Pa231_pg[blank_index]
atmoic_mass_g_mol = 231
Avogadro=6e23
lambda_1_year = 2.12e-5
lambda_1_min = lambda_1_year/525960
Pa231_atoms = Pa231_pg_minus_blank*1e-12/atmoic_mass_g_mol*Avogadro
Pa231_dpm = Pa231_atoms * lambda_1_min
if sample_info_type == 'txt':
    sed_mass_g = sample_info['f2'] / 1000
elif sample_info_type == 'xlsx':
    sed_mass_g = sample_info[2] / 1000
Pa231_dpm_g = Pa231_dpm / sed_mass_g
Pa231_dpm_g_2_sigma = Pa231_dpm_g * Pa231_pg_RSD * 2
export = np.zeros((num_samples,2))
export[:,0] = Pa231_dpm_g
export[:,1] = Pa231_dpm_g_2_sigma

# since numpy array can't have both string and float, converting to pandas dataframe and add sample name as the first column in export
export_data_df = pd.DataFrame(data=export,index=np.arange(num_samples),columns=['231Pa dpm/g','231 Pa (dpm/g) 2 sigma'])
if sample_info_type == 'txt':
    sample_name_df = pd.DataFrame({'Sample name':sample_info['f0']})
elif sample_info_type == 'xlsx':
    sample_name_df = pd.DataFrame({'Sample name':sample_info[0]})
if nochem_mixPa_flag:
    avg_233Pa_231Pa_df = pd.DataFrame({'avg_233Pa_231Pa for Marty':[decay_days,avg_233Pa_231Pa,nochem_mixPa_233231_avg]},index=[0,1,2])
else:
    avg_233Pa_231Pa_df = pd.DataFrame({'avg_233Pa_231Pa for Marty':[decay_days,avg_233Pa_231Pa]},index=[0,1])
export_df = pd.concat([sample_name_df,export_data_df,avg_233Pa_231Pa_df],axis=1)
#%% save to csv
output_file_name = asksaveasfilename(title='Save the output file as')
if 'xlsx' not in output_file_name:
    output_file_name = output_file_name + '.xlsx'
export_df.to_excel(output_file_name)