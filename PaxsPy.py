# -*- coding: utf-8 -*-
'''
This script reads ICP-MS output and calculates Pa activities, which can then be plugged into the 'master' spreadsheet
Yuxin Zhou
yzhou@ldeo.columbia.edu
'''
import numpy as np
import numpy.ma as ma
from scipy import stats # for linear regression
import tkinter as tk
from tkinter import filedialog
import sys
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pprint

spike_answer = str(input("Are you using 2006-2 UTh spike and 2022-1a Pa spike? If not, click no and search \'MixedPa' in script and change its values. [y] or n:") or 'y')
if spike_answer == 'n':
    sys.exit()
figure_answer = str(input("Do you want to inspect ICPMS raw output in figures?[y] or n:") or 'y')

root = tk.Tk()
root.withdraw() # we don't want a full GUI, so keep the root window from appearing
file_names = filedialog.askopenfilenames(title="Select all the ICPMS output files and a \'sample_info' excel file") # show an "Open" dialog box and return the path to the selected file

def return_five_point_avg(file_name):
    # read txt as csv, using tab as separator
    txt_handle = pd.read_csv(file_name,sep='\t',header=None)
    txt_handle.dropna(how='all',axis='index',inplace=True) # drop the rows where all elements are missing
    txt_handle.dropna(how='all',axis='columns',inplace=True) # drop the columns where all elements are missing
    txt_handle.reset_index(drop=True,inplace=True) # index at this point doesn't start with 0,1,2 etc because some rows were deleted. This will reset index, and drop=True will prevent a new column named "index" be created
    txt_handle.drop([0,1,2],inplace=True) # remove the first three rows 
    txt_handle.reset_index(drop=True,inplace=True)
    txt_handle = txt_handle.astype(float)
    if figure_answer == 'y':
        txt_handle_r = txt_handle.transpose() # create a transposed version
        txt_handle_r.columns = txt_handle_r.iloc[0] # name the columns with the first row (mass)
        txt_handle_r.drop(txt_handle_r.index[0],inplace=True) # drop the first row (mass)
        txt_handle_r.plot(sharex=True,title=file_name)
        plt.savefig(file_name+'.png')
    txt_handle.set_index(txt_handle[0],inplace=True) # set index as mass
    txt_handle.drop(columns=[0],inplace=True) # drop the mass column
    txt_handle = reject_outliers(txt_handle)
    # average accros multiple masses of the same element
    txt_handle.set_index(np.floor(txt_handle.index),inplace=True)
    five_point_avg = txt_handle.groupby(txt_handle.index).mean()
    five_point_avg_2nd_outlier_detection = reject_outliers(five_point_avg)
    masekd_array = ma.masked_invalid(five_point_avg_2nd_outlier_detection.values)
    print(file_name + ' # outliers: ' + str(np.count_nonzero(~np.isnan(masekd_array))))
    return masekd_array
    
def reject_outliers(data, m = 2.):
    '''
    from https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
    data is expected to be a pandas dataframe, where each row is a number of measurements on the same mass
    '''
    d = data.subtract(data.median(axis=1),axis='index').abs()
    mdev = d.median(axis=1)
    s = d.divide(mdev,axis='index')
    return data.mask(s>m)

#%% process blanks and stds. Calculate tailcrxn slope and intercept
names = [name for name in file_names if ('Blank' in name or 'blank' in name or 'BLANK' in name or 'Th_std' in name) and 'SRM' not in name]
if not names:
    raise RuntimeError('No blank or std files found!')
print("Identified the following files as either blank or Th_std:")
pprint.pprint(names)
print('\n')
# set up lists for tail corrections
# the three lists are for 231, 232, 233
Th_std_tailCrxn = [[],[],[]]
blank_Th_tailCrxn = [[],[],[]]

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    two_hundred_run_avg = ma.mean(five_point_avg, axis=1)
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

#%% SRM blank
names = [name for name in file_names if ('SRM' in name and 'blank' in name) or ('SRM' in name and 'BLANK' in name) or ('SRM' in name and 'Blank' in name)]
if not names:
    SRM_blank_flag = False
else:
    SRM_blank_flag = True
    print("Identified the following files as SRM blanks:")
    pprint.pprint(names)
    print('\n')
# set up lists to store the 3 SRM_a
SRM_238_blank_avg = []
SRM_235_blank_avg = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_238_avg = ma.mean(five_point_avg[2,:])
    two_hundred_run_235_avg = ma.mean(five_point_avg[1,:])
    SRM_238_blank_avg.append(two_hundred_run_238_avg)
    SRM_235_blank_avg.append(two_hundred_run_235_avg)

#%% SRM
names = [name for name in file_names if 'SRM' in name and 'blank' not in name and 'Blank' not in name and 'BLANK' not in name]
if not names:
    raise RuntimeError('No SRM files found!')
print("Identified the following files as SRM:")
pprint.pprint(names)
print('\n')
# set up lists to store the 3 SRM_a
SRM_238_avg = []
SRM_235_avg = []
SRM_238235_std = []
SRM_238235_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_238_avg = ma.mean(five_point_avg[2,:])
    two_hundred_run_235_avg = ma.mean(five_point_avg[1,:])
    SRM_238_avg.append(two_hundred_run_238_avg)
    SRM_235_avg.append(two_hundred_run_235_avg)
    two_hundred_run_238235_avg = ma.mean(five_point_avg[2
                                                    ,:]/five_point_avg[1,:])
    two_hundred_run_238235_std = ma.std(five_point_avg[2
                                                   ,:]/five_point_avg[1,:])/ma.sqrt(five_point_avg.shape[1])
    SRM_238235_std.append(two_hundred_run_238235_std)
    two_hundred_run_238235_RSD = two_hundred_run_238235_std/two_hundred_run_238235_avg
    SRM_238235_RSD.append(two_hundred_run_238235_RSD)
    
if SRM_blank_flag:
    SRM_238235_avg = (SRM_238_avg - ma.mean(SRM_238_blank_avg)) / (SRM_235_avg - ma.mean(SRM_235_blank_avg))
else:
    SRM_238235_avg = ma.array(SRM_238_avg) / ma.array(SRM_235_avg)
#%% MixPa
names = [name for name in file_names if 'zmix' in name]
if not names:
    raise RuntimeError('No zmix files found!')
print("Identified the following files as zmix:")
pprint.pprint(names)
print('\n')
# set up lists to store the 3 MixPa
MixPa_233231_avg = []
MixPa_233231_std = []
MixPa_233231_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)

    # tail correction
    five_point_avg[0,:] -= slopes_tailCrxn[0] * five_point_avg[1,:] + intercepts_tailCrxn[0]
    five_point_avg[2,:] -= slopes_tailCrxn[1] * five_point_avg[1,:] + intercepts_tailCrxn[1]
    
    two_hundred_run_233231_avg = ma.mean(five_point_avg[2
                                                    ,:]/five_point_avg[0,:])
    MixPa_233231_avg.append(two_hundred_run_233231_avg)
    two_hundred_run_233231_std = ma.std(five_point_avg[2
                                                   ,:]/five_point_avg[0,:])/ma.sqrt(five_point_avg.shape[1])
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
    print("Identified the following files as no chem mixPa:")
    pprint.pprint(names)
    print('\n')
    
if nochem_mixPa_flag:
    file_name = names[0]
    five_point_avg = return_five_point_avg(file_name)
    
    # tail correction
    five_point_avg[0,:] -= slopes_tailCrxn[0] * five_point_avg[1,:] + intercepts_tailCrxn[0]
    five_point_avg[2,:] -= slopes_tailCrxn[1] * five_point_avg[1,:] + intercepts_tailCrxn[1]
    
    nochem_mixPa_233231_avg = ma.mean(five_point_avg[2
                                                    ,:]/five_point_avg[0,:])
    nochem_mixPa_233231_std = ma.std(five_point_avg[2
                                                   ,:]/five_point_avg[0,:])/ma.sqrt(five_point_avg.shape[1])
    nochem_mixPa_233231_RSD = two_hundred_run_233231_std/two_hundred_run_233231_avg

#%% sample results

# if this is UTh data file
names = [name for name in file_names if '_Pa.txt' in name or '_Pa_re.txt' in name]
if not names:
    raise RuntimeError('No Pa files found!')
names.sort()
print("Identified the following files as sample files:")
pprint.pprint(names)
print('\n')
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
    
    master[i,0] = ma.mean(five_point_avg[0,:]/five_point_avg[2,:])
    two_hundred_run_231233_std = ma.std(five_point_avg[0
                                                   ,:]/five_point_avg[2,:])/ma.sqrt(five_point_avg.shape[1])
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
SRM = ma.mean(SRM_238235_avg)
SRM_RSD = ma.sqrt((ma.sum((ma.array(SRM_238235_avg) * ma.array(SRM_238235_RSD))**2)))/3/SRM
accepted_238235 = 137.55
accepted_238235_RSD = 0.50*0.01
mass_bias_per_amu = (SRM/accepted_238235-1)/3
mass_bias_per_amu_RSD = ma.sqrt((SRM_RSD**2+accepted_238235_RSD**2))

# MixedPa
avg_233Pa_231Pa = ma.mean(MixPa_233231_avg)
avg_233Pa_231Pa_RSD = ma.sqrt(ma.sum((ma.array(MixPa_233231_avg) * ma.array(MixPa_233231_RSD))**2))/3/avg_233Pa_231Pa
Pa231_cncn_pg_g = 38
Pa231_soln_mass_in_spike_g = 0.1924
Pa233_soln_mass_in_spike_g = 26.5163
Pa233_mass_in_spike_pg_g = avg_233Pa_231Pa*Pa231_cncn_pg_g*Pa231_soln_mass_in_spike_g/Pa233_soln_mass_in_spike_g
decay_days = int(input('Enter the number of days from Pa233 spike preparation to Pa clean up column: '))

## write avg_233Pa_231Pa to txt file for Marty
#with open("avg_233_231Pa.txt", "w") as text_file:
#    text_file.write("Purchase Amount: {}".format(avg_233Pa_231Pa))
#%%
##UnSpike
Pa231_Pa233_mass_bias_crtd = (1+(-2*mass_bias_per_amu))*master[:,0]
Pa231_Pa233_mass_bias_crtd_RSD = ma.sqrt(master[:,1]**2+(2*mass_bias_per_amu_RSD)**2)
if sample_info_type == 'txt':
    Pa233_pg = Pa233_mass_in_spike_pg_g * (sample_info['f4']/1000)
elif sample_info_type == 'xlsx':
    Pa233_pg = Pa233_mass_in_spike_pg_g * (sample_info[4]/1000)
Pa231_pg = Pa233_pg * Pa231_Pa233_mass_bias_crtd * 231 / 233
Pa231_pg_RSD = ma.sqrt(Pa231_Pa233_mass_bias_crtd_RSD**2 + avg_233Pa_231Pa_RSD**2)
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
    avg_233Pa_231Pa_df = pd.DataFrame({'avg_233Pa_231Pa':[decay_days,avg_233Pa_231Pa,nochem_mixPa_233231_avg]},index=[0,1,2])
else:
    avg_233Pa_231Pa_df = pd.DataFrame({'avg_233Pa_231Pa':[decay_days,avg_233Pa_231Pa]},index=[0,1])
export_df = pd.concat([sample_name_df,export_data_df,avg_233Pa_231Pa_df],axis=1)

#%% save to excel
output_file_name = filedialog.asksaveasfilename(title='Save the output file as')
if 'xlsx' not in output_file_name:
    output_file_name = output_file_name + '.xlsx'
export_df.to_excel(output_file_name)