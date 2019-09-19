import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.preprocessing import Imputer
from difflib import SequenceMatcher
import seaborn as sns
from statistics import mean 
from ast import literal_eval
from scipy import stats


def extract_boar_teloFISH_as_list(path):
    """
    FUNCTION FOR PULLING KELLY'S TELOFISH DATA FOR 40 BOARS into a LIST.. TO BE MADE INTO A DATAFRAME & JOINED W/ MAIN DATAFRAME if possible
    These excel files take forever to load.. the objective here is to synthesize all the excel files for 
    telomere FISH data into one dataframe, then save that dataframe to csv file to be retrieved later
    loading one whole csv file containing all the data will be much, much faster than loading the parts of the whole 
    Along the way, we'll normalize the teloFISH data using controls internal to each excel file
    """
    boar_teloFISH_list = []
    for file in os.scandir(path):
        if 'Hyb' in file.name:
            print(f'Handling {file.name}...')
            full_name = path + file.name
            telo_excel_dict = pd.read_excel(full_name, sheet_name=None, skiprows=4, usecols=[3], nrows=5000)
            if 'Telomere Template' in telo_excel_dict.keys():
                del telo_excel_dict['Telomere Template']

            excel_file_list = []
            for sample_id, telos in telo_excel_dict.items():
                telos_cleaned = clean_individ_telos(telos)

                if sample_id != 'Control': 
                    excel_file_list.append([sample_id, telos_cleaned.values, np.mean(telos_cleaned)]) 
                elif sample_id == 'Control':
                    control_value = np.mean(telos_cleaned)

            
            #normalize teloFISH values by control value 
            for sample in excel_file_list:
                #normalize individual telos
                sample[1] = np.divide(sample[1], control_value)
                
                #normalize telo means
                sample[2] = np.divide(sample[2], control_value)
                boar_teloFISH_list.append(sample)
#                 print(f'{sample[0]} finished..')
            
    print('Finished collecting boar teloFISH data')
    return boar_teloFISH_list

def gen_missing_values_andimpute_or_randomsampledown(n_cells, telosPercell, df):
    
    max_telos = n_cells * telosPercell
    half_telos = (n_cells * telosPercell) / 2

    if df.size > max_telos:
        df_sampled = df.sample(max_telos)
        return df_sampled

    if df.size > 25 and df.size <= half_telos:
        missing_data_difference = abs( (n_cells * telosPercell) - df.size )
        rsampled = df.sample(missing_data_difference, replace=True, random_state=28)
        concat_ed = pd.concat([rsampled, df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        return concat_ed

    if df.size > 25 and df.size < max_telos:
        missing_data_difference = abs( (n_cells * telosPercell) - df.size )
        rsampled = df.sample(missing_data_difference, random_state=28)
        concat_ed = pd.concat([rsampled, df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        return concat_ed
    
    else:
        return df
    
    
def clean_individ_telos(telo_data):
    
    labels=[6, 172, 338, 504, 670, 836, 1002, 1168, 1334, 1500, 1666, 1832, 
    1998, 2164, 2330, 2496, 2662, 2828, 2994, 3160, 3326, 3492, 3658, 3824,
    3990, 4156, 4322, 4488, 4654, 4820]

    labels_offset_by6 = [(x-6) for x in labels]
    
    telo_data = telo_data.drop(labels_offset_by6)
    telo_data = pd.to_numeric(telo_data.iloc[:,0], errors='coerce')
    telo_data = telo_data.dropna(axis=0, how='any')
    telo_data = telo_data.to_frame(name=None)
    telo_data = telo_data[(np.abs(stats.zscore(telo_data)) < 3).all(axis=1)]
    telo_data = pd.Series(telo_data.iloc[:,0])
    telo_data = gen_missing_values_andimpute_or_randomsampledown(30, 160, telo_data)
    telo_data.reset_index(drop=True, inplace=True)
    return telo_data


def remove_dashes_space_sampleIDs(row):

    if '-' in str(row):
        row = str(row).replace('-', '').replace(' ', '')
    
    if ' ' in str(row):
        row = str(row).replace(' ', '')
    
    if 'gps' in str(row):
        row = str(row).replace('gps', '')

    if 'GPS' in str(row):
        row = str(row).replace('GPS', '')
    
    if 'collar' in (row):
        row = str(row).replace('collar', '')
    
    if 'COLLAR' in str(row):
        row = str(row).replace('COLLAR', '')
    
    return row
    
def readable_snake_df_dummy_variables(snake_df):

    Exposure_Status = []
    for row in snake_df['Sample ID']:
        if row.startswith('C'):
            Exposure_Status.append('Control')
        elif row.startswith('E'):
            Exposure_Status.append('Exposed')
    snake_df['Exposure Status'] = Exposure_Status

    ### making dummy variables for snake exposure status
    snake_dum = pd.get_dummies(snake_df['Exposure Status'], prefix='Encoded', drop_first=True)
    snake_df['Encoded Exposed'] = snake_dum
    
    return snake_df


def count_shared_sample_IDs(df1, df2, print_names=None):
    df1_IDs = set(df1['Sample ID'].unique())
    df2_IDs = set(df2['Sample ID'].unique())
    
#     common_IDs = df1_list - (df1_list - df2_list)
    common_IDs = list(df1_IDs & df2_IDs)
    
    print(f'The number of sample IDs in common are: {len(common_IDs)}')
    
    if print_names == 'yes' or print_names == 'Yes':
        print(f'The sample IDs in common are:\n{common_IDs}')
        
        
def average_age_weeks(row):
    
    if '-' in str(row):
        numbers = str(row).split('-')
        average = (int(numbers[1]) + int(numbers[0])) / len(numbers)
        return int(average)
    else:
        return int(row)
        

        
def quartile_cts_rel_to_df1(df1, df2):
    df1 = pd.DataFrame(df1)
    df2 = pd.DataFrame(df2)
    
    # count how many instances in df2 are below the 0.25 quantile of df1
    quartile_1 = df2[df2 <= df1.quantile(0.25)].count()
    
    # count how many instances in df2 are within the 0.25 - 0.75 range quantile of df1
    quartile_2_3 = df2[(df2 > df1.quantile(0.25)) & (df2 < df1.quantile(0.75))].count()

     # count how many instances in df2 are above 0.75 range quantile of df1
    quartile_4 = df2[df2 >= df1.quantile(0.75)].count()
    
    # return counts of values
    return int(quartile_1.values), int(quartile_2_3.values), int(quartile_4.values)



def make_quartiles_columns(total_boar_telos, df):
    
    pos_1, pos_2, pos_3 = 17, 18, 19
    sample_id, telo_data = 0, 1

    for i, row in df.iterrows():
        
        boar_sample_telos = row[telo_data]
        
        df.iat[i, pos_1], df.iat[i, pos_2], df.iat[i, pos_3] = (quartile_cts_rel_to_df1(total_boar_telos, boar_sample_telos))
            
    return df