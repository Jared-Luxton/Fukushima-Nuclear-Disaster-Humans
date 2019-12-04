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
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from pygam import LinearGAM, s, l, f


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
        
    if '_' in str(row):
        row = str(row).replace('_', '')
    
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
    
    return str(row)
    
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


def linear_regression_graphs_between_variables(x=None, y=None, data=None, 
                                               hue=None, col=None,
                                               hue_order=None, col_order=None,
                                               snake=False):
    
    if 'Binary' in y:
        ax=sns.lmplot(x=x, y=y, hue=hue, col=col, data=data, logistic=True, 
        height=5.5, aspect=1, scatter_kws={"s": 175, "edgecolor":'black'})
    else:
        ax=sns.lmplot(x=x, y=y, hue=hue, col=col, data=data,
        height=5.5, aspect=1, scatter_kws={"s": 175, "edgecolor":'black'})

    fig = ax.fig 

    ax.set_xlabels(x, fontsize=18)
    ax.set_xticklabels(fontsize=14)
    ax.set_ylabels(y, fontsize=18)
    ax.set_yticklabels(fontsize=14)
    ax.set_titles(size=14)
    
#     if 'Cortisol' in y:
#         ax.set(ylim=(0, 40))

    plt.subplots_adjust(top=0.88)

    if hue == None and col == None:
        fig.suptitle(f'{x} vs.\n {y} in Fukushima Wild Boar', fontsize=18, 
                    )
#         ax.savefig(f"../graphs/{x} vs {y}.png", dpi=400)

    if snake:
        fig.suptitle(f'{x} vs.\n {y} in Fukushima Wild Snake', fontsize=18, 
                    )
            
#     elif hue == 'Sex' and col == 'Sex':
#         fig.suptitle(f'{x} vs. {y}\nper Sex in Fukushima Wild Boar', fontsize=16, weight='bold')
#         fig.legend(fontsize='large')
#         ax.savefig(f"../graphs/{x} vs {y} per sex.png", dpi=400)


def graph_dose_age_vs_telos(df=None, x=None, x2=None, y=None):
    f, axes = plt.subplots(1, 2, figsize=(12,5), sharey=False, sharex=False)
    # dose vs. telomeres
    sns.regplot(x=x, y=y, data=df, ax=axes[0], 
                scatter_kws={'alpha':0.8, 'linewidth':1, 'edgecolor':'black', 's':df['Age (months)']*12, })
    axes[0].set_xlabel(x, fontsize=14)
    axes[0].set_ylabel(y, fontsize=14)
    axes[0].tick_params(labelsize=12)


    # age vs. telomeres
    sns.regplot(x=x2, y=y, data=df, ax=axes[1], 
                scatter_kws={'alpha':0.8, 'linewidth':1, 'edgecolor':'black', 's':175, })
    axes[1].set_xlabel(x2, fontsize=14)
    axes[1].set_xlim(-4,55)
    axes[1].set_ylabel(y, fontsize=14)
    if y == 'teloFISH means':
        axes[1].set_ylim(0.2,1.6)
    if y == 'Mean Telomere Length (qPCR)':
        axes[1].set_ylim(0.6,1.8)
    axes[1].tick_params(labelsize=12)
            
            
def score_linear_regressions(x=None, y=None, data=None):
    
#     sexes = ['Male', 'Female', 'Overall']
    sexes = ['Overall']

    for sex in sexes:
        if sex == 'Overall':
            X_r = data[x].values.reshape(-1, len(x))
            y_r = data[y].values.reshape(-1, 1)
            regression = LinearRegression().fit(X_r,y_r)
            print(f'Linear regression for {x} vs. {y}:\nOverall R2 is {regression.score(X_r, y_r):.4f}\n')

        else:
            X_r = data[data['Sex'] == sex][x].values.reshape(-1, len(x))
            y_r = data[data['Sex'] == sex][y].values.reshape(-1, 1)
            regression = LinearRegression().fit(X_r,y_r)
            print(f"Linear regression for {x} vs. {y}:\nR2 for {sex}s is {regression.score(X_r, y_r):.4f}")


            
def eval_number(x):
    if x > 15:
        x = 1
        return x
    elif x < 15:
        x = 0
        return x
        

def score_logistic_regressions(x=None, y=None, data=None):
    
#     for y in y_cols:
    
    sexes = [
#         'Male', 
#         'Female', 
        'Overall']

    for sex in sexes:
        if sex == 'Overall':
            X_r = data[x].values.reshape(-1, 1)
            y_r = data[y].values.reshape(-1, )
            log_reg = LogisticRegression(solver='lbfgs')
            regression = log_reg.fit(X_r,y_r)
            print(f'Logistic regression for {x} vs. {y}:\nOverall R2 is {regression.score(X_r, y_r):.4f}\n')

        else:
            X_r = data[data['Sex'] == sex][x].values.reshape(-1, 1)
            y_r = data[data['Sex'] == sex][y].values.reshape(-1,  )
            regression = LinearRegression().fit(X_r,y_r)
            print(f"Logistic regression for {x} vs. {y}:\nR2 for {sex}s is {regression.score(X_r, y_r):.4f}")
            
            
def encode_sex(row):
    if row == 'Male':
        return 1
    elif row == 'Female':
        return 0
    else:
        print(f'ERROR.. row == {row}')
        
        
def merge_return_df_cols_interest(dose_df, cortisol_df, cols_of_interest):
    merge_dose_cortisol = dose_df.merge(cortisol_df, on=['Sample ID'])
    trim_dose_cortisol = merge_dose_cortisol[cols_of_interest].copy()
    return trim_dose_cortisol


def enforce_col_types(df):
    for col in df.columns:
        if col == 'Sample ID' or col == 'Sex':
            df[col] = df[col].astype('str')
        elif col == 'Age (months)':
            df[col] = df[col].astype('int64')
        else:
            df[col] = df[col].astype('float64')
            
            
def male_or_female(row):
    if row == 'M' or row == 'm':
        return 'Male'
    elif row == 'F' or row == 'f':
        return 'Female'
    else:
        print(f'error... row == {row}')
        return np.NaN
            
            
def linear_regression_scores_X_y(df, y, y_name, dose_types):
    """
    specifically for EDA
    """
    for Xn in dose_types:
        features_list = [[Xn], [Xn, 'Age (months)'], [Xn, 'Age (months)', 'encoded sex']]
        for features in features_list:
            X = df[features].values.reshape(-1, len(features))
            fit_lm = LinearRegression().fit(X, y)
            print(f'OLS | {features} vs. {y_name} --> R2: {fit_lm.score(X, y):.4f}')
        print('')
        
            
def fit_gam_plot_dependencies(df=None, features=None, target=None, 
                              basis_1=s, basis_2=False, summary=False):
    X = df[features]
    y = df[target]
    
    if basis_1 and basis_2:
        gam = LinearGAM(basis_1(0, lam=60) + basis_2(1, lam=60), fit_intercept=True).fit(X, y)
    
    elif basis_1:
        gam = LinearGAM(basis_1(0, lam=60), fit_intercept=True).fit(X, y)
        
    else:
        print('no basis called for features.. error')
    
    if summary:
        print(gam.summary())
    plot_gam_partial_dependencies(gam, features, target)
    
    
def plot_gam_partial_dependencies(gam, features, target):
    for i, term in enumerate(gam.terms):
        if term.isintercept:
            continue
            
        XX = gam.generate_X_grid(term=i)
        pdep, confi = gam.partial_dependence(term=i, X=XX, width=0.95)

        plt.figure()
        plt.plot(XX[:, term.feature], pdep)
        plt.plot(XX[:, term.feature], confi, c='r', ls='--')
        plt.xlabel(f'{features[i]}', fontsize=14)
        plt.ylabel(f'{target}', fontsize=14)
        plt.title(f'Functional dependence of Y on X', fontsize=14)
        plt.show()