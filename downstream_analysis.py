import os
import math
import zipfile
import pickle
# Essential Libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import datetime
import geopandas as gpd
from shapely.geometry import Point
# Preprocessing
from sklearn.preprocessing import MinMaxScaler
# Algorithms
# from minisom import MiniSom
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from zipfile import ZipFile

# Load in tslearn
from tslearn.barycenters import dtw_barycenter_averaging

# Load in another of my libraries hopefully in same dir
from MSMC_clustering import Msmc_clustering

def plot_cluster_table(cluster_table=False,
                       Msmc_clustering_obj=False,
                       cluster_label_dict=False,
                       cluster_ts_dict=False,
                       basic_cluster_ts_dict=False,
                       clusts = 8,
                       x_aspect = 64,
                       y_aspect = 8,
                       cols = 4,
                       plot_barycenters = True,
                       max_iter=5):
#     from tslearn.barycenters import dtw_barycenter_averaging

    '''
    Given a cluster table and a Msmc_clustering object (loaded with time series
    data which has been transformed into real time via coal. theory), plot out
    time series curves in their respective clusters. Plotting barycenters is
    optional.

    3 methods of plotting depending on which parameter is used:
    - cluster_label_dict: Dict mapping of sample ID to assigned cluster label
    - cluster_ts_dict: Dict mapping of sample ID to corresponding time series; 
        time series is in the format of the dataframe:
        MSMC_clustering.Msmc_clustering.name2series attribute
    - basic_cluster_ts_dict: Dict mapping of sample ID to corresponding
        time series; timeseries is in the basic format of a numpy array 
        {sample ID: timeseries np.array}
    '''
    if cluster_table and Msmc_clustering_obj:
        # iterate over items in sample2label_series and plot them into axs
        method = 1
        labels = cluster_table.Labels.unique()
        indices = cluster_table.index
    elif cluster_label_dict:  # If cluster assignments are known
        if cluster_ts_dict: # If given mapping of sample to ts in my custom format
            method = 2
        elif basic_cluster_ts_dict: # If given in basic format 
            method = 3
        labels = set(list(cluster_label_dict.values()))
        indices = list(cluster_label_dict.keys())
        
    rows = int(np.ceil(clusts/cols)) # always have 4 columns, row = clusts//4
    fs_x = x_aspect
    fs_y = y_aspect * rows
    fig, axs = plt.subplots(rows, cols, figsize=(fs_x, fs_y))
    plt.suptitle(f"Time Series Clusters", fontsize=40)
    label2series = {label:[] for label in labels} # Maps label to timeseries
    for idx in indices:
        if method == 1: # Given MSMC_clustering obj
            label = cluster_table.loc[idx].Labels 
            df = Msmc_clustering_obj.name2series[cluster_table.loc[idx].Sample].drop(["right_time_boundary"], axis=1)
            label2series[label].append(df.to_numpy())
            NE = df["lambda"].to_numpy()
            time = df["left_time_boundary"].to_numpy()
        elif method == 2: # Given pd.df from an attribute of my MSMC_clustering obj
            label = cluster_label_dict[idx]
            df = cluster_ts_dict[idx].drop(["right_time_boundary"], axis=1)
            label2series[label].append(df.to_numpy())
            NE = df["lambda"].to_numpy()
            time = df["left_time_boundary"].to_numpy()
        elif method == 3: # Given basic sample ID to np.array dict mapping
            label = cluster_label_dict[idx]
            if len(basic_cluster_ts_dict[idx].shape) == 1: # If we have 1D ts
                NE = basic_cluster_ts_dict[idx]
                time = list(range(len(basic_cluster_ts_dict[idx])))
            elif len(basic_cluster_ts_dict[idx].shape) == 2: # If we have 2D ts
                # I'm assuming ts will be [(x_val aka time, y_val aka Ne), ...]
                time, NE = zip(*basic_cluster_ts_dict[idx].tolist())
            label2series[label].append(basic_cluster_ts_dict[idx])
            
        if rows <= 1:
            ax = axs[label%cols] # subplot ax only requires col
        else:
            ax = axs[label//cols, label%cols]
        ax.step(time, NE, "+-", c="grey", alpha=0.4)
        ax.set_title(f"Cluster {label}", fontsize=30)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)

    # After a pass over all data, we can plot barycenters for each cluster
    if plot_barycenters:
        for item in label2series.items():
            label = item[0]
            data = item[1] # This is a df of the MSMC_clustering attribute format again
            barycenter = dtw_barycenter_averaging(data, max_iter=max_iter)
            if len(barycenter.shape) == 3:
                x, y = zip(*barycenter)
            elif len(barycenter.shape) == 2:
                y = barycenter
                x = list(range(len(y)))
            if rows <= 1:
                ax = axs[label%cols] # subplot ax only requires col
            else:
                ax = axs[label//cols, label%cols]
            ax.step(x, y, "+-", c="red", alpha=0.8)
    
    
def grab_silhouette_test_clusterlabels(df_avonet3, translate2avonet3=True, verbose=False, verboser=False):
    '''
    Function returns a dict of dfs containing the clustering labels of each
    combination of gamma and clusternumber

    Also adapts clustering label dataframe's Latin name column to whatever
    Latin names are used in AVONET3. For some reason there are name
    inconsistencies which I had to manually fix
    '''
    path = "../"
    clusteringListPath = path + "results/silhouette-tests/lists/"
    gamma_clustnum_df_holder = dict()
    gammas = [np.round(i, 1) for i in np.arange(0.1, 1, 0.2)]
    gammas.insert(0, 0.0)
    clustnums = np.arange(2, 21, 1)
    for gamma in gammas:
        for clustnum in clustnums:
            tsv_path = clusteringListPath + f"gamma_{gamma}/clusternum_{clustnum}/gamma_{gamma}_clusternum_{clustnum}.tsv"
            gamma_clustnum_df_holder[f"gamma_{gamma}_clustnum_{clustnum}"] = pd.read_csv(tsv_path, sep="\t")
            if f"clusternum_{clustnum}" in os.listdir(clusteringListPath + f"gamma_{gamma}/") and translate2avonet3:
                # Begin adaptation of dataframe indices to known AVONET3 dataframe indices
                # e.g. Edolisoma coerulescens in B10K labeling is AKA Coracina coerulescens on Google and Avonet3
                scinames = ["Gallicolumba beccarii", "Pycnonotus atriceps", "Callaeas cinereus", "Larus maculipennis",
                        "Platysteira castanea", "Coracina coerulescens", "Cettia vulcania", "Oceanodroma tethys",
                        "Larus argentatus", "Nectarinia aspasia", "Eupodotis_ruficrista", "Rostratula semicollaris", 
                        "Parus atricapillus", "Megalaima haemacephala", "Phylloscopus sibilatrix", "Dendroica kirtlandii",
                        "Paradoxornis webbianus", "Stachyris dennistouni", "Porzana atra"] # Hand picked/found scientific names (found in alphabetical order of OG names)
                ogNames = sorted([i for i in gamma_clustnum_df_holder[f"gamma_{gamma}_clustnum_{clustnum}"]["Latin name"] if i not in df_avonet3.index]) # Dict with "coloquial" names as keys and scientific names as vals
                # print(ogNames)
                dict_og2sci = {key : scinames[idx].replace(" ", "_") for idx, key in enumerate(ogNames)}
                for idx, name in enumerate(ogNames):
                    name_to_add = scinames[idx].replace(" ", "_")
                    if verboser:
                        print(idx, ":", name, "->", scinames[idx])
                gamma_clustnum_df_holder[f"gamma_{gamma}_clustnum_{clustnum}"]["Latin name"] = [i if i not in dict_og2sci else dict_og2sci[i] for i in gamma_clustnum_df_holder[f"gamma_{gamma}_clustnum_{clustnum}"]["Latin name"]] # Makes changes to df_clusterLabel
                gamma_clustnum_df_holder[f"gamma_{gamma}_clustnum_{clustnum}"] = gamma_clustnum_df_holder[f"gamma_{gamma}_clustnum_{clustnum}"].set_index("Latin name")
            if verbose or verboser:
                print(f"loaded gamma_{gamma}_clusternum_{clustnum}.tsv")
    return gamma_clustnum_df_holder

def grab_and_translate_cluster_labels(df, df_avonet3, translate2avonet3=True, verbose=False, verboser=False):
    '''
    Meant to work with MC.clusterTable type dataframes (df). Slightly more
    general than function above this.
    
    Function returns a dict of dfs containing the clustering labels of each
    combination of gamma and clusternumber.
    Also adapts dataframe's Latin name to whatever Latin names are used in AVONET3
    '''
    df = df.copy()
    if translate2avonet3:
        # Begin adaptation of dataframe indices to known AVONET3 dataframe indices
        # e.g. Edolisoma coerulescens in B10K labeling is AKA Coracina coerulescens on Google and Avonet3
        scinames = ["Gallicolumba beccarii", "Pycnonotus atriceps", "Callaeas cinereus", "Larus maculipennis",
                "Platysteira castanea", "Coracina coerulescens", "Cettia vulcania", "Oceanodroma tethys",
                "Larus argentatus", "Nectarinia aspasia", "Eupodotis_ruficrista", "Rostratula semicollaris", 
                "Parus atricapillus", "Megalaima haemacephala", "Phylloscopus sibilatrix", "Dendroica kirtlandii",
                "Paradoxornis webbianus", "Stachyris dennistouni", "Porzana atra"]  # Hand picked/found scientific names (found in alphabetical order of OG names)
        ogNames = sorted([i for i in df.index if i not in df_avonet3.index]) # Dict with "coloquial" names as keys and scientific names as vals
        dict_og2sci = {key : scinames[idx].replace(" ", "_") for idx, key in enumerate(ogNames)}
        for idx, name in enumerate(ogNames):
            name_to_add = scinames[idx].replace(" ", "_")
            if verboser:
                print(idx, ":", name, "->", scinames[idx])
        df.index = [i if i not in dict_og2sci else dict_og2sci[i] for i in df.index] # Makes changes to df_clusterLabel
    return df

def B10K2AVONET_translation(B10K_latin_names, df_avonet3):
    '''
    B10K_latin_names: list of latin names from the B10K set, " " replaced with "_"
    df_avonet3: the avonet df indexed by latin names with " " replaced with "_"
    scinames: A list of names in avonet3 that might be B10K synonyms. A product 
              of human rummaging
              
    Function returns a dict which assigns the original names in
    B10K_latin_names a potential substitute name which works in AVONET3
    
    '''
    scinames = ["Gallicolumba beccarii", "Pycnonotus atriceps", "Callaeas cinereus", "Larus maculipennis",
                "Platysteira castanea", "Coracina coerulescens", "Cettia vulcania", "Oceanodroma tethys",
                "Larus argentatus", "Nectarinia aspasia", "Eupodotis_ruficrista", "Rostratula semicollaris", 
                "Parus atricapillus", "Megalaima haemacephala", "Phylloscopus sibilatrix", "Dendroica kirtlandii",
                "Paradoxornis webbianus", "Stachyris dennistouni", "Porzana atra"] # Hand picked/found scientific names (found in alphabetical order of OG names)
    B10K2AVONET_translation_dict = dict()
    for lat_name in B10K_latin_names:
        name_parts = lat_name.split("_") # split latin name
        lat_name_surrogates = []
        match_found = False
        for avonet_name in list(df_avonet3.index):
            if lat_name == avonet_name:
                match_found = True
            for name_part in name_parts:
                if name_part in avonet_name:
                    lat_name_surrogates.append(avonet_name)
        if not match_found:
            B10K2AVONET_translation_dict[lat_name] = []
            if len(lat_name_surrogates) == 1: # If only one surrogate is found, might as well use it
                print(f"{lat_name} can use surrogate")
                B10K2AVONET_translation_dict[lat_name] = lat_name_surrogates[0]
            else: # Else, use the EcoNameTranslator library to find it
                print(f"{lat_name} might use a synonym")
                potential_synonyms = synonyms([lat_name])[lat_name][-1]
                print(f"Synonyms: {potential_synonyms}")
    #             print(f"Surrogates: {lat_name_surrogates}")
                for ps in potential_synonyms:
                    if len(lat_name_surrogates)==0:
                        B10K2AVONET_translation_dict[lat_name] = ps.replace(" ", "_")
                    for lns in lat_name_surrogates:
                        if lns == ps.replace(" ", "_"):
                            print("Bingo!")
                            B10K2AVONET_translation_dict[lat_name] = lns
            if len(B10K2AVONET_translation_dict[lat_name]) == 0: # Last ditch effort, check my rummaged list
                print(f"{lat_name} needs a last ditch sub")
                for sciname in scinames:
                    for name_part in name_parts:
                        if name_part in sciname:
                           B10K2AVONET_translation_dict[lat_name] = sciname.replace(" ", "_")
    return B10K2AVONET_translation_dict

from matplotlib import pyplot as plt
from sklearn.decomposition import PCA

def plot_ts_data_set(list_dfs, title, ax):
#     plt.figure(figsize=(12, 6))
    for df in list_dfs:
        time, ne = zip(*df.drop(labels=['right_time_boundary'], axis=1).to_numpy())
        ax.step(time, ne, color="grey")
        ax.set_title(title)
    
def plot_hi_lo_s(avonet3_df, feature, nums):
    '''
    Take in avonet df
    '''
    sorted_list =  list(avonet3_df.sort_values(by=feature).index)
    smallest_few = sorted_list[:num]
    smallest_few = [translated_ts_dict_latin_keyed[i] for i in smallest_few]
    largest_few = sorted_list[-num:]
    largest_few = [translated_ts_dict_latin_keyed[i] for i in largest_few]
    plot_ts_data_set(smallest_few, title=f"Smallest {num}")
    plot_ts_data_set(largest_few, title=f"Largest {num}")
        
import statsmodels.api as sm
from statsmodels.formula.api import ols
pd.set_option("display.max_columns", None)

def anova_on_df(df, avonet_df, traits, formula):
    '''
    Joins your clustering table df (latin name, sample name, labels) with 
    an avonet dataframe which a 1-way ANOVA can be performed on. 
    
    Independent Variable: Everything
    Dependent Variable  : Label
    '''
    df_avonet3.loc[df.index] # Make selection on avonet using df indices
    df_joined = avonet_df.join(df, on=avonet_df.index).dropna(subset=["Labels"]) # Join avonet3 df with chosen df (default df) on "Latin name"
    # Load your data into a pandas dataframe
    df = df_joined.copy()
    df.columns = [i.replace('.','_').replace("-", "_") for i in df.columns]
    # Use statsmodel ols function to fit model using ordinary least squares
    model = ols(formula, df).fit()
    # Perform ANOVA
    aov_table = sm.stats.anova_lm(model, typ=2)
    return (model, aov_table, df)

# Performs translations on all clusters, like one done on new_944
def translate_df_for_avonet(ARI_tested_cluster_table_dict,
                            df_avonet3):
    '''
    Takes in dict keyed by clustering seed. Each val is a cluster table for the
    seed. Translates the cluster table's indices to match the species available
    in avonet3. Outputs a cluster table of the same type as the one taken as
    input, but possibly with different indices (Species Latin names).
    '''
    translated_ARI_cluster_tables = dict()
    for key in ARI_tested_cluster_table_dict.keys():
        translated_ARI_cluster_tables[key] = grab_and_translate_cluster_labels(df=ARI_tested_cluster_table_dict[key],
                                                                               df_avonet3=df_avonet3,
                                                                               translate2avonet3 = "MSMC-Exploratory-Analysis/results/consistency-tests/ARI_testing/",
                                                                               verbose=False,
                                                                               verboser=False)
    return translated_ARI_cluster_tables

def perform_anovas(translated_ARI_cluster_tables,
                   df_avonet3,
                   traits,
                   formula):
    '''
    Using a dict of cluster tables, performs an ANOVA for each cluster table
    (df). Outputs a dict keyed by the seeds used to make each cluster table
    with their corresponding values being results of an ANOVA.
    
    Notes on output:
    Get a dict keyed by random seeds to hold ANOVA results
    Each item in ARI_tested_cluster_ANOVA_dict is a tuple with:
    [0] - model
    [1] - ANOVA results
    [2] - AVONET & cluster labels Joined Table
    '''
    ARI_tested_cluster_ANOVA_dict = dict()
    for key in translated_ARI_cluster_tables.keys():
        df = translated_ARI_cluster_tables[key]
        ARI_tested_cluster_ANOVA_dict[key] = anova_on_df(df=df,
                                                         avonet_df=df_avonet3,
                                                         traits=traits,
                                                         formula=formula)
    return ARI_tested_cluster_ANOVA_dict

def boxplots_for_mass_anova(ARI_tested_cluster_ANOVA_dict):
    '''
    Takes in dicts produced by perform_anovas().
    
    Create boxplots for each feature to see distributions of F and PR(>F)
    Make a df for each feature and compile all F and PR(>F) values into 2 columns.
    Then each feature will have its own 2 boxplots
    '''
    features = list(list(ARI_tested_cluster_ANOVA_dict.values())[0][1].index)[:-1]
    result = []
    for key in ARI_tested_cluster_ANOVA_dict.keys():
        result.append(ARI_tested_cluster_ANOVA_dict[key][1])   

    main_df = pd.concat(result)
    fig, ax = plt.subplots(2, len(features), figsize=(40, 20))

    field="F"
    # ax.set_title(f"Distribution of {field}")
    for idx, feature in enumerate(features):
        df = main_df.loc[feature]
        sns.boxplot(y=df[field], ax=ax[0, idx])
        ax[0,idx].set_xlabel(feature)
        if idx==0:
            ax[0, idx].set_ylabel(field)
        else:
            ax[0, idx].set_ylabel("")
    # ax.set_title(f"Distribution of {field}")
    field="PR(>F)"
    for idx, feature in enumerate(features):
        df = main_df.loc[feature]
        sns.boxplot(y=df[field], ax=ax[1, idx])
        ax[1,idx].set_xlabel(feature)
        if idx==0:
            ax[1, idx].set_ylabel(field)
        else:
            ax[1, idx].set_ylabel("")
    