import os
import pickle
from sklearn.metrics import adjusted_rand_score
# Local clustering class with a bunch of methods
os.chdir("/scratch/nick/MSMC-Curve-Analysis")
# Data thingy libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import silhouette_samples
from sklearn.metrics import adjusted_rand_score
from MSMC_clustering import Msmc_clustering

# Load in precomputed list of files to filter out (Saves time)
# Consider making the sharpness filter on percent too
regular_path = "lenient-curve-filter"
filter_path = "MSMC-Exploratory-Analysis/results/figures/filtered-out-curves"
kept_path = "MSMC-Exploratory-Analysis/results/figures/kept-curves"

omit_list_path = "MSMC-Exploratory-Analysis/results/lists"
omit_test_lenient_file = "omit_test_lenient.txt"
kept_test_lenient_file = "kept_test_lenient.txt"

if omit_test_lenient_file not in os.listdir(omit_list_path):
    omit_test_lenient = filter_flattness_sharpness(cluster_rt,
                        identical_val_threshold = 0.7, err = 0.01,  # Args for filtering flattness
                        magnitude_jump_threshold = 0.95, points_per_jump = 3,  # Args for filtering magnitude jumps
                        plot = True, min_unique_vals=999, ignore_low_mag=0, 
                        fs_x=10, fs_y=9, save_to_filter=filter_path, save_to_kept=kept_path)  # These are some pretty good settings
    with open(omit_list_path+"/"+omit_test_lenient_file, "w") as myFile:
        for fileName in omit_test_lenient:
            myFile.write(f"{fileName}\n")

else:  # If list exists
    omit_test_lenient = []  # List of names to files which were filtered out due to being elboq shaped
    with open(omit_list_path+"/"+omit_test_lenient_file, "r") as myFile:
        for line in myFile:
            omit_test_lenient.append(line.rstrip())

kept_test_lenient = []  # List of names to files which were kept
for jpg_name in os.listdir(kept_path):
    jpg_name = jpg_name[:-4]  # assuming that all files in "kept_path" are .jpgs and end with the specifier ".jpg"
    kept_test_lenient.append(jpg_name)

# List of latin names is useful for indexing on metadata dfs later on
omit_test_lenient_latin = [name[:name.index("_GC")] for name in omit_test_lenient] 
kept_test_lenient_latin = [name[:name.index("_GC")] for name in kept_test_lenient]

'''
Helper functions
'''


def ARITestSeedGenerator(oldSeedNum, newSeedNum):
    '''
    oldSeedMum specifies number of seeds (batches) to use as reference in
    comparisons between a reference and newSeedNum amount of other clusterings.
    '''
    oldSeeds = list(np.random.randint(low=1, high=9999, size=oldSeedNum))
    seedsList = [list(np.random.randint(low=1, high=9999, size=newSeedNum))
                 for i in range(len(oldSeeds))]
    return oldSeeds, seedsList


def writeList(save_to, inlist):
    with open(save_to, "w") as f:
        for x in inlist:
            f.write(x+"\n")
    return


def readList(path):
    with open(path, "r") as f:
        output = []
        for x in f:
            output.append(x)
    return output


def label2subtable(table, label):
    res = table[table["Labels"] == label]
    return res


def ARI_consistency_test(algo: "str",
                         to_omit: "list<str>",
                         iters: "int",
                         gammas: "list<int>",
                         clusts: "list<int>",
                         oldSeed: "int",
                         newSeeds: "list<int>",
                         save_to: "bool/str" = False,
                         cutoff_last=70) -> "dict":
    '''
    Does pairwise comparisons between a reference clustering (using oldSeed)
    and various other clusterings (using newSeeds). Includes option of taking
    in list of sample names (full file name of sample) as items to omit.
    
    Clusterings can be done over a range of gammas (if using softdtw) and over
    a range of manual clustering sizes
    '''

    seeds = newSeeds
    by_gamma_dict = dict() # Contains [old_seed_dict, new_seed_dict, randi_dict]
    for gamma in gammas:
        old_seed_dict = {clust : [] for clust in clusts}
        new_seed_dict = {clust : [] for clust in clusts}
        randi_dict =    {clust : [] for clust in clusts}
        new_clusters = {clust : dict() for clust in clusts}
        old_clusters = {clust : dict() for clust in clusts}
        for clust in clusts:
            cluster_rt_norm_lenient_og0 = Msmc_clustering(directory="msmc_curve_data/",
                                                    mu=1.4e-9,
                                                    generation_time_path='generation_lengths/',
                                                    real_time=True,
                                                    to_omit=to_omit,
                                                    normalize_lambda=True,
                                                    log_scale_time=True,
                                                    plot_on_log_scale=True,
                                                    exclude_subdirs=["Archive", "mammals_part_1"], 
                                                    manual_cluster_count=clust,
                                                    algo=algo) # cluster count by sqrt method is 14
            # NEW
            new_to_omit = cluster_rt_norm_lenient_og0.namesofMySeries.copy()
            random_selection = list(np.random.choice(new_to_omit, size=cutoff_last, replace=False))
            # print(to_omit)
            # print()
            # print(random_selection)
            finalOmits = to_omit + random_selection
            
            cluster_rt_norm_lenient_og = Msmc_clustering(directory="msmc_curve_data/",
                                                    mu=1.4e-9,
                                                    generation_time_path='generation_lengths/',
                                                    real_time=True,
                                                    to_omit=finalOmits,
                                                    normalize_lambda=True,
                                                    log_scale_time=True,
                                                    plot_on_log_scale=True,
                                                    exclude_subdirs=["Archive", "mammals_part_1"], 
                                                    manual_cluster_count=clust,
                                                    algo=algo) # cluster count by sqrt method is 14
            cluster_rt_norm_lenient_og.cluster_curves(omit_front=0,
                                                      omit_back=0,
                                                      cols=4, 
                                                      fs_x=60,
                                                      fs_y=30,
                                                      metric_params={"gamma" : gamma},
                                                      metric="softdtw",
                                                      random_state=oldSeed,
                                                      plot_everything=False,
                                                      iter=iters)
            for seed in seeds:
                new_to_omit_seed = cluster_rt_norm_lenient_og0.namesofMySeries.copy()
                random_selection = list(np.random.choice(new_to_omit_seed, size=cutoff_last, replace=False))
                finalOmits = to_omit + random_selection
                cluster_rt_norm_lenient = Msmc_clustering(directory="msmc_curve_data/",
                                                        mu=1.4e-9, 
                                                        generation_time_path='generation_lengths/', 
                                                        real_time=True,
                                                        to_omit=finalOmits,
                                                        normalize_lambda=True, 
                                                        log_scale_time=True, 
                                                        plot_on_log_scale=True, 
                                                        exclude_subdirs=["Archive", "mammals_part_1"], 
                                                        manual_cluster_count=clust,
                                                        algo=algo) # cluster count by sqrt method is 14

                cluster_rt_norm_lenient.cluster_curves(omit_front=0, 
                                                       omit_back=0, 
                                                       cols=4,  
                                                       fs_x=60, 
                                                       fs_y=30,
                                                       metric_params={"gamma" : gamma},
                                                       metric="softdtw",
                                                       random_state=seed,
                                                       plot_everything=False,
                                                       iter=iters)

                old = cluster_rt_norm_lenient_og.dtw_labels
                old_seed_dict[clust].append(old)

                new = cluster_rt_norm_lenient.dtw_labels
                new_seed_dict[clust].append(new)
                randi_dict[clust].append(adjusted_rand_score(old, new))
                new_clusters[clust][seed] = cluster_rt_norm_lenient
                print(f"[RAND INDEX]: Random seed comparison of {clust} clusters: {adjusted_rand_score(old, new)}")
            old_clusters[clust][oldSeed] = cluster_rt_norm_lenient_og
        by_gamma_dict[gamma] = [old_seed_dict, new_seed_dict, randi_dict, old_clusters, new_clusters]
        if save_to:
            pickle.dump(by_gamma_dict, open(save_to, 'wb'))
    return by_gamma_dict

def badVoteCounter(ARI_consistency_test_dict: "dict",
                   clust: "int",
                   oldSeed: "int",
                   newSeeds: "list<int>",
                   gamma=0.0,
                   savefig_to=False,
                   savepkl_to=False)->"dict":
    '''
    I have seriesDict which is {name: series}

    We have a constant clustering and compare it to N number of random clusterings.

    I will need to separate samples by cluster.
    I add votes to samples for each difference that their constant clustering has
    with each new clustering.

    For each sample in old cluster
        Find clustermates of old sample
        For same sample in new cluster
            Find clustermates of old sample in new cluster
            Count the number of differences between old and new cluster
            Map count of differences to sample name
    '''
    seeds = newSeeds
    old_clusters = ARI_consistency_test_dict[gamma][-2]
    new_clusters = ARI_consistency_test_dict[gamma][-1]
    old = old_clusters[clust][oldSeed]
    oldSampleNames = old.namesofMySeries
    badVoteCount = {sample: 0 for sample in oldSampleNames}
    for oldSampleName in oldSampleNames:
        oldSampleKey = oldSampleName[:oldSampleName.index("_GC")] # Latin name
        oldSampleLabel = old.clusterTable.loc[oldSampleKey]["Labels"]
        oldMates = label2subtable(old.clusterTable, oldSampleLabel) # Clustermates of old sample

        for seed in seeds:
            new = new_clusters[clust][seed]
            newSampleNames = new.namesofMySeries
            newSampleLabel = new.clusterTable.loc[oldSampleKey]["Labels"] # w cluster assignment of old sample
            newMates = label2subtable(new.clusterTable, newSampleLabel)

            difference = list(set(newMates.index) - set(oldMates.index))
            badVoteCount[oldSampleName] += len(difference)

    x, y = zip(*sorted(badVoteCount.items(), key = lambda x:x[1]))
    fig = plt.figure(figsize = (40, 40))
    plt.barh(x, y)
    if savefig_to and savepkl_to:
        plt.savefig(savefig_to)
        plt.close()
        pickle.dump(badVoteCount, open(savepkl_to, 'wb'))
    return badVoteCount


def voteOutMovingSample(badVoteCount: "dict",
                        lower_cutoff: "float" = 0.0,
                        cutoff_last: "bool/int" = False,
                        save_to: "str" = False) -> "list":
    # Sort vote dict items by votes (bad samples have many votes)
    x, y = zip(*sorted(badVoteCount.items(), key=lambda x: x[1]))
    y = np.array(y)
    if cutoff_last:
        newOmit = list(x[-cutoff_last:])
    else:
        newOmit = list(x[-len(y[y>lower_cutoff]):])  # Samples with more votes
    newOmit = list(set(newOmit))
    if save_to:
        writeList(save_to, newOmit)
    return newOmit


# BEGIN TESTING
algo = "kmeans"
iters = 50
gamma = 0.0
gammas = [gamma]
clust = 8
clusts = [clust]
# due to nature of sampling, len(oldSeeds) == len(seedsList) 
cutoff_last = 70  # Number of samples to cutoff for each badVoteCount iteration

path = "MSMC-Exploratory-Analysis/results/consistency-tests/ARI_testing/"
specifier = f"random_dropping_cutoff_last_{cutoff_last}"
if specifier not in os.listdir(path):
    os.mkdir(path + specifier)
path = path + specifier + "/"
print("PATH: ", path)

# cutoff_last = 10  # Number of samples to cutoff for each badVoteCount iteration
# oldSeeds = [502, 111, 234, 454, 222, 323,
#             456, 766, 776, 453, 256, 771]
# seedsList = [[810, 564, 588, 155,  72],
#              [134, 555, 665, 345, 123],
#              [808, 667, 788, 551,  23],
#              [335, 213,  34, 444, 443],
#              [202, 101, 303, 505, 506],
#              [321, 322, 909, 797, 777],
#              [811, 561, 581, 151,  71],
#              [114, 515, 615, 315, 113],
#              [838, 627, 718, 511, 423],
#              [345, 215, 554, 455, 446],
#              [232, 131, 333, 535, 536],
#              [341, 342, 949, 747, 737]]
# finalOmits = omit_test_lenient.copy()  # finalOmits should add cutoff_last*len(oldSeeds)
# for idx, oldSeed in enumerate(oldSeeds):
#     seeds = seedsList[idx]
#     ARI_test_dict_path = path + f"ARI_test_dict_oldSeed{oldSeed}.pkl"
#     ARI_consistency_test_dict = ARI_consistency_test(algo="kmeans",
#                                                      to_omit=omit_test_lenient,
#                                                      iters=50,
#                                                      gammas=gammas,
#                                                      clusts=clusts,
#                                                      oldSeed=oldSeed,
#                                                      newSeeds=seeds,
#                                                      save_to=ARI_test_dict_path)
#     badVoteCount_path = path + f"badVote_oldSeed{oldSeed}.pkl"
#     badVoteCountFig_path = path + f"badVoteFig_oldSeed{oldSeed}.png"
#     badVoteCount = badVoteCounter(ARI_consistency_test_dict,
#                                   clust,
#                                   oldSeed,
#                                   seeds,
#                                   gamma=gamma,
#                                   savefig_to=badVoteCountFig_path,
#                                   savepkl_to=badVoteCount_path)
#     if "additionalOmissions" not in os.listdir(path):
#         os.mkdir(path+"additionalOmissions")
#     additionalOmitsPath = path + f"additionalOmissions/ARIomits_oldSeed{oldSeed}.txt"
#     newOmit = voteOutMovingSample(badVoteCount,
#                                   cutoff_last=cutoff_last,
#                                   save_to=additionalOmitsPath)
#     finalOmits += newOmit

# Now we need to test the ARIs of finalOmits as a to_omit list
finalOmits = omit_test_lenient.copy()
finalOmits = list(set(finalOmits))

oldSeed = 999
testSeeds, _ = ARITestSeedGenerator(1, 100)
seeds = _[0]
# print(seeds)

ARI_test_dict_path = path + f"finalARI_test_dict_oldSeed{oldSeed}.pkl"
ARI_consistency_test_dict = ARI_consistency_test(algo="kmeans",
                                                 to_omit=finalOmits,
                                                 iters=50,
                                                 gammas=gammas,
                                                 clusts=clusts,
                                                 oldSeed=oldSeed,
                                                 newSeeds=seeds,
                                                 save_to=ARI_test_dict_path)
# badVoteCount_path = path + f"finalbadVote_oldSeed{oldSeed}.pkl"
# badVoteCountFig_path = path + f"finalbadVoteFig_oldSeed{oldSeed}.png"
# badVoteCount = badVoteCounter(ARI_consistency_test_dict,
#                               clust,
#                               oldSeed,
#                               seeds,
#                               gamma=gamma,
#                               savefig_to=badVoteCountFig_path,
#                               savepkl_to=badVoteCount_path)
additionalOmitsPath = path + f"additionalOmissions/finalARIomits_oldSeed{oldSeed}.txt"
# newOmit = voteOutMovingSample(badVoteCount,
#                               cutoff_last=cutoff_last,
#                               save_to=additionalOmitsPath)
writeList(save_to=path+"finalOmitsList.txt", inlist=finalOmits)
ARI_consistency_test_dict[gamma][2]
randi_dict = ARI_consistency_test_dict[gamma][2]  # list of ARIs sampled
df = pd.DataFrame.from_dict(randi_dict)
df = df.rename(columns={clust : "ARIs"})
df.to_csv(path + "final_ARIs.csv")
box_plot = sns.boxplot(x=df["ARIs"])
fig = box_plot.get_figure()
fig.savefig(path + "final_ARI_boxplot.png")
