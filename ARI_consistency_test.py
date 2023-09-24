import os
import pickle
import argparse
# Data thingy libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import adjusted_rand_score
# Local clustering class with a bunch of methods
os.chdir("/scratch/nick/MSMC-Curve-Analysis")
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
            x = x.rstrip("\n")
            output.append(x)
    return output


def label2subtable(table, label):
    res = table[table["Labels"]==label]
    return res


def ARI_consistency_test(
                         algo: "str",
                         to_omit: "list<str>",
                         iters: "int",
                         gammas: "list<int>",
                         clusts: "list<int>",
                         oldSeed: "int",
                         newSeeds: "list<int>",
                         directory: "str"="msmc_curve_data/",
                         save_to: "bool/str"=False,
                         plot_everything: "bool"=False,
                         **Msmc_clustering_kwargs,
                         )->"dict":
    '''
    [!] DOES A LOT OF MSMC TS CLUSTERING
    Does pairwise comparisons between a reference clustering (using oldSeed)
    and various other clusterings (using newSeeds). Includes option of taking
    in list of sample names (full file name of sample) as items to omit.
    
    Clusterings can be done over a range of gammas (if using softdtw) and over
    a range of manual clustering sizes
    '''
    seeds = newSeeds
    by_gamma_dict = dict()  # Contains [old_seed_dict, new_seed_dict, randi_dict]
    for gamma in gammas:
        old_seed_dict = {clust: [] for clust in clusts}
        new_seed_dict = {clust: [] for clust in clusts}
        randi_dict = {clust: [] for clust in clusts}
        new_clusters = {clust: dict() for clust in clusts}
        old_clusters = {clust: dict() for clust in clusts}
        for clust in clusts:
            cluster_rt_norm_lenient_og = Msmc_clustering(directory=directory,
                                                         mu=1.4e-9,
                                                         generation_time_path='generation_lengths/',
                                                         to_omit=to_omit,
                                                         exclude_subdirs=["Archive", "mammals_part_1"], 
                                                         manual_cluster_count=clust,
                                                         algo=algo,
                                                         **Msmc_clustering_kwargs) # cluster count by sqrt method is 14
            cluster_rt_norm_lenient_og.cluster_curves(omit_front=0, 
                                                      omit_back=0, 
                                                      cols=4,  
                                                      fs_x=60, 
                                                      fs_y=30,
                                                      metric_params={"gamma" : gamma},
                                                      metric="softdtw",
                                                      random_state=oldSeed,
                                                      plot_everything=False,
                                                      iter=iters,
                                                      )
            for seed in seeds:
                # for seed_idx, seed in enumerate(seeds):
                # print(f"Child clustering {seed_idx}/{len(seeds)}")
                cluster_rt_norm_lenient = Msmc_clustering(directory=directory, 
                                                          mu=1.4e-9, 
                                                          generation_time_path='generation_lengths/', 
                                                          to_omit=to_omit,
                                                          exclude_subdirs=["Archive", "mammals_part_1"], 
                                                          manual_cluster_count=clust,
                                                          algo=algo,
                                                          **Msmc_clustering_kwargs) # cluster count by sqrt method is 14

                cluster_rt_norm_lenient.cluster_curves(omit_front=0, 
                                                       omit_back=0, 
                                                       cols=4,  
                                                       fs_x=60, 
                                                       fs_y=30,
                                                       metric_params={"gamma" : gamma},
                                                       metric="softdtw",
                                                       random_state=seed,
                                                       plot_everything=False,
                                                       iter=iters,
                                                       )

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
                   newSeeds: "list",
                   gamma=0.0,
                   savefig_to=False,
                   savepkl_to=False,
                   plot_everything: "bool" = False) -> "dict":
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
    if plot_everything:
        x, y = zip(*sorted(list(badVoteCount.items()), key = lambda x:x[1]))
        fig = plt.figure(figsize = (40, 40))
        plt.barh(x, y)
        plt.title("Incremental Votes for Samples")
        plt.xlabel("Votes")
        plt.ylabel("Samples")
        if savefig_to:
            plt.savefig(savefig_to)
        plt.close()
    if savepkl_to:
        pickle.dump(badVoteCount, open(savepkl_to, 'wb')) 
    return badVoteCount


def voteOutMovingSample(toAddDict: "dict",
                        badVoteCount: "dict",
                        lower_cutoff: "float" = 0.0,
                        cutoff_last: "bool/int" = False,
                        save_to: "str"=False,
                        accumulate: "bool" = False,
                        both=False)-> "list":
    '''
    toAddDict: like badVoteCount dict but this will be something you want to 
               add to and eventually return at the end of the function. The
               idea is to accumulate votes into this dict.
               
    badVoteCount: {sample_name : votes for being the worst in a clustering}
    
    both: Acts as a flag which will determine if you do both things
    
    Returns either a list of samples to omit or a dict acting as a frequency
    counter for each sample.
    
    Instead of cutting off the last {cutoff_last} number of points (ones with
    the most votes) we accumulate votes in a dict and use some other param to
    cutoff some number of points from the dict (also cutting off the ones with
    the most votes). This might do something different... 
    '''
    # Sort vote dict items by votes (bad samples have many votes)
    if both:
        accumulate = False
    if not accumulate:
        x, y = zip(*sorted(list(badVoteCount.items()), key = lambda x:x[1]))
        y = np.array(y)
        if cutoff_last:
            newOmit = list(x[-cutoff_last:])
        else:
            newOmit = list(x[-len(y[y>lower_cutoff]):]) # Samples with more votes
        newOmit = list(set(newOmit))
        if save_to:
            writeList(save_to, newOmit)
        if not both:
            return newOmit # List of samples to cutoff/omit due to high votes
    if accumulate or both:
        if len(toAddDict.items()) == 0:
            toAddDict = badVoteCount.copy()
        else:
            toAddDict = toAddDict.copy()
            for item in badVoteCount.items():
                sample, votes = item
                toAddDict[sample] += votes
        if not both:
            return toAddDict  # Hist of samples and their badVotes
    return [newOmit, toAddDict]  # Return both thingies


'''
Reliability testing: now a script on mermaid (ariTesting.py)
(MAKE ANOTHER SCIRPT)

REWORK: This uses accumulateVotes function which replaces voteOutMovingSample. 
Doing this makes it so that you don't remove samples with every
oldSeed, newSeed comparison. Instead you accumulate votes over all comparisons
and cutoff some number of the worst votes at the end.
'''


def findVotedOutSamples(sep,
                        suffix,
                        directory,
                        path,
                        algo,
                        iters,
                        gammas,
                        clusts,
                        oldSeeds,
                        seedsList,
                        cutoff_last,
                        cutoff_lastAccumulated,
                        finalOmits,
                        finalOmitsAccumulated,
                        accumulate,
                        bothIncrementalAccumulative,
                        plot_everything,
                        save_accumulation_fig_to,
                        real_time=True,
                        normalize_lambda=True, 
                        log_scale_time=True, 
                        plot_on_log_scale=True, 
                        plot_lightly=True,
                        save_all=True):
    '''
    save_accumulation_fig_to - FIGURE OUT A SUBDIR THAT THIS CAN GO INTO WHEN 
                               SCRIPTIFYING THIS THANG

    Essentially appends samples to finalOmits/finalOmitsAccumulated. If either
    of these are empty lists, output should just be a list of voted-out samples
    since nothing was given to start off.

    Basically a main function that calls on:
    - ARI_consistency_test: Does this to do pairwise clustering comparisons
    - badVoteCounter: Using samples and labels found with ARI_consistency_test,
                      perform badVote counts and return a histogram-like dict
                      containing votes for each sample
    - voteOutMovingSample: This depends on bothIncrementalAccumulative and accumulate
        If not accumulate and not both: (Incremental)
            Returns a list of samples to omit based on incremental cutoff on votes
        If accumulate and not both: (Accumulative)
            toAddDict will be fed into function, have votes added to it for this
            iteration/increment, then return toAddDict for the next iteration.
        If Both:
            Does both Incremental and Accumulative stuff

    After calling these functions, return a pair of lists. First list will 
    correspond to the incremental method of voting out samples. Second list
    will correspond to the accumulative method of voting out samples.
    '''

    if bothIncrementalAccumulative:  # Perform both types of voting
        accumulate = False
    both = bothIncrementalAccumulative  # Pseudonym for bothIncrementalAccumulative
    toAddDict = dict()
    for idx, oldSeed in enumerate(oldSeeds):
        print(f"Batch: {idx}/{len(oldSeeds)}")
        seeds = seedsList[idx]
        ARI_test_dict_path = path + f"ARI_test_dict_oldSeed{oldSeed}.pkl"
        ARI_consistency_test_dict = ARI_consistency_test(sep=sep,
                                                         suffix=suffix,
                                                         directory=directory,
                                                         algo=algo,
                                                         to_omit=omit_test_lenient,
                                                         real_time=real_time,
                                                         normalize_lambda=normalize_lambda, 
                                                         log_scale_time=log_scale_time, 
                                                         plot_on_log_scale=plot_on_log_scale, 
                                                         iters=iters,
                                                         gammas=gammas,
                                                         clusts=clusts,
                                                         oldSeed=oldSeed,
                                                         newSeeds=seeds,
                                                         save_to=ARI_test_dict_path,
                                                         plot_everything=plot_everything)
        badVoteCount_path = path + f"badVote_oldSeed{oldSeed}.pkl"
        badVoteCountFig_path = path + f"badVoteFig_oldSeed{oldSeed}.png"
        clust = clusts[0]
        gamma = gammas[0]
        badVoteCount = badVoteCounter(ARI_consistency_test_dict,
                                      clust,
                                      oldSeed,
                                      seeds,
                                      gamma=gamma,
                                      savefig_to=badVoteCountFig_path,
                                      savepkl_to=badVoteCount_path,
                                      plot_everything=plot_everything)
        if save_all:
            save_to = False
        else:
            save_to = path + f"additionalOmissions/finalARIomits_oldSeed{oldSeed}.txt"
        if not accumulate and not both:
            newOmit = voteOutMovingSample(toAddDict=toAddDict,
                                          badVoteCount=badVoteCount,
                                          cutoff_last=cutoff_last,
                                          save_to=save_to)
            finalOmits += newOmit
        elif accumulate and not both:
            toAddDict = voteOutMovingSample(toAddDict=toAddDict,
                                            badVoteCount=badVoteCount,
                                            accumulate=accumulate)
        elif both:
            pair = voteOutMovingSample(toAddDict=toAddDict,
                                       badVoteCount=badVoteCount,
                                       cutoff_last=cutoff_last,
                                       save_to=save_to,
                                       both=both)
            newOmit, toAddDict = pair
            finalOmits += newOmit
    if accumulate or both:
        x, y = zip(*sorted(list(toAddDict.items()), key=lambda x:x[1]))
        if plot_everything or plot_lightly:
            plt.figure(figsize=(40, 40))
            plt.title("Accumulated Votes for Samples")
            plt.xlabel("Votes")
            plt.ylabel("Samples")
            plt.barh(x, y, color="orange")
            plt.show()
            if save_accumulation_fig_to:
                plt.savefig(save_accumulation_fig_to)
            plt.close()
        finalOmitsAccumulated += x[-cutoff_lastAccumulated:]
        if not both:
            return ([], finalOmitsAccumulated)
        if both:
            return (finalOmits, finalOmitsAccumulated)
    if not accumulate:
        return (finalOmits, [])


def runARItesting(sep,
                  suffix,
                  directory,
                  algo,
                  path,
                  iters,
                  gammas,
                  clusts,
                  oldSeed,
                  seeds,
                  omitListToTest,
                  real_time=True,
                  normalize_lambda=True, 
                  log_scale_time=True, 
                  plot_on_log_scale=True, 
                  accumulative=False):
    '''
    omitListToTest is ideally made from iterative or accumulative voting.

    Runs |seeds| number of clusters for comparison against oldSeed cluster to
    find some distribution of the ARIs given when using an omit list
    (omitListToTest).
    '''
    status = "ACCUMULATIVE" if accumulative else "ITERATIVE"
    print("TESTING", status, "OMIT LIST")
    gamma = gammas[0]
    clust = clusts[0]
    omitListToTest = list(set(omitListToTest))
    if accumulative:
        ARI_test_dict_path = path + f"finalARI_test_dict_oldSeed{oldSeed}_Accumulative.pkl"
        badVoteCount_path = path + f"finalbadVote_oldSeed{oldSeed}_Accumulative.pkl"
        badVoteCountFig_path = path + f"finalbadVoteFig_oldSeed{oldSeed}_Accumulative.png"
        # additionalOmitsPath = path + f"additionalOmissions/omitListToTest_oldSeed{oldSeed}_Accumulative.txt"
        writeListPath = path+"omitListTested_Accumulative.txt"
        final_ARI_df2csvPath = path + "final_ARIs_Accumulative.csv"
        boxplotSavePath = path + "final_ARI_boxplot_Accumulative.png"
    else:
        ARI_test_dict_path = path + f"finalARI_test_dict_oldSeed{oldSeed}.pkl"
        badVoteCount_path = path + f"finalbadVote_oldSeed{oldSeed}.pkl"
        badVoteCountFig_path = path + f"finalbadVoteFig_oldSeed{oldSeed}.png"
        # additionalOmitsPath = path + f"additionalOmissions/omitListToTest_oldSeed{oldSeed}.txt"
        writeListPath = path+"omitListTested.txt"
        final_ARI_df2csvPath = path + "final_ARIs.csv"
        boxplotSavePath = path + "final_ARI_boxplot.png"
    ARI_consistency_test_dict = ARI_consistency_test(sep=sep,
                                                     suffix=suffix,
                                                     directory=directory,
                                                     algo=algo,
                                                     to_omit=omitListToTest,
                                                     real_time=real_time,
                                                     normalize_lambda=normalize_lambda, 
                                                     log_scale_time=log_scale_time, 
                                                     plot_on_log_scale=plot_on_log_scale, 
                                                     iters=iters,
                                                     gammas=gammas,
                                                     clusts=clusts,
                                                     oldSeed=oldSeed,
                                                     newSeeds=seeds,
                                                     save_to=ARI_test_dict_path)
    # Call this, but I don't think anything needs to be saved
    badVoteCounter(ARI_consistency_test_dict,
                   clust,
                   oldSeed,
                   seeds,
                   gamma=gamma,
                   savefig_to=badVoteCountFig_path,
                   savepkl_to=badVoteCount_path)

    writeList(save_to=writeListPath, inlist=omitListToTest)
    ARI_consistency_test_dict[gamma][2]
    randi_dict = ARI_consistency_test_dict[gamma][2]  # list of ARIs sampled
    df = pd.DataFrame.from_dict(randi_dict)
    df = df.rename(columns={clust : "ARIs"})
    df.to_csv(final_ARI_df2csvPath)
    box_plot = sns.boxplot(x=df["ARIs"])
    if accumulative:
        plt.title("ARI Distribution of Accumulative Omit List")
    else:
        plt.title("ARI Distribution of Iterative Omit List")
    figgy = box_plot.get_figure()
    figgy.savefig(boxplotSavePath)
    plt.close()
    return df
# BEGIN TESTING
parser = argparse.ArgumentParser()
parser.add_argument("--spec", default="", type=str, help="Add a custom string to end the directory name")
parser.add_argument("--sep", default="\t", type=str, help="Choose sep or delimeter of data used as input")
parser.add_argument("--suffix", default=".tsv", type=str, help="Choose suffix of data used as input")
parser.add_argument("--directory", default="data/", type=str, help="Choose directory containing clustering data")
parser.add_argument("--algo", default="kmeans", type=str, help="Choose clustering algorithm")
parser.add_argument("--iters", default=50, type=int, help="Iterations to run for clustering")
parser.add_argument("--gamma", default=0.0, type=float,  help="Choose gamma for soft-dtw if using it")
parser.add_argument("--clust", default=8, type=int, help="Number of cluster labels to generate")
parser.add_argument("--cutoff_last", default=5,type=int, help="After each iteration of vote counting from cluster comparison, cutoff this many samples")
parser.add_argument("--cutoff_lastAccumulated", default=70, type=int, help="After accumulating votes over all iterations of cluster comparison, cutoff this many samples")
parser.add_argument("--batches", default=40, type=int, help="Number of oldSeeds to cluster on when voting")
parser.add_argument("--trials", default=20, type=int, help="Number of newSeeds to compare to each oldSeed in terms of clustering when voting")
parser.add_argument("--test_runs", default=30, type=int, help="Number of tests to perform once omit lists are computed")
parser.add_argument("--path", default="MSMC-Exploratory-Analysis/results/consistency-tests/ARI_testing/", type=str, help="Path to save files to")
parser.add_argument("--accumulative", default=False,  action="store_true", help="True if you want to use accumulative voting, False if you want to use iterative voting")
parser.add_argument("--both", default=False,  action="store_true", help="True if you want to use both iterative and accumulative voting")
parser.add_argument("--plot_everything", default=False,  action="store_true", help="True if you want to plot everything")
parser.add_argument("--plot_lightly", default=False,  action="store_true", help="True if you want to plot just final voting histogram")
parser.add_argument("--real_time", default=False,  action="store_true", help="True if you want to transform time field of input data into real time via coalescent scaling (Choose False if data is already in real time)")
parser.add_argument("--normalize_lambda", default=False,  action="store_true", help="True if you want to normalize input data lambda/Ne field onto [0, 1] (Choose False if data is already [0, 1] normalized])")
parser.add_argument("--log_scale_time", default=False,  action="store_true", help="True if you want to logscale time field of input data")
parser.add_argument("--plot_on_log_scale", default=False,  action="store_true", help="True if you want to plot on logscale")


parser.add_argument("--save_all", default=False,  action="store_true", help="True if you want to save all figures and data structures made after clustering")

args = parser.parse_args()
print("args.real_time", args.real_time)
print("accumulate", args.accumulative)
print("bothIncrementalAccumulative", args.both)
gammas = [args.gamma]
clusts = [args.clust]
oldSeeds, seedsList = ARITestSeedGenerator(args.batches, args.trials)  # Used for making voting lists
seeds, _ = ARITestSeedGenerator(args.test_runs, 1)  # Used for evaluating voting lists
print(f"Testing runs to perform: {len(seeds)}")
specifier = f"clusts_{args.clust}_batches_{args.batches}_trials_{args.trials}_iter_cutoff_{args.cutoff_last}_acc_cutoff_{args.cutoff_lastAccumulated}_test_runs_{len(seeds)}_"+args.spec
if specifier not in os.listdir(args.path):  # Specifier makes path specific to settings
    os.mkdir(args.path + specifier)
path = args.path + specifier + "/"

if args.accumulative or args.both:
    accumulative_specifier = "accumulative/"
    if accumulative_specifier not in os.listdir(path): 
        os.mkdir(path + accumulative_specifier)

print("PATH: ", path)
# Params for omission lists
finalOmits = omit_test_lenient.copy()  # finalOmits should add cutoff_last*len(oldSeeds)
finalOmitsAccumulated = omit_test_lenient.copy()
if args.accumulative:
    save_accumulation_fig_to = path + accumulative_specifier + "accumulative_voting_fig.png"
else:
    save_accumulation_fig_to = path + "accumulative_voting_fig.png"
    
finalOmits, finalOmitsAccumulated = findVotedOutSamples(sep=args.sep,
                                                        suffix=args.suffix,
                                                        directory=args.directory,
                                                        path=path,
                                                        algo=args.algo,
                                                        iters=args.iters,
                                                        gammas=gammas,
                                                        clusts=clusts,
                                                        oldSeeds=oldSeeds,
                                                        seedsList=seedsList,
                                                        cutoff_last=args.cutoff_last,
                                                        cutoff_lastAccumulated=args.cutoff_lastAccumulated,
                                                        finalOmits=finalOmits,
                                                        finalOmitsAccumulated=finalOmitsAccumulated,
                                                        save_accumulation_fig_to=save_accumulation_fig_to,
                                                        real_time=args.real_time,
                                                        normalize_lambda=args.normalize_lambda, 
                                                        log_scale_time=args.log_scale_time, 
                                                        plot_on_log_scale=args.plot_on_log_scale, 
                                                        accumulate=args.accumulative,
                                                        bothIncrementalAccumulative=args.both,
                                                        plot_everything=args.plot_everything,
                                                        plot_lightly=args.plot_lightly,
                                                        save_all=args.save_all)

print("FINISHED VOTING OUT SAMPLES")
print(f"len finalOmits: {len(finalOmits)}")
print(f"len finalOmitsAccumulated: {len(finalOmitsAccumulated)}")
# Now we need to test the ARIs of finalOmits as a to_omit list
finalOmits = list(set(finalOmits))
finalOmitsAccumulated = list(set(finalOmitsAccumulated))

oldSeed = 999

# seeds = [100]

# PROBABLY CALL runARItesting here
if args.both:
    it = runARItesting(sep=args.sep,
                       suffix=args.suffix,
                       directory=args.directory,
                       algo=args.algo,
                       path=path,
                       iters=args.iters,
                       gammas=gammas,
                       clusts=clusts,
                       oldSeed=oldSeed,
                       seeds=seeds,
                       omitListToTest=finalOmits,
                       real_time=args.real_time,
                       normalize_lambda=args.normalize_lambda, 
                       log_scale_time=args.log_scale_time, 
                       plot_on_log_scale=args.plot_on_log_scale, 
                       accumulative=False)
    ac = runARItesting(sep=args.sep,
                       suffix=args.suffix,
                       directory=args.directory,
                       algo=args.algo,
                       path=path + accumulative_specifier,
                       iters=args.iters,
                       gammas=gammas,
                       clusts=clusts,
                       oldSeed=oldSeed,
                       seeds=seeds,
                       omitListToTest=finalOmitsAccumulated,
                       real_time=args.real_time,
                       normalize_lambda=args.normalize_lambda, 
                       log_scale_time=args.log_scale_time, 
                       plot_on_log_scale=args.plot_on_log_scale, 
                       accumulative=True)
if not args.both and args.accumulative:
    ac = runARItesting(sep=args.sep,
                       suffix=args.suffix,
                       directory=args.directory,
                       algo=args.algo,
                       path=path + accumulative_specifier,
                       iters=args.iters,
                       gammas=gammas,
                       clusts=clusts,
                       oldSeed=oldSeed,
                       seeds=seeds,
                       omitListToTest=finalOmitsAccumulated,
                       real_time=args.real_time,
                       normalize_lambda=args.normalize_lambda, 
                       log_scale_time=args.log_scale_time, 
                       plot_on_log_scale=args.plot_on_log_scale, 
                       accumulative=True)
if not args.both and not args.accumulative:
    it = runARItesting(sep=args.sep,
                       suffix=args.suffix,
                       directory=args.directory,
                       algo=args.algo,
                       path=path,
                       iters=args.iters,
                       gammas=gammas,
                       clusts=clusts,
                       oldSeed=oldSeed,
                       seeds=seeds,
                       omitListToTest=finalOmits,
                       real_time=args.real_time,
                       normalize_lambda=args.normalize_lambda, 
                       log_scale_time=args.log_scale_time, 
                       plot_on_log_scale=args.plot_on_log_scale, 
                       accumulative=False)
print("FINISHED TESTING ARI TESTING")
print("args.both", args.both)
print("args.accumulative", args.accumulative)
print(f"seeds: {seeds}")