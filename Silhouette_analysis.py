# import stuff for silhouette analysis, REQUIRED for next code block
from tslearn.metrics import cdist_dtw, cdist_soft_dtw_normalized
from sklearn.metrics import silhouette_samples as sklearn_silhouette_samples
from tslearn.utils import to_time_series_dataset, to_time_series
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np


def hello():
    return "Hi from Silhouette_analysis.py"


def bootleg_silhouette_samples(X, labels, metric=None, sample_size=None,
                               metric_params=None, n_jobs=None, verbose=0,
                               random_state=None, **kwds):
    '''
    Basically the same as how tslearn made silhouette_score i.e. by computing
    a list of pairwise distances I think.
    '''
    sklearn_metric = None
    if metric_params is None:
        metric_params_ = {}
    else:
        metric_params_ = metric_params.copy()
    for k in kwds.keys():
        metric_params_[k] = kwds[k]
    if "n_jobs" in metric_params_.keys():
        del metric_params_["n_jobs"]
    if metric == "precomputed":
        sklearn_X = X
    elif metric == "dtw" or metric is None:
        sklearn_X = cdist_dtw(X, n_jobs=n_jobs, verbose=verbose,
                              **metric_params_)
    elif metric == "softdtw":
        sklearn_X = cdist_soft_dtw_normalized(X, **metric_params_)
    # elif metric == "euclidean":
    #     X_ = to_time_series_dataset(X)
    #     X_ = X_.reshape((X.shape[0], -1))
    #     sklearn_X = cdist(X_, X_, metric="euclidean")
    else:
        X_ = to_time_series_dataset(X)
        n, sz, d = X_.shape
        sklearn_X = X_.reshape((n, -1))

        def sklearn_metric(x, y):
            return metric(to_time_series(x.reshape((sz, d)),
                                         remove_nans=True),
                          to_time_series(y.reshape((sz, d)),
                                         remove_nans=True))
    metric = "precomputed" if sklearn_metric is None else sklearn_metric
    return sklearn_silhouette_samples(X=sklearn_X,
                                      labels=labels,
                                      metric=metric,
                                      sample_size=sample_size,
                                      random_state=random_state,
                                      **kwds)


def plot_sample_silhouette_spectrum(sample_silhouette_values,
                                    clusternum: 'int',
                                    labels,
                                    avg_silhouette_score,
                                    save_to=False, save_name=False):
    y_lower = 10
    for i in range(clusternum):  # iterate over current num of clusters
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[labels == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / clusternum)
        plt.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )
        # Label the silhouette plots with their cluster numbers at the middle
        plt.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
    plt.title("The silhouette plot for the various clusters.")
    plt.xlabel("The silhouette coefficient values")
    plt.ylabel("Cluster label")
    # The vertical line for average silhouette score of all the values
    plt.axvline(x=avg_silhouette_score, color="red", linestyle="--")
    plt.yticks([])  # Clear the yaxis labels / ticks
    plt.xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    if save_to:
        plt.savefig(fname=save_to + save_name, dpi=300)
