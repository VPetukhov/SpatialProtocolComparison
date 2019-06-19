import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame, Series


def get_expression_value_fracs(expression: np.ndarray, max_val: int = 10, round: int = 3):
    expression = np.array(expression)
    n_transcripts_per_gene = expression.flatten().astype(int)
    n_transcripts_per_gene = n_transcripts_per_gene[n_transcripts_per_gene > 0]
    n_transcripts_per_gene[n_transcripts_per_gene > max_val] = max_val
    n_transcripts_per_gene = n_transcripts_per_gene.astype(str)
    n_transcripts_per_gene[n_transcripts_per_gene == str(max_val)] = ">={}".format(max_val)
    return (Series.value_counts(n_transcripts_per_gene).sort_index() / n_transcripts_per_gene.size * 100).round(
        round).map("{}%".format)


def paired_hist(vals1, vals2, ax, xlabel, ylabel, xlim=None, legend_loc='upper right', do_paired=True):
    if xlim is not None:
        vals1 = vals1[vals1 < xlim]
        vals2 = vals2[vals2 < xlim]

    if do_paired:
        ax.hist(np.array(vals1), bins=50, label=">1")

    ax.hist(np.array(vals2), bins=50, alpha=0.5, label="All")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if legend_loc is not None and do_paired:
        ax.legend(loc=legend_loc)

    if xlim is not None:
        ax.set_xlim(0, xlim)


def plot_expression_metrics(expression: np.ndarray, xlim1=None, xlim2=None, xlim3=None, figsize=(15, 3), do_paired=True):
    expression = np.array(expression)
    fig, axes = plt.subplots(ncols=3, figsize=figsize)

    edf = np.array(expression)
    edf[edf == 1] = 0
    paired_hist(edf.sum(axis=1), expression.sum(axis=1),
                axes[0], xlabel="#Transcripts per cell", ylabel="#Cells", xlim=xlim1, do_paired=do_paired)

    paired_hist((expression > 1).sum(axis=1), (expression > 0).sum(axis=1),
                axes[1], xlabel="#Genes per cell", ylabel="#Cells", xlim=xlim2, do_paired=do_paired)
    paired_hist((expression > 1).sum(axis=0), (expression > 0).sum(axis=0),
                axes[2], xlabel="#Cells per gene", ylabel="#Genes", xlim=xlim3, do_paired=do_paired)

    plt.tight_layout()
    return fig


def get_scalar_metrics(expression: np.ndarray):
    expression = np.array(expression)
    return DataFrame({"Sparsity Level": [(expression == 0).mean()],
                      "#Cells": [expression.shape[0]],
                      "#Genes": [expression.shape[1]],
                      "#Transcripts": [expression.sum()]})
