
import matplotlib.pyplot as plt
import pandas as pd

def plot_celltypes(celltype_scores,out):
    names=list(celltype_scores.keys())
    vals=list(celltype_scores.values())
    plt.bar(names,vals)
    plt.ylabel("fraction")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(out)
