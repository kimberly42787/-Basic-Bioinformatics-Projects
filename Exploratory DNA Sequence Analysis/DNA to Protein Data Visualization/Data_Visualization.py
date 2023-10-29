# Import the python files 
from main_Biopython import *
from Toolkit_Biopython import *

# Import all of the packages needed for data visualization
import seaborn as sns
import matplotlib.pyplot as plt

# Graph nucleotide composition
def nucl_bar(data_df):
    sns.set(style="white", font="sans serif")
    pal = sns.color_palette(n_colors=3)
    fig_dims = (7, 5.5)

    f, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True, figsize=fig_dims)
    ax1 = sns.barplot(x=data_df[data_df.columns[0]], 
                      y=data_df[data_df.columns[1]],
                      hue=data_df[data_df.columns[2]], 
                      data=merge_df, 
                      palette=pal, 
                      ax=ax1)
    # we basically do the same thing again for the second plot
    ax2 = sns.barplot(x=data_df[data_df.columns[0]], 
                      y=data_df[data_df.columns[1]],
                      hue=data_df[data_df.columns[2]], 
                      data=merge_df, 
                      palette=pal, 
                      ax=ax2)

    # Set Parameters for y- axis
    ax1.set_ylim(4600, 5000)
    ax2.set_ylim(0, 4500)
    # the upper part does not need its own x axis as it shares one with the lower part
    ax1.get_xaxis().set_visible(False)

    # by default, each part will get its own y-axis label, so remove the y label for both subplots
    ax1.set_ylabel("")
    ax2.set_ylabel("")
    # Add in a y-axis label
    f.text(0.04, 0.50, "Composition", va="center", rotation="vertical")

    # by default, seaborn also gives each subplot its own legend, so remove both default legends 
    ax1.get_legend().remove()
    ax2.get_legend().remove()
    # create a new legend and put it to the side of the figure (also requires trial and error)
    ax2.legend(loc=(.68,1.75), title="Sequence")

    # Add some ticks on the top of the upper part and bottom of the lower part for style
    ax1.xaxis.tick_top()
    ax2.xaxis.tick_bottom()

    f.subplots_adjust(left=0.15, right=0.85, bottom=0.15, top=0.90) # Adjust plot to for aesthetics
    for container in ax1.containers:
        ax1.bar_label(container, size = 8)
    plt.show()

nucl_bar(merge_df)

# Make a bar graph of the Kmers
def kmers_bar(data_df):
    sns.set_context("paper")
    fig, ax = plt.subplots(figsize = (14,7.5))
    sns.barplot(
        x="Kmer", 
        y="# of Occurrence", 
        hue = "Sequence",
        ax=ax, 
        data=data_df,
        palette = "Set2")
    plt.xticks(rotation=45,horizontalalignment='right',fontweight='light')
    plt.show()

kmers_bar(merge_kmer)

# Make a bar graph of Amino Acids
def codon_bar(data_df1, data_df2):
    sns.set(style="white", font="sans serif")
    pal = sns.color_palette(n_colors=3)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6.5))

    ax1 = sns.barplot(x=data_df1[data_df1.columns[0]], 
                      y=data_df1[data_df1.columns[1]],
                      data=data_df1, 
                      palette=pal, 
                      ax=ax1)

    ax1.grid(axis='y')
    ax1.set_xlabel('Amino Acids')
    ax1.set_ylabel('# of Occurrence')
    ax1.set_title('Sequence 1')
    for container in ax1.containers:
        ax1.bar_label(container, size = 6.5)
    plt.setp(ax1.get_xticklabels(), rotation=45);

    ax2 = sns.barplot(x=data_df2[data_df2.columns[0]], 
                      y=data_df2[data_df2.columns[1]],
                      data=data_df2, 
                      palette=pal, 
                      ax=ax2)

    ax2.grid(axis='y')
    ax2.set_xlabel('Amino Acids')
    ax2.set_ylabel('# of Occurrence')
    ax2.set_title('Sequence 2')
    for container in ax2.containers:
        ax2.bar_label(container, size = 6.5)
    plt.tight_layout()
    plt.setp(ax2.get_xticklabels(), rotation=45)
    plt.show()

codon_bar(protein_count1_df,protein_count2_df)
