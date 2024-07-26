import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


pd.options.mode.chained_assignment = None

sns.set_theme(style='whitegrid', rc={
    'figure.figsize': (6.4, 6.4),  # figure size in inches
    'axes.labelsize': 14,           # font size of the x and y labels
    'xtick.labelsize': 12,          # font size of the tick labels
    'ytick.labelsize': 12,          # font size of the tick labels
    'legend.fontsize': 12,          # font size of the legend
})

  # set the resolution to 300 DPI

def parse_line(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
    AS_score = None
    ms_score = None
    tp_score = None
    for field in fields[12:]:
        if field.startswith('ms'):
            ms_score = int(field.split(':')[2])
        if field.startswith('AS'):
            AS_score = int(field.split(':')[2])
        if field.startswith('tp'):
            tp_score = field.split(':')[2]
            break
    return read_id, AS_score, ms_score, tp_score


def paf_to_dataframe(file):
    print(file)
    # with open(file, 'r') as f:
    #     data = [parse_line(line) for line in list(f)[:10000]]
    data = [parse_line(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id', 'AS_score', 'ms_score', 'tp_score'])
    # return a test df
    #
    #df = df.head(1000)
    return df

def main():
    parser = argparse.ArgumentParser(description='Plot AS and MS scores from minimap2 PAF file')
    parser.add_argument('paf_file', type=str, help='Path to the PAF file')
    parser.add_argument('output_prefix', type=str, help='Prefix for the output file')
    args = parser.parse_args()

    df = paf_to_dataframe(args.paf_file)

    # Plot AS and MS scores
    with PdfPages(f'{args.output_prefix}_AS_ms_scores.pdf') as pdf:
        fig, ax = plt.subplots()
        sns.scatterplot(data=df, x='AS_score', y='ms_score', ax=ax, size = 0.01, hue = 'tp_score')
        ax.set_xlabel('AS score')
        ax.set_ylabel('ms score')
        # dont show ax labels
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        # add a legend
        ax.legend(title='tp score', loc='upper right')

        pdf.savefig(fig)
        plt.close()

if __name__ == '__main__':
    main()



# python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/plot_ms_AS_scores.py  out/aligned/align_competetive_AG06213_PAB_N200_allhaps.paf  Orangutan