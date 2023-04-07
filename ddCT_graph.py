import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Take ddCT file and graph means and std")
# Define and use input argments for file, mean, std columns. Set defaults for columns if no input
parser.add_argument("-i", "--input", type=str, help="Name of input file")
parser.add_argument("-y", "--yaxis", type=str, help="Name of column with y values")
parser.add_argument("-s", "--stdev", type=str, help="Name of column with stdev")
parser.add_argument("-n", "--samplename", type=str, help="Name of Sample Name column")
parser.add_argument("-o", "--outfile", type=str, help="Destination file to save")

args = parser.parse_args()

if args.input:
    infile = args.input

if args.yaxis:
    mean = args.yaxis
else:
    mean = "foldchange"

if args.stdev:
    std = args.stdev
else:
    std = "std"

if args.samplename:
    name = args.samplename
else:
    name = "Sample Name"

if args.outfile:
    outfile = args.outfile
else:
    outfile = "ddCT.png"

infile = pd.read_excel(infile)

sns.set(style="ticks")
plot = sns.barplot(data=infile, x=name, y=mean, errorbar="sd", palette="Paired")
plt.savefig(outfile, dpi="300")
