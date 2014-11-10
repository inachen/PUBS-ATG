import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
import csv
import numpy as np
import math
import matplotlib.cm as cm
from mpltools import style
from mpltools import color
import cPickle as pic
import pandas as pd

style.use('ggplot')

def plot_gc():
    # dat_gc = np.genfromtxt("data_files/growthcurve.csv", delimiter=',', names=True)
    df_gc=pd.read_csv("data_files/growthcurve.csv", sep=',',header=0)

    fig = plt.figure(0)
    plt.clf()
    plt.plot(df_gc["Time Point"], df_gc["1X"], label="1X")
    plt.plot(df_gc["Time Point"], df_gc["0.25X"], label="0.25X")
    plt.plot(df_gc["Time Point"], df_gc["Control"], label="Control")
    plt.legend()
    plt.xlabel("Time (h)")
    plt.ylabel("Growth Rate (OD-600)")
    plt.title("Growth of Yeast with Perturbation")
    fig.savefig('gc_plot.png')

def plot_exp_corr():
    df_exp=pd.read_csv("data_files/exp_corr.csv", sep=',',header=0)

    print df_exp

    m, b = np.polyfit(df_exp["Day 1"], df_exp["Day 2"], 1)

    fig = plt.figure(1)
    plt.clf()
    plt.plot(df_exp["Day 1"], df_exp["Day 2"], '.')
    plt.plot(df_exp["Day 1"], m*df_exp["Day 1"] + b)
    plt.xlabel("Day 1 Fitness Values")
    plt.ylabel("Day 2 Fitness Values")
    plt.title("Correlation Between Day 1 and Day 2 Fitness Values")
    fig.savefig('plot_exp_corr.png')

def demult_viz():
    index = ['TTGACT', 'GGAACT', 'TGACAT', 'GGACGG', 'CTCTAC', 'GCGGAC']
    counts = [2730894, 2785106, 2348801, 3294492, 3413462, 4608969]
    # index = ['GGAACT', 'CTCTAC', 'TGACAT', 'GCGGAC', 'GGACGG', 'TTGACT']
    # counts = [2785106, 3413462, 2348801, 4608969, 3294492, 2730894]

    x = range(6)

    fig = plt.figure(0)
    plt.clf()
    plt.bar(x, counts)
    plt.ylabel("Counts")
    plt.title("Demultiplexing Counts")
    plt.xticks([i+0.5 for i in x], index)
    fig.savefig('index_counts.png')

def qual_viz():
    index = ['Total Count', '- Bad Qual', 'Identical', 'Near', 'Ambiguous']
    counts = [19181724, 18240351, 16633987, 343010, 109343]
    # index = ['GGAACT', 'CTCTAC', 'TGACAT', 'GCGGAC', 'GGACGG', 'TTGACT']
    # counts = [2785106, 3413462, 2348801, 4608969, 3294492, 2730894]

    x = range(5)

    fig = plt.figure(0)
    plt.clf()
    plt.bar(x, counts, color ='#E24A33')
    plt.ylabel("Counts")
    plt.title("Barcode Counts")
    plt.xticks([i+0.5 for i in x], index)
    fig.savefig('qual_counts.png')

    # print plt.rcParams['axes.color_cycle']

demult_viz()