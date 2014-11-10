import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
from collections import namedtuple
import csv
import numpy as np
import math
import matplotlib.cm as cm
import copy
import os.path
# from mpltools import style
# from mpltools import color
from Bio import SeqIO
import cPickle as pickle
import pandas as pd
import json

# Constants

# file paths
IN_DIR = "input_jsons"
FSEP = "/"
FP_MONO = "uby_1ubq.json"

# parsing mutations A X1Y, 1 is position, Y is mutation
POS_POS = 3
MUT_POS = 4

def json_to_dic(jfile, param):

    '''convert json file of Rosetta DDG values to dictionary'''

    json_data = open(jfile, 'r')
    raw_ddg = json.load(json_data)

    ddg_dic = {}

    for dat in raw_ddg['data']:

        allel = (int(dat['Mutation'][POS_POS]), str(dat['Mutation'][MUT_POS]))
        ddg_dic[allel] = dat[param]

    return ddg_dic


def run():
    # parameter value to extract
    param = "global_DDG"

    ddg_dic = json_to_dic(IN_DIR + FSEP + FP_MONO, param)

    pickle.dump(ddg_dic, open('ddg_dic.pkl', 'wb'))

run()