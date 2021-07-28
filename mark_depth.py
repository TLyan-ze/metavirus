# -*- coding: utf-8 -*-
from __future__ import division
 
import csv
import sys
import numpy as np
import pandas as pd
from pandas import DataFrame
import openpyxl


merge1 = sys.argv[1]

df = pd.DataFrame(pd.read_csv(merge1,sep="\t",low_memory = False,header=None))
print (df)
nor=df.shape[0]
dict1 = {}
list_family=[]
list_samples=[]
for i in range(0, nor):
    #a_key= str(df.loc[i, 'Lab_sample_id'])
    if df.loc[i, 0] not in dict1.keys():
        dict1[df.loc[i, 0]] = str(df.loc[i, 2])
    elif df.loc[i, 0] in dict1.keys():
        dict1[df.loc[i, 0]] = dict1[df.loc[i, 0]]+","+str(df.loc[i, 2])

#print(dict1)
dict2={}
dict3={}
for key,value in dict1.items():
    value1=value.split(',')
    a=len(value1)
    p1=value1.count("0")
    p2=a - p1
    cover=p2/a
    value1 = list(map(int,value1))
    average_a = np.mean(value1)
    median_a = np.median(value1)
    dict2[key] = str(a)+","+str(average_a) + ","+str(median_a)

print(dict2)