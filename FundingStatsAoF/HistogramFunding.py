#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 09:55:12 2023

@author: alexanderjung
"""

import os, sys
import pandas as pd
import matplotlib.pyplot as plt

path = "/Users/alexanderjung/FundingStatsAoF"
dirs = os.listdir(path)
cntr=0 

for file in dirs:
    filename, file_extension = os.path.splitext(file)
    if file_extension==".csv": 
        if cntr== 0 :
            df1 = pd.read_csv(file,sep=";",thousands='.')
            cntr=1
        else : 
            df1.append(pd.read_csv(file,sep=";",thousands='.') )        
        print(filename)
        
df1 = df1[["Applicant / Contact person","Funding"]]        
histdata= df1.groupby(by=["Applicant / Contact person"], dropna=True).sum()
sorteddf= histdata.sort_values(by="Funding", ascending=False)

array = sorteddf["Funding"].to_numpy()

plt.plot(range(len(array)),array)
plt.xlabel("a researcher")
plt.ylabel("total funding")
plt.title("Funding granted during 2008 -2023")