import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns


os.chdir(os.getcwd())

def Pearson_Correlation(x, y=None):
    
  def mean(x):
    sum = 0
    for a in range(len(x)):
      sum += x[a]
    return sum/len(x)
  
  def adj(c):
    adj = []
    Mean = mean(c)
    for i in c:
      adj.append((i-Mean))
    return(adj)

  def sum(d):
    add = 0
    t = adj(d)
    for i in t:
      add+=i**2
    return add

  def sum_xy(e,f):
    X = adj(e)
    Y = adj(f)
    add_xy = 0
    for i in range(len(X)):
      add_xy += (X[i]*Y[i])
    return add_xy

  def pearson_corr(x,y):
    pearsoncorr =[]
    if type(y) != pd.core.frame.DataFrame:
      for i in x.columns:
        corr = [] 
        for j in x.columns: 
          a = sum_xy(x[i],x[j])/((sum(x[i])*sum(x[j]))**0.5)
          corr.append(a)
        pearsoncorr.append(corr)
    else:
      if len(x.columns) == len(y.columns):
        for i in x.columns:
          corr = [] 
          for j in y.columns: 
            a = sum_xy(x[i],y[j])/((sum(x[i])*sum(y[j]))**0.5)
            corr.append(a)
          pearsoncorr.append(corr)
      else:
        raise ValueError('Both the matrix have different column length')
    return pearsoncorr
  return pearson_corr(x,y)
