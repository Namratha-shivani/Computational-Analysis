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

# QUESTION 1:

matrix1 = pd.read_csv('matrix1.txt',sep='\t')
matrix1.set_index('miRNA',inplace = True)
matrix1.drop('Unnamed: 13',axis = 1, inplace = True)

corr_matrix1 = pd.DataFrame(np.array(Pearson_Correlation(matrix1)), index = matrix1.columns , columns = matrix1.columns)
print("")
print('Pearson Correlated 12 X 12 matrix of matrix1')
print(corr_matrix1)


# QUESTION 2:



matrix2 = pd.read_csv("matrix2.txt",sep = '\t')
matrix2.set_index('miRNA',inplace = True)
matrix2.drop('Unnamed: 13',axis = 1, inplace = True)

corr_matrix2 = pd.DataFrame(np.array(Pearson_Correlation(matrix2)), index = matrix2.columns, columns = matrix2.columns)
print("")
print('Pearson Correlated 12 X 12 matrix of matrix2')
print(corr_matrix2)



mat1_2corr = pd.DataFrame(np.array(Pearson_Correlation(matrix1,matrix2)), index = matrix1.columns, columns = matrix2.columns)
print("")
print('Pearson Correlated 12 X 12 matrix of matrix1 and matrix2')
print(mat1_2corr)


f,(ax1,ax2,ax3) = plt.subplots(1,3,sharey=True, figsize= (15,5))
g1 = sns.heatmap(corr_matrix1, ax= ax1).set(title = 'Matrix1 Corr' )
g2 = sns.heatmap(corr_matrix2, ax = ax2).set(title = 'Matrix2 Corr' )
g3 = sns.heatmap(mat1_2corr, ax = ax3).set(title = 'Matrix1 & Matrix2 Corr' )

plt.show()

