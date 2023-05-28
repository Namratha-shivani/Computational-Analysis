import pandas as pd
import numpy as np
import os

os.chdir(os.getcwd())

# using the count matrix calculate the weight matrix

counts_matrix = pd.read_csv('read the Matrix text file', sep = '\t', header = None, index_col=0)
counts_matrix.drop(1, inplace = True, axis = 1)
counts_matrix.columns = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]

Frequency_matrix = counts_matrix/sum(counts_matrix[1])

print('The Frequency Matrix is : ')
print('\n')
print(Frequency_matrix)
print('\n')

Pseudocount_matrix = (counts_matrix+1)/sum((counts_matrix+1)[1])

print('The Pseudocount Matrix is : ')
print('\n')
print(Pseudocount_matrix)
print('\n')

Weight_matrix = np.log2(Pseudocount_matrix / 0.25)
print('The Weight Matrix is : ')
print(Weight_matrix)

#Using the wieght matrix identify the top binding site

file = pd.read_csv('fasta file to find the motifs', header= None, sep= '\\')
# drop the rows with all NA's and replace the column names
file.dropna(axis = 1, how = 'all', inplace = True)
file.columns = ['gene_id','sequence']

def top_genes(data):
  position_matrix = Weight_matrix.transpose()
  W_len = position_matrix.shape[0]

  def calculate_weight_score(subseq):
    score = 0
    for i in range(len(subseq)):
      if subseq[i] in position_matrix.columns:
        score+= position_matrix[subseq[i]][i+1]
    return score

  def motif_finding(fasta):
    max_score = {}
    for index, rows in E_coli.iterrows():
      gene_id = rows['gene_id']
      sequence = rows['sequence'].replace(' ','')
      score = []
      for j in range(len(sequence)-W_len):
        sub_seq = sequence[j:(j+W_len)]
        weight_score = calculate_weight_score(sub_seq)
        score.append(weight_score)
      max_score[gene_id]  = max(score)
    return max_score

  return list(dict(sorted(motif_finding(data).items(),key = lambda item: item[1],reverse = True)).items())[:30]


print('Top 30 gene_ids : ')
for i in top_genes(file):
    print(i)
