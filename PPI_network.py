#installation of igraph
pip install igraph

import pandas as pd
import numpy as np
import os
import igraph
import matplotlib.pyplot as plt
from scipy.stats import ranksums


os.chdir(os.getcwd())
# read the files 

# Humann-PPI.txt file contains rows of proteins one row of unique proteins with their connections in the other rows
# Protein-list1 and Protein-list2.txt contain list of proteins to find the shortest paths among them

PPI = pd.read_table('Human-PPI.txt')
p1 = pd.read_csv('protein-list1.txt', header = None) 
p2 = pd.read_csv('protein-list2.txt',header = None) 

# change the name of the columns and convert the network to dataframe
PPI_network = {'Protein_A': PPI.index, 'Protein_B': PPI.columns[0]}
PPI_network = pd.DataFrame(PPI_network)
PPI_network.index = [i for i in range(len(PPI_network))]



# Degree Calculation

degree = {}
for i in PPI_network.columns:
  for j in range(len(PPI_network[i])):
    if PPI_network[i][j] not in degree:
      degree[PPI_network[i][j]] = 1
    else:
      degree[PPI_network[i][j]]+=1  

Degree = (pd.DataFrame(degree, index = ['Degree']).transpose())

print('')
print('The degrees for each protein: ')
print(Degree)

# connections dictionary

connections1 = {}
for name in Degree.index:
  selected_rows = (PPI_network[(PPI_network['Official_symbol_A']==name)| ( PPI_network['Official_symbol_B'] == name)])
  other_values = [x for x in list(selected_rows['Official_symbol_A']) + list(selected_rows['Official_symbol_B']) if x != name]
  connections1[name] = other_values


# Clustering Co-efficient

# M =number of connections seen
# N = total nodes connected to the node of interest
# clustering Co-efficient = (2*M)/ (N*(N-1))

Clustering_coef= {}
for i in connections1.keys():
    M = 0
    N = len(connections1[i])
    for j in range(N-1):
      for k in range(j+1, N):
          if connections1[i][k] in connections1[connections1[i][j]]:
             M+=1
          else:
             continue
    if N*(N-1) == 0:
      Clustering_coef[i] = 0
    else:
      Clustering_coef[i] = (2*M)/(N*(N-1))    
  
  

clustering_coef = pd.DataFrame(Clustering_coef, index = ['clustering coef']).transpose()

print('Clustering Coefficients for each protein in the network: ')
print(clustering_coef)


#Average clustering coefficient

print('')
print('The average clustering coefficient is : ',np.mean(clustering_coef['clustering coef']))


# Scale Free Network

degree_freq = np.bincount(Degree['Degree'])
degree_freq = degree_freq[1:]

plt.loglog(range(len(degree_freq)), degree_freq, 'o', markersize=8,color = 'black')
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.title('Degree distribution of the network')
plt.show()


# Shortest path
graph = igraph.Graph.TupleList(PPI_network.itertuples(index=False), directed=False, vertex_name_attr='name')

def shortest_path_calculation(p, adj_matrix):
    shortest_path = []
    graph = igraph.Graph.Adjacency((adj_matrix.values > 0).tolist())
    graph.vs['name'] = adj_matrix.index.tolist()
    for i in range(len(p[0])-1):
        for j in range(i+1, len(p[0])):
            protein1 = p[0][i]
            protein2 = p[0][j]
            if protein1 in adj_matrix.index and protein2 in adj_matrix.index:
                start = adj_matrix.index.get_loc(protein1)
                end = adj_matrix.index.get_loc(protein2)
                path = graph.get_shortest_paths(start, to=end, output='vpath')
                if path:
                    shortest_path.append(len(path[0])-1)
            else:
                pass
    return(shortest_path)

adj_df = pd.DataFrame(np.array(graph.get_adjacency().data), index = graph.vs['name'], columns=graph.vs['name'])
path_list1 = shortest_path_calculation(p1, adj_df)
path_list2 = shortest_path_calculation(p2, adj_df)

wilcox = ranksums(path_list1,path_list2)
print('')
print('Wilcox Results: ')
print('statistic =', wilcox.statistic)
print('Pvalue =', wilcox.pvalue)
