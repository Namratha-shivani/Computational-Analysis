import pandas as pd
import numpy as np
import os
from sklearn.linear_model import LinearRegression

os.chdir(os.getcwd())
# data contains all the 3 time courses in a single file, first read the file
data = pd.read_csv('DecayTimecourse.txt',sep= '\t')

# Splitting the data and natural log transformation of the data
t1 = pd.DataFrame(data[data.columns[1:10]])
t1.set_index(data['Time course #'], inplace= True)
t1.dropna(how='all', inplace=True)
t1.set_axis([0,5,10,15,20,30,40,50, 60], axis = 'columns', inplace=True)
t1.drop('YORF', axis= 0, inplace=True)
t1 = np.log(t1)


t2 = pd.DataFrame(data[data.columns[10:19]])
t2.set_index(data['Time course #'], inplace= True)
t2.dropna(how='all', inplace=True)
t2.set_axis([0,5,10,15,20,30,40,50, 60], axis = 'columns', inplace=True)
t2.drop('YORF', axis= 0, inplace=True)
t2 = np.log(t2)

t3 = pd.DataFrame(data[data.columns[19:28]])
t3.set_index(data['Time course #'], inplace= True)
t3.dropna(how='all', inplace=True)
t3.set_axis([0,5,10,15,20,30,40,50, 60], axis = 'columns', inplace=True)
t3.drop('YORF', axis= 0, inplace=True)
t3 = np.log(t3)

# Regression fit to find the slope
def regression(a):
    lm = LinearRegression()
    slope = []
    for i in range(len(a)):
        a.iloc[i].replace([np.inf, -np.inf], np.nan, inplace=True)
        y = a.iloc[i].dropna(axis = 0)
        x = np.array(y.index).reshape(-1,1)
        lm.fit(x,y)
        slope.append(lm.coef_[0])
    return slope
        
slope_t1 = pd.DataFrame(regression(t1)).set_index(t1.index).rename(columns={0:'Slope'})
slope_t2 = pd.DataFrame(regression(t2)).set_index(t2.index).rename(columns={0:'Slope'})
slope_t3 = pd.DataFrame(regression(t3)).set_index(t3.index).rename(columns={0:'Slope'})
print(' ')
print('The slopes for each Transcrupt in the time course data are :')
print('TimeCourse1')
print(slope_t1)
print('TimeCourse2')
print(slope_t2)
print('TimeCourse3')
print(slope_t3)

# Calculating the half life
Half_life1 = np.log(2)/slope_t1.rename(columns={'Slope':'Half life'})
Half_life2 = np.log(2)/slope_t2.rename(columns={'Slope':'Half life'})
Half_life3 = np.log(2)/slope_t3.rename(columns={'Slope':'Half life'})
print(' ')
print('The Half lifes for each Transcrupt in the time course data are :')
print('Half lifes 1')
print(Half_life1)
print('Half lifes 2')
print(Half_life2)
print('Half lifes 3')
print(Half_life3)

# Join all the 3 time courses half lives to calculate the Average T 1/2 across all the data
me_12 = Half_life1.join(Half_life2,how = 'outer', lsuffix = 'Time course 1', rsuffix = 'Time course 2')
Half_life = me_12.join(Half_life3,how = 'outer', rsuffix = 'Time course 3')

Avg_Halflife = pd.DataFrame(Half_life.mean(skipna = True, axis = 1).sort_values(ascending=False))
Avg_Halflife.columns = ['Halflife']
print(' ')
print('The average Half-life over all the three time courses for each transcript :')
print(Avg_Halflife)
print(' ')


# considering only the top and bottom 10% of transcripts
Top_Halflifes = pd.concat([Avg_Halflife.iloc[0:int(((len(Avg_Halflife)*10)/100))],
                           Avg_Halflife.iloc[len(Avg_Halflife)-int(((len(Avg_Halflife)*10)/100)) : len(Avg_Halflife)]])

print('The top 10% and the bottom 10% ranked transcripts based on Half-lives : ')
print(Top_Halflifes)












































