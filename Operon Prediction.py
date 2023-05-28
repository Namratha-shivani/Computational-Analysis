import pandas as pd
import os 
import gzip
import shutil


os.chdir(os.getcwd())

# The function first unzips the files if needed and then calculates the operons with intergenic distance less than 50

def predictOperon(y):
    
    def fileread(x):
        
        def unzip(x):
            f_in = gzip.open(file, 'rb')
            f_out = open(file.split('.gz')[0], 'wb')
            shutil.copyfileobj(f_in, f_out)
            f_in.close()
            f_out.close()
            return (os.getcwd()+'/'+x.split('.gz')[0])
        
        if '.ptt' in x:
            if x in os.listdir() :
                return (pd.read_csv(x, sep = '\t', skiprows = 2))
            
            else:
                file = x +'.gz'
                return(pd.read_csv(unzip(file), sep = '\t', skiprows=2)) 
            
        else:
            return (pd.read_csv(x, sep='\t', comment='#', header=None, usecols=[0,3,4,6], names=['Gene', 'start', 'stop', 'Strand']))
        
    def operons(z):
        
        Operons = [] 
        operon = [] 
        i =0
        
        while i != len(file):
            if i ==0:
                operon.append(file['Gene'].iloc[i])
                i +=1
                
            else:
                prev_gene = file.iloc[i-1]
                gene = file.iloc[i]
                
                if prev_gene['Strand'] != gene['Strand'] or gene['start'] - prev_gene['stop'] >= 50 :
                    Operons.append(operon)
                    operon = [gene['Gene']]
                    i +=1  
                    
                else:
                    operon.append(gene['Gene'])
                    i +=1
                    
        Operons.append(operon)
        return (Operons,len(Operons))
  
    if '.ptt' in y:
        file = fileread(y)
        start, stop = [[],[]]
        
        for Location in file['Location']:
            start.append(int(Location.split('..')[0]))
            stop.append(int(Location.split('..')[1]))
            
        file['start'] = start
        file['stop'] = stop     
            
        file.sort_values(by = 'start', inplace = True)
        return operons(file)
    
    else:
        file = fileread(y)
        file.sort_values(by = ['Gene','start'], inplace = True)
        return operons(file)
    

