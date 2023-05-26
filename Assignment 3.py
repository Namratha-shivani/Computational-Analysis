import pandas as pd
import os 
import gzip
import shutil


os.chdir(os.getcwd())

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
    
print('QUESTION 1')
print('')  
print('E.COLI : ')
print('')
print("Operons in E.coli are:\n\n{}\n\nNumber of Operons in E.coli are: {}".format(predictOperon('E_coli_K12_MG1655.ptt')[0], predictOperon('E_coli_K12_MG1655.ptt')[1])) #2669
print('')
print("HALOBACTERIUM : ")
print('')  
print("Operons in Halobacterium are:\n\n{}\n\nNumber of Operons in Halobacterium are: {}".format(predictOperon('Halobacterium_NRC1.ptt')[0], predictOperon('Halobacterium_NRC1.ptt')[1])) #1464
print('')
print('SYNECHOCYSTIS : ')
print('')  
print("Operons in Synechocystis are:\n\n{}\n\nNumber of Operons in Synechocystis are: {}".format(predictOperon('Synechocystis_PCC6803_uid159873.ptt')[0], predictOperon('Synechocystis_PCC6803_uid159873.ptt')[1])) #2521
print('')
print('B_SUBTILIS : ')
print('')  
print("Operons in B_subtilis are:\n\n{}\n\nNumber of Operons in B_subtilis are: {}".format(predictOperon('B_subtilis_168.ptt')[0], predictOperon('B_subtilis_168.ptt')[1])) #2662  
print('')  
print('QUESTION 2')
print('')
print('CROP_MICROBIOME : ')
print('')  
print("Operons in Crop_Microbiome are:\n\n{}\n\nNumber of Operons in Crop_Microbiome are: {}".format(predictOperon('2088090036.gff')[0], predictOperon('2088090036.gff')[1])) #12108
 
    
    
    
    
    
    
    
    
    
    
    
    