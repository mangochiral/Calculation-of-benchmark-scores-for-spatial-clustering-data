# Modules required --------------------------------------------------
import os
import pandas as pd
import numpy as np
from sklearn.metrics.cluster import normalized_mutual_info_score, rand_score,fowlkes_mallows_score,adjusted_rand_score,davies_bouldin_score

# Folder path to the datasets --------------------------------------------

file_path = os.getcwd()
print(file_path)
path = input("Enter file path for counts, labels, and spatial_clustering_results csv files: "):
if file_path != path
    file_path = os.chdir(path)
file_path = os.getcwd()
print(file_path)
# List of all files and directories in a give folder
dest = os.listdir(file_path)

# List of all files and directories where clustering and labels are kept
files = []
for i in dest:
    if i.startswith('V'):
        files.append(i)
print(files)


class convert:
    def __init__(self,diction,wh):
        for i in files:
            # An empty dictionary to store filenames
            tmp = {i: 0}

            # Path to clustering file and saving it as pandas dataframe
            cluster = os.path.join(i, 'spatial_clustering_results.csv')
            spa = pd.read_csv(cluster)
            spa_pd = pd.DataFrame(spa)
            # Selecting cluster column
            spa_pd = spa_pd.loc[:, "Cluster"]

            # Path to labels file and saving it as pandas dataframe
            label = os.path.join(i, 'labels.csv')
            name = pd.read_csv(label)
            name_pd = pd.DataFrame(name)

            # Selecting class column
            name_pd = name_pd.loc[:, "class"]

            # This block to make arguments for DBI scores------------------------------------------

            # Path to counts file which contains genes and spots and saving it as pandas dataframe
            X = os.path.join(i, 'counts.csv')
            X = pd.read_csv(X)
            # Transposing the X dataframe for DBI Scoring
            X = X.transpose()
            # Converting the cluster information into labels for DBI scoring
            labels = np.array(spa_pd)

            # This block to calculate the scores----------------------------------------------
            # All the scores are stored in tmp dictionary and are
            # further updated to final dictionary diction arugment
            if wh == 'nmi' or wh == 'NMI':
                tmp[i] = normalized_mutual_info_score(spa_pd, name_pd)
                diction.update(tmp)
            elif wh == 'rand' or wh == 'RAND':
                tmp[i] = rand_score(spa_pd, name_pd)
                diction.update(tmp)
            elif wh == 'fow' or wh == 'FOW':
                tmp[i] = fowlkes_mallows_score(spa_pd, name_pd)
                diction.update(tmp)
            elif wh == 'ari' or wh == 'ARI':
                tmp[i] = adjusted_rand_score(spa_pd, name_pd)
                diction.update(tmp)
            elif wh == 'dbi':
                tmp[i] = davies_bouldin_score(X, labels)
                diction.update(tmp)
        # Converting the diction dictionary to pandas dataframe
        self.score = pd.DataFrame(diction, index=[0])
        # self.rand = args[1]
        # self.fow = args[2]
    def cal(self):
            # Returning diction data frame
            return self.score

# For Saving the data in csv file format
file_gen = input("enter the scores you want to calculate with space in between(ari/dbi/rand/fow/nmi): ").split()
for t in file_gen:
    if t == 'ari' or t == 'ARI':
        data = convert(diction= {}, wh = t)
        data = pd.melt(data.cal(),var_name='Dataset', value_name='ARI', ignore_index=True)
        if (not os.path.exists('BayesSpaceARI.csv')):
            data.to_csv('BayesSpacesARI.csv')
    elif t == 'fow' or t == 'FOW':
        data = convert(diction= {}, wh = t)
        data = pd.melt(data.cal(),var_name='Dataset', value_name='FOWL', ignore_index=True)
        if (not os.path.exists('BayesSpaceFOWL.csv')):
            data.to_csv('BayesSpaceFOWL.csv')
    elif t == 'nmi' or t == 'NMI':
        data = convert(diction= {}, wh = t)
        data = pd.melt(data.cal(),var_name='Dataset', value_name='NMI', ignore_index=True)
        if (not os.path.exists('BayesSpaceNMI.csv')):
            data.to_csv('BayesSpaceNMI.csv')
    elif t == 'rand' or t == 'RAND':
        data = convert(diction= {}, wh = t)
        data = pd.melt(data.cal(),var_name='Dataset', value_name='RAND', ignore_index=True)
        if (not os.path.exists('BayesSpaceRAND.csv')):
            data.to_csv('BayesSpaceRAND.csv')
    elif t == 'dbi':
        data = convert(diction= {}, wh = t)
        data = pd.melt(data.cal(),var_name='Dataset', value_name='DBI', ignore_index=True)
        if (not os.path.exists('BayesSpaceDBI.csv')):
            data.to_csv('BayesSpaceDBI.csv')

