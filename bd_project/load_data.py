import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from Bio.Align import substitution_matrices
from Bio.Data import IUPACData
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from pathlib import Path
import os
import re
import ast
import h5py

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def load_embeddings(embeddings_path):
    data_list = []

    with h5py.File(embeddings_path, "r") as f:
        for dataset_name in f.keys():
            dataset = f[dataset_name][:]
            data_list.append([dataset_name, dataset])

    return pd.DataFrame(data_list, columns=["ID", "embeddings"])


def load_dataset(training_data_path, test_data_path, sub_ontology):

    # Train location
    train_set_path = training_data_path / 'train_set.tsv'
    train_ids_path = training_data_path / 'train_ids.txt'
    train_fasta_path = training_data_path / 'train.fasta'
    train_embeddings_path = training_data_path / 'train_embeddings.h5'
    train_protein2ipr_path = training_data_path / 'train_protein2ipr.dat'
    go_basic_path = training_data_path / 'go-basic.obo'

    # Test location
    test_ids_path = test_data_path / 'test_ids.txt'
    test_fasta_path = test_data_path / 'test.fasta'
    test_embeddings_path = test_data_path / 'test_embeddings.h5'
    test_protein2ipr_path = test_data_path / 'test_protein2ipr.dat'






#                                  ######################
#    ~~~~~~~~~~~~~~~~~~~~~~~~~~    #    TRAINING SET    #   ~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  ######################





#   ••••••••••••••••••••••••••••••••••••
#   •    Extracting `train_set.tsv`    •  
#   ••••••••••••••••••••••••••••••••••••

    if not train_set_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{train_set_path}' doesn't exist.")
        return None 
    
    train_set = pd.read_csv(train_set_path, sep='\t')

    # Rename Protein_ID and aspect columns
    train_set.rename(columns={'Protein_ID': 'ID', 'aspect' : 'sub_ontology'}, inplace=True)




#   ••••••••••••••••••••••••••••••••••••
#   •    Extracting `train_ids.txt`    •  
#   ••••••••••••••••••••••••••••••••••••

    if not train_ids_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{train_ids_path}' doesn't exist.")
        return None
    
    # Extracting train_ids.txt
    with open(training_data_path / 'train_ids.txt', 'r') as file:
        train_ids = file.read().splitlines()
    



#   ••••••••••••••••••••••••••••••••••
#   •    Extracting `train.fasta`    •  
#   ••••••••••••••••••••••••••••••••••    

    if not train_fasta_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{train_fasta_path}' doesn't exist.")
        return None
    
    # If the file exists, read the FASTA file
    train_fasta_list = list(SeqIO.parse(train_fasta_path, 'fasta')) 
    

    # Extract relevant information from SeqRecord
    train_fasta_dict = [{
        'ID': record.id,
        'name': record.name,
        'description': record.description,
        'num_features': len(record.features),
        'sequence': record.seq,
    } for record in train_fasta_list]

    # Create a DataFrame from the extracted data
    train_fasta = pd.DataFrame(train_fasta_dict)


    # Removing 'name', 'description', and 'num_features' columns
    dropped = train_fasta.drop(columns=['name', 'description', 'num_features'], inplace=True)




#   ••••••••••••••••••••••••••••••••••••••••••
#   •    Extracting `train_embeddings.h5`    •  
#   ••••••••••••••••••••••••••••••••••••••••••
    
    if not train_embeddings_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{train_embeddings_path}' doesn't exist.")
        return None
    

    data_list = []

    with h5py.File(train_embeddings_path, "r") as f:
        for dataset_name in f.keys():
            dataset = f[dataset_name][:]
            data_list.append([dataset_name, dataset])

    train_embeddings = pd.DataFrame(data_list, columns=["ID", "embeddings"])





#   ••••••••••••••••••••••••••••••••••••••••••••
#   •    Extracting `train_protein2ipr.dat`    •  
#   ••••••••••••••••••••••••••••••••••••••••••••

    if  not train_protein2ipr_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{train_set_path}' doesn't exist.")
        return None


    train_protein2ipr = pd.read_csv(train_protein2ipr_path, sep='\t')

    # Rename Protein_ID and aspect columns
    train_protein2ipr.columns = ['ID', 'ipr', 'domain', 'familyID', 'start', 'end']

    # Group by 'ID' and aggregate other columns into lists
    train_protein2ipr_grouped = train_protein2ipr.groupby('ID').agg(lambda x: tuple(x)).reset_index()




#   •••••••••••••••••••••••••••••••••••
#   •    Extracting `go-basic.obo`    •
#   •••••••••••••••••••••••••••••••••••

    
    if  not go_basic_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{go_basic_path}' doesn't exist.")
        return None
    

    # Step 1: Initialize storage for GO terms
    go_terms = []

    # Step 2: Parse the .obo file
    with open(go_basic_path, 'r') as file:
        current_term = {}
        for line in file:
            line = line.strip()
            
            # Start of a new term
            if line == "[Term]":
                if current_term:  # Save the previous term
                    go_terms.append(current_term)
                current_term = {}  # Start a new term
                
            elif line.startswith("id:"):
                current_term['ID'] = line.split("id: ")[1]
                
            elif line.startswith("alt_id:"):
                alt_id = line.split("alt_id: ")[1]
                current_term.setdefault('alt_ids', []).append(alt_id)
                
            elif line.startswith("name:"):
                current_term['name'] = line.split("name: ")[1]
                
            elif line.startswith("namespace:"):
                current_term['namespace'] = line.split("namespace: ")[1]
                
            elif line.startswith("is_a"):
                match = re.search(r"GO:\d+", line)  # Search for GO ID
                if match:  # Check if a match was found
                    is_a_id = match.group()
                    current_term.setdefault('is_a', []).append(is_a_id)
                    
            elif line.startswith("relationship: part_of"):
                match = re.search(r"GO:\d+", line)  # Search for GO ID
                if match:  # Check if a match was found
                    part_of_id = match.group()
                    current_term.setdefault('part_of', []).append(part_of_id)

                
        # Add the last term
        if current_term:
            go_terms.append(current_term)

    # Step 3: Create a unified list of all IDs (primary and alt_ids)
    expanded_terms = []
    for term in go_terms:
        primary_id = term['ID']
        alt_ids = term.get('alt_ids', [])
        all_ids = [primary_id] + alt_ids
        
        for term_id in all_ids:
            expanded_terms.append({
                'ID': term_id,
                'name': term.get('name'),
                'namespace': term.get('namespace'),
                'is_a': term.get('is_a', []),
                'part_of': term.get('part_of', [])
            })

    # Step 4: Convert to a DataFrame
    #  = pd.DataFrame(expanded_terms)


#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Merging the training set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    def group_and_combine(df, sub_ontology_value):
        return df[df['sub_ontology'] == sub_ontology_value].groupby('ID')['GO_term'].apply(tuple).reset_index()

    
    # Extract the relevant sub-ontology

    df_training = group_and_combine(train_set, sub_ontology)

    
    # Merge the DataFrames

    # combined_train = pd.merge(train_embeddings, train_fasta, on='ID')
    # combined_train = pd.merge(combined_train, train_protein2ipr_grouped, on='ID', how='left')

    df_training_full= pd.merge(train_embeddings, df_training, on='ID', how='right')


    X_df = df_training_full.iloc[:, :-1]
    y_df = df_training_full.iloc[:, -1]




#                                       ##################
#    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    #    TEST SET    #   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       ##################




#   •••••••••••••••••••••••••••••••••••
#   •    Extracting `test_ids.txt`    •  
#   •••••••••••••••••••••••••••••••••••

    if not test_ids_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{test_ids_path}' doesn't exist.")
        return None


    with open(test_ids_path, 'r') as file:
        test_ids = file.read().splitlines()

    

    
#   •••••••••••••••••••••••••••••••••
#   •    Extracting `test.fasta`    •  
#   •••••••••••••••••••••••••••••••••

    if not test_fasta_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{test_fasta_path}' doesn't exist.")
        return None


    test_fasta_list = list(SeqIO.parse(test_fasta_path, 'fasta'))


    # Extract relevant information from SeqRecord
    test_fasta_dict = [{
        'ID': record.id,
        'name': record.name,
        'description': record.description,
        'num_features': len(record.features),
        'sequence': record.seq,
    } for record in test_fasta_list]

    # Create a DataFrame from the extracted data
    test_fasta = pd.DataFrame(test_fasta_dict)

    # Removing 'name', 'description', and 'num_features' columns
    dropped = test_fasta.drop(columns=['name', 'description', 'num_features'], inplace=True)




#   •••••••••••••••••••••••••••••••••••••••••
#   •    Extracting `test_embeddings.h5`    •  
#   •••••••••••••••••••••••••••••••••••••••••

    if not test_embeddings_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{test_embeddings_path}' doesn't exist.")
        return None
    

    data_list = []

    with h5py.File(test_embeddings_path, "r") as f:
        for dataset_name in f.keys():
            dataset = f[dataset_name][:]
            data_list.append([dataset_name, dataset])

    test_embeddings = pd.DataFrame(data_list, columns=["ID", "embeddings"])




#   •••••••••••••••••••••••••••••••••••••••••••
#   •    Extracting `test_protein2ipr.dat`    •  
#   •••••••••••••••••••••••••••••••••••••••••••

    if not test_protein2ipr_path.exists():
        # If the file doesn't exist, return None
        print(f"Error: File path '{test_protein2ipr_path}' doesn't exist.")
        return None


    test_protein2ipr = pd.read_csv(test_protein2ipr_path, sep='\t')

    # Rename Protein_ID and aspect columns
    test_protein2ipr.columns = ['ID', 'ipr', 'domain', 'familyID', 'start', 'end']

    # Remove 'domain' 
    dropped = test_protein2ipr.drop('domain', axis=1)

    # Group by 'ID' and aggregate other columns into lists
    test_protein2ipr_grouped = test_protein2ipr.groupby('ID').agg(lambda x: tuple(x)).reset_index()




#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Merging the test set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 



    combined_test = pd.merge(test_embeddings, test_fasta, on='ID')
    df_test_full = pd.merge(combined_test, test_protein2ipr_grouped, on='ID', how='left')


    X_test_df = df_test_full.iloc[:, :-1]
    y_test_df = df_test_full.iloc[:, -1]






#                                  ############################
#    ~~~~~~~~~~~~~~~~~~~~~~~~~~    #    PCA for embeddings    #   ~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  ############################



    def pca_col(df, name, column_name='embeddings', variance_threshold=0.90):
        embeddings = np.array(df[column_name].tolist())
        
        scaler = StandardScaler()
        embeddings_standardized = scaler.fit_transform(embeddings)
        
        pca = PCA()
        pca.fit(embeddings_standardized)
        explained_variance = pca.explained_variance_ratio_

        # Components to retain (in our case 90% of variance)
        n_components = np.argmax(np.cumsum(explained_variance) >= variance_threshold) + 1
            
        # plt.figure(figsize=(8, 4.5))
        # plt.axhline(y=variance_threshold, color='r', linestyle='--')
        # x_value = np.argmax(np.cumsum(explained_variance) >= variance_threshold) + 1
        # plt.axvline(x=x_value, color='g', linestyle='--')
        # plt.text(x_value + 10, 0.025, f'{x_value}', color='green', size='large', weight='bold')
        # plt.text(- 50, variance_threshold + 0.025, f'{variance_threshold:.2f}', color='red', size='large', weight='bold')

        # plt.plot(np.cumsum(explained_variance), marker='o')
        # plt.title(f'Cumulative Explained Variance by PCA Components for {name}')
        # plt.xlabel('Number of Components')
        # plt.ylabel('Cumulative Explained Variance')
        # plt.grid(True)
        # plt.show()
        
        pca = PCA(n_components=n_components)
        reduced_embeddings = pca.fit_transform(embeddings_standardized)
        
        # Integrate reduced embeddings into the DataFrame
        df[f'reduced_{column_name}'] = reduced_embeddings.tolist()
        return df



    # Apply PCA and add column for each dataset
    X_df = pca_col(X_df, sub_ontology)
    #x_df_CC = pca_col(x_df_CC, 'CC')
    #x_df_MF = pca_col(x_df_MF, 'MF')
    ##### Save the datasets
    # Ensure the directory exists
    datasets_path = Path('../data/datasets')
    datasets_path.mkdir(parents=True, exist_ok=True)


    return X_df, y_df, X_test_df, y_test_df 