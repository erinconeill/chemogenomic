#!/usr/bin/env python
# coding: utf-8

# In[74]:


import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem import AllChem
from rdkit import DataStructs


# In[75]:


# Import Data 
fda_approved = pd.read_csv('FDA_Approved.csv', header = None, names = ['ID', 'Drug Name']) # 2331 rows
dc_compounds = pd.read_csv('DC_Compounds.csv') # 4099 rows
dc_drug_tar = pd.read_csv('Drug_Target.csv') # 19378 rows

# Merge files on 'ID' columns with 'left' join
merged_df = pd.merge(dc_compounds, fda_approved, on = 'ID', how = 'left') 

# Create column 'Approved' indicating whether the drug has FDA approval or not (0 or 1) by checking for ID in fda_approved
# Retain results without FDA approval
merged_df['Approved'] = merged_df['ID'].isin(fda_approved['ID']).replace({True: 1, False: 0})
merged_df.loc[merged_df['Approved'] != 1, 'Approved'] = 0

# Compare merged_df with Drug_Target CSV, create new resulting df 
drug_central = pd.merge(dc_drug_tar, merged_df, left_on = 'STRUCT_ID', right_on = 'ID', how = 'left')
drug_central = drug_central.drop(columns = ['Drug Name', 'INN'])

# Count the occurrences of '1' in the 'Approved' column
counts = merged_df['Approved'].value_counts().get(1, 0)

# Display and export
display(drug_central)
drug_central.to_csv('drug_central.csv', index = False)
print(f"'1' appears in the 'Approved' column: {counts}") # FDA approved


# ## SGC Donated Chemical Probes

# In[76]:


sgc_compounds = pd.read_csv('SGC_Compounds.csv', skiprows = 1) #123 rows x 24 columns
sgc_compounds = sgc_compounds.dropna(subset=['SMILES (unique cis trans)']) # none to drop
sgc_compounds


# ## Chemical Probes Portal

# In[77]:


chem_probes = pd.read_csv('ChemicalProbesPortal.csv') 
chem_probes


# In[78]:


def smiles_to_inchi(smiles):
    if isinstance(smiles, float):
        return None

    mol = Chem.MolFromSmiles(str(smiles))
    
    if mol is not None:
        inchi_string = inchi.MolToInchi(mol)
        return inchi_string
    else:
        return None


# In[79]:


# Create 'InChI' column using the 'SMILES' column
drug_central['InChI (cleaned)'] = drug_central['SMILES'].apply(smiles_to_inchi)

# Display
display(drug_central)


# In[80]:


# Create 'InChI' column using the 'SMILES (unique cis trans)' column
sgc_compounds['InChI (cleaned)'] = sgc_compounds['SMILES (unique cis trans)'].apply(smiles_to_inchi)

# Display
display(sgc_compounds)


# In[81]:


chem_probes['InChI (cleaned)'] = chem_probes['SMILES'].apply(smiles_to_inchi)

# Display
display(chem_probes)
chem_probes.to_csv('chem_probes.csv')


# In[82]:


drug_central_inchi = drug_central[['InChI (cleaned)', 'DRUG_NAME']]
sgc_inchi = sgc_compounds[['InChI (cleaned)', 'Compound name']]
chem_probes_inchi = chem_probes[['InChI (cleaned)', 'Probe name']]

display(drug_central_inchi)
display(sgc_inchi)
display(chem_probes_inchi)


# ## Concatination to determine compound overlap

# In[83]:


df = pd.concat([drug_central_inchi, sgc_inchi, chem_probes_inchi], ignore_index=True)

# Rename the columns
df.rename(columns={'DRUG_NAME': 'Drug Central', 'Compound name': 'SGC', 'Probe name': 'Chemical Probes'}, inplace=True)
df


# In[84]:


melted_df = pd.melt(df, id_vars=['InChI (cleaned)'], var_name='Source', value_name='Compound')
melted_df


# In[85]:


# Groupby and aggregate
first_result = melted_df.groupby('InChI (cleaned)', as_index=False).agg({'Source': ', '.join, 'Compound': 'first'}) #3271 rows
# combines source values within each group (same inchi)
first_result


# In[86]:


# List result
list_result = melted_df.groupby('InChI (cleaned)', as_index=False).agg({'Source': ', '.join, 'Compound': list}) # list is a more optimal structure
list_result


# In[87]:


def remove_nan_sources(row):
    sources = row['Source'].split(', ')
    compounds = row['Compound']
    
    updated_sources = [source for source, compound in zip(sources, compounds) if not pd.isna(compound)]
    
    return ', '.join(updated_sources)

# Apply to update 'Source' column
list_result['Source'] = list_result.apply(remove_nan_sources, axis=1)

# Drop NaN values from the 'Compound' column
list_result['Compound'] = list_result['Compound'].apply(lambda x: [comp for comp in x if not pd.isna(comp)])

# Display
display(list_result)

#3,271


# In[88]:


# Create copies to avoid SettingWithCopyWarning
drug_central_inchi_copy = drug_central_inchi.copy()
sgc_inchi_copy = sgc_inchi.copy()
chem_probes_inchi_copy = chem_probes_inchi.copy()

# Remove duplicates based on 'InChI' column for drug_central_inchi_copy
drug_central_inchi_no_dupes = drug_central_inchi_copy.drop_duplicates(subset=['InChI (cleaned)'])

# Remove duplicates based on 'InChI' column for sgc_inchi_copy
sgc_inchi_no_dupes = sgc_inchi_copy.drop_duplicates(subset=['InChI (cleaned)'])

# Remove duplicates based on 'InChI' column for chem_probes_inchi_copy
chem_probes_inchi_no_dupes = chem_probes_inchi_copy.drop_duplicates(subset=['InChI (cleaned)'])

# Display the DataFrames without duplicates based on 'InChI' column
print("drug_central_inchi without duplicates:")
display(drug_central_inchi_no_dupes)

print("\nsgc_inchi without duplicates:")
display(sgc_inchi_no_dupes)

print("\nchem_probes_inchi without duplicates:")
display(chem_probes_inchi_no_dupes)


# In[89]:


df = pd.concat([drug_central_inchi_no_dupes, sgc_inchi_no_dupes, chem_probes_inchi_no_dupes], ignore_index=True)

# Rename the columns
df.rename(columns={'DRUG_NAME': 'Drug Central', 'Compound name': 'SGC', 'Probe name': 'Chemical Probes'}, inplace=True)
df


# In[90]:


result_df = melted_df.groupby('InChI (cleaned)', as_index=False).agg({'Source': ', '.join, 'Compound': list})

def remove_nan_sources(row):
    sources = row['Source'].split(', ')
    compounds = row['Compound']
    
    updated_sources = [source for source, 
                       compound in zip(sources, compounds)
                       if not pd.isna(compound)]
    
    return ', '.join(updated_sources)

# Apply to update 'Source' column
result_df['Source'] = result_df.apply(remove_nan_sources, axis=1)

# Drop NaN values from the 'Compound' column
result_df['Compound'] = result_df['Compound'].apply(lambda x: [comp for comp in x if not pd.isna(comp)])

# Convert 'Compound' column to string format
result_df['Compound'] = result_df['Compound'].apply(lambda x: ', '.join(x))

# Reorganize columns for aethetics
result_df = result_df[['InChI (cleaned)', 'Compound', 'Source']]
all_inchi = result_df[['InChI (cleaned)']]

# Select ONLY the overlapping InChI
selected_inchi = result_df[(result_df['Source'] == 'Drug Central, Chemical Probes') | (result_df['Source'] == 'SGC, Chemical Probes')]
selected_inchi = selected_inchi[['InChI (cleaned)']]

# Display
display(result_df)
display(all_inchi)
display(selected_inchi)
result_df.to_csv('result_df.csv')


# In[91]:


def get_inchikey(inchi_string):
    try:
        mol = Chem.MolFromInchi(inchi_string)
        inchikey = Chem.InchiToInchiKey(inchi_string)
        return inchikey
    except (TypeError, ValueError):
        return None

def inchi_to_inchikey(inchi_dataframe):
    inchi_copy = inchi_dataframe.copy() 
    inchi_copy['InChIKey'] = inchi_copy['InChI (cleaned)'].apply(get_inchikey)
    return inchi_copy

all_inchikey = inchi_to_inchikey(all_inchi)
all_inchikey = all_inchikey[['InChIKey']]

selected_inchikey = inchi_to_inchikey(selected_inchi)
selected_inchikey = selected_inchikey[['InChIKey']]
adjusted_inchikey = selected_inchikey.copy()
adjusted_inchikey['InChIKey']= 'InChIKey=' + adjusted_inchikey['InChIKey']

display(all_inchikey) 
display(selected_inchikey)
display(adjusted_inchikey)
selected_inchi.to_csv('selected_inchi.csv')

#174


# In[92]:


result_df


# In[93]:


compounds_cv = pd.read_csv("compounds_cv.csv", encoding='utf-8')
compounds_cv


# In[94]:


eubopen_merge = pd.merge(compounds_cv, result_df, left_on='Compound InChi', right_on='InChI (cleaned)')
eubopen_merge


# In[101]:


def get_similar_compounds(smiles, dataset, top_n=5):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

    similarities = []
    for _, row in dataset.iterrows():
        mol_other = Chem.MolFromSmiles(row['Compound SMILES'])
        if mol_other:
            fp_other = AllChem.GetMorganFingerprintAsBitVect(mol_other, 2, nBits=1024)
            similarity = DataStructs.TanimotoSimilarity(fp, fp_other)
            similarities.append((row['Virtual Compound Preferred Name'], row['Compound SMILES'], similarity))

    similarities.sort(key=lambda x: x[2], reverse=True)
    top_similar = similarities[:top_n]

    return top_similar

# Example usage
smiles_input = "CCO"
similar_compounds = get_similar_compounds(smiles_input, eubopen_merge)
for name, smiles, similarity in similar_compounds:
    print(f"Name: {name}, SMILES: {smiles}, Similarity: {similarity}")


# In[102]:


eubopen_merge

