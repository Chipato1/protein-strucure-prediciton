#!/usr/bin/env python
# coding: utf-8

### More Advanced Stuff
#### Author: Carl Winkler


import pandas as pd
import itertools

def parse_sequences(fileName):
    colums_pred = ['ID','length','Sequence','Target']
    data = []
    
    with open(fileName) as f:
        lines = f.readlines() # list containing lines of file
        i = 0
        A_row = []

        for line in lines:
            line = line.strip()
            elems = line.split('\t')
            if i == 0:
                # ID
                A_row.append(elems[0])
            elif i == 1:
                #length
                A_row.append(int(elems[0]))
            elif i == 2:
                # Whole list of tokens - sequence
                A_row.append(elems)
            elif i == 3:
                # Whole list of tokens - target
                A_row.append(elems)
                data.append(A_row) 

            #Logic for parsing every 4 things as a row
            i = i + 1
            if i > 3:
                i = 0
                A_row = []
            
    return pd.DataFrame(data, columns = colums_pred) 


# In[75]:


def parse_triplet_frequency(filename):
    # First we get the dataframe from the file parser
    data_df = parse_sequences(filename)
    
    # We add '-' here for the borders, Meaning the '---' combination will be created which is never used
    acid_types = ['-','A', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L', 'M', 'N', 'P', 'Q','R', 'S', 'T', 'V', 'W', 'X', 'Y']
    
    # Create all possible combinations length 3
    comb_ls = list(itertools.product(acid_types,repeat=3))
    dict_cnt = dict((el,{"H": 0, "E": 0, "C": 0}) for el in comb_ls)
    
    data_df = data_df.reset_index()  # make sure indexes pair with number of rows
    for index, row in data_df.iterrows():
        seq = row['Sequence']
        tar = row['Target']
        
        for idx,elem in enumerate(seq):
            triple = ['-','','-']
            cnt_p = tar[idx]

            if idx == 0:
                triple[1] = elem
                triple[2] = seq[idx+1]

            elif idx == len(seq)-1:
                triple[0] = seq[idx-1]
                triple[1] = elem
            else:
                triple = seq[idx-1:idx+2]

            triple = tuple(triple)   
            dict_cnt[triple][cnt_p] += 1
            
    return dict_cnt


# In[113]:


import random
def creat_dec_dict(freq_dict, random_unkown = True):
    acid_types = ['-','A', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L', 'M', 'N', 'P', 'Q','R', 'S', 'T', 'V', 'W', 'X', 'Y']
    allowed_tok = ['H','E','C']
    
    # Create all possible combinations length 3
    comb_ls = list(itertools.product(acid_types,repeat=3))
    dict_dec = dict((el,'') for el in comb_ls)
    
    for elem in comb_ls:
        if freq_dict[elem]['H'] > freq_dict[elem]['E'] and freq_dict[elem]['H'] > freq_dict[elem]['C']:
            dict_dec[elem] = 'H'
        elif freq_dict[elem]['E'] > freq_dict[elem]['C']:
            dict_dec[elem] = 'E'
        elif freq_dict[elem]['C'] > freq_dict[elem]['E']:
            dict_dec[elem] = 'C'
        # Havent seen that combo... Or all equal -> Choose randomly
        else:
            if random_unkown:
                dict_dec[elem] = random.choice(allowed_tok)
            else:
                dict_dec[elem] = 'C'
                
    return dict_dec


# In[136]:


# Always predicts randomly for not known
class tripletPredict:
    
    def __init__(self, dict_p): 
        self.dict_p = dict_p
        
    def get_name(self):
        return "tripletPredict"
    
    def predict(self, sequence):
        prediction = []
        
        for idx, elem in enumerate(sequence):
            triple = ['-','','-']
            
            if idx == 0:
                triple[1] = elem
                triple[2] = sequence[idx+1]

            elif idx == len(sequence)-1:
                triple[0] = sequence[idx-1]
                triple[1] = elem
            else:
                triple = sequence[idx-1:idx+2]
                #print("elem", elem)
                #print("triple", triple)
            triple = tuple(triple)
            #print("Triple: ", triple)
            #print("Prediction: ", self.dict_p[triple])
            prediction.append(self.dict_p[triple])
            
        return prediction 


# In[137]:


def number_of_correct_pred(prediction, target_seq):
    count = 0
    for idx, tok in enumerate(prediction):
        if tok == target_seq[idx]:
            count = count + 1
    return count
    
def evalute_basic_predictor(predictor, filenames):
    for filename in filenames:
        print("-------------------------------------------------")
        print("Evaluating predictior:", predictor.get_name(), " on: ", filename, "\n")
        
        data_df = parse_sequences(filename)
        
        n_of_predicted_tokens = 0
        n_of_correct_predictions = 0
        
        data_df = data_df.reset_index()
        
        for index, row in data_df.iterrows():
            prediction = predictor.predict(row['Sequence'])
          
            n_of_predicted_tokens += len(prediction)
            n_of_correct_predictions += number_of_correct_pred(prediction, row['Target'])

        acc = n_of_correct_predictions / n_of_predicted_tokens
        print("The accuracy is:", acc, "\n")


# In[99]:


#Init Predictor
freq_dict = parse_triplet_frequency('trainSS.txt')
decision_dict = creat_dec_dict(freq_dict)


# In[138]:


# Evaluate this approach
predictor = tripletPredict(decision_dict)
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[140]:


#Predict C for unknown
freq_dict = parse_triplet_frequency('testSS.txt')
decision_dict = creat_dec_dict(freq_dict, False)


# In[141]:


predictor = tripletPredict(decision_dict)
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[147]:


#Here I evaluate the number of C's again just to be sure
def eval_c_pred(filenames):
    for filename in filenames:
        print("-------------------------------------------------")
        print("Evaluating predictior:", predictor.get_name(), " on: ", filename, "\n")
        
        data_df = parse_sequences(filename)
        
        total_length = 0#
        n_of_c = 0#
        
        n_of_predicted_tokens = 0
        n_of_correct_predictions = 0
        
        data_df = data_df.reset_index()
        
        for index, row in data_df.iterrows():
            
            total_length+= len(row['Sequence'])
            n_of_c += row['Target'].count('C')

        acc = n_of_c / total_length
        print("The accuracy for all c is:", acc, "\n")


# In[148]:


freq_dict = eval_c_pred(["trainSS.txt", "testSS.txt"])


# In[ ]:




