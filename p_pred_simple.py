#!/usr/bin/env python
# coding: utf-8

# Basic Stuff Implementations
# Author: Carl Winkler

# In[362]:


import pandas as pd


# In[363]:


def parse_counts(fileName, absolutes = False):
    # We always expect the order C, E, H
    colums_pred = ['AminoAcid','Coil_p','Strand_p','Helix_p']
    data = []
    with open(fileName) as f:
        lines = f.readlines() # list containing lines of file
        i = 0
        A_row = []

        for line in lines:
            elems = line.split('\t')
            if i == 0:
                # Amino is first column
                A_row.append(line[0])
                #Count for Coil
                A_row.append(int(elems[1]))
            elif i == 1:
                #Count for Strand
                A_row.append(int(elems[1]))
            elif i == 2:
                #Count for Helix
                A_row.append(int(elems[1]))
                data.append(A_row) 

            i = i + 1
            if i > 2:
                i = 0
                A_row = []

    # Calculate the probabilites
    if absolutes == False:
        for elem in data:
            total_Count = sum(elem[1:])
            elem[1] = elem[1] / total_Count
            elem[2] = elem[2] / total_Count
            elem[3] = elem[3] / total_Count

    # Work with pd data frame here
    return pd.DataFrame(data, columns = colums_pred)  


# In[359]:


#These are the different ACIDs we use for the last part
parse_counts("testCounts.txt",False)


# In[364]:


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


# In[366]:


def getPredDict(dataFrame):
    dict_p = {}
    dataFrame = dataFrame.reset_index()  # make sure indexes pair with number of rows
    
    for index, row in dataFrame.iterrows():
        token = ""
        if row['Coil_p'] > row['Helix_p'] and row['Coil_p'] > row['Strand_p']:
            token = "C"
        elif row['Strand_p'] > row['Helix_p']:
            token = "E"
        else:
            token = "H"
        dict_p[row['AminoAcid']] = token
        
    return dict_p


# In[198]:


df = parse_counts("trainCounts.txt")
print("The naive approach maps the tokens as follows:", getPredDict(df))


# In[365]:


def number_of_correct_pred(prediction, target_seq):
    count = 0
    for idx, tok in enumerate(prediction):
        if tok == target_seq[idx]:
            count = count + 1
    return count
    


# In[367]:


# Predicts the secondary structure which the acid is seen in mostly
class basic_naive_predict:
    
    def __init__(self, dict_p): 
        self.dict_p = dict_p
    def get_name(self):
        return "Naive-Greedy"
    
    def predict(self, sequence):
        prediction = []
        for token in sequence:
            prediction.append(self.dict_p[token])
        return prediction  


# In[223]:


# Always predicts c
class always_c_predict:
        
    def get_name(self):
        return "Always-C"
    
    def predict(self, sequence):
        prediction = []
        for token in sequence:
            prediction.append('C')
        return prediction  


# In[224]:


import random
# Predicts randomly
class random_predict:
        
    def get_name(self):
        return "Random-Token"

    def predict(self, sequence):
        allowed_tok = ["H","E","C"]
        prediction = []
        for token in sequence:
            prediction.append(random.choice(allowed_tok))
        return prediction  


# In[340]:


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
        


# In[233]:


counts_train_df = parse_counts("trainCounts.txt")
counts_test_df = parse_counts("testCounts.txt")
pred_dict_train = getPredDict(counts_train_df)
pred_dict_test = getPredDict(counts_test_df)


# In[234]:


# Evaluate naive approach based on "trainCounts.txt"
predictor = basic_naive_predict(pred_dict_train)
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[235]:


# Evaluate naive approach based on "testCounts.txt"
predictor = basic_naive_predict(pred_dict_test)
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[ ]:





# Now we evaluate the two other baselines.

# In[236]:


# Evaluate always_h_predict, It doesnt even need testCounts however its just a baseline so I didn't optimize that
predictor = always_c_predict()
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[237]:


# Evaluate random_predict, It doesnt even need testCounts however its just a baseline so I didn't optimize that as well
predictor = random_predict()
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# ## More advanced things

# In[361]:


# Always predicts c
class window_predict():
    
    def __init__(self, count_df, window_size, weights = -1): 
        self.count_df = count_df
        #self.dict_p = dict_p
        self.window_size = window_size
        countdict = {}
        
        for index, row in self.count_df.iterrows():
            #format C E H
            countdict[row['AminoAcid']] = [row['Coil_p'],row['Strand_p'],row['Helix_p']]

        self.dict_decision = countdict
        
        # Use uniform window if no window is given
        if weights == -1:
            self.weights = list = [1/window_size] * window_size 
        else:
            self.weights = weights
            
    def get_name(self):
        return "Window-Predictor"
    
    def predict(self, sequence):
        prediction = []
        for idx, token in enumerate(sequence):
            # Check if we are on the edges and predict with window_size 1 then
            
            if idx < self.window_size/2 or idx > len(sequence)-self.window_size/2: 
                row = self.count_df.loc[self.count_df['AminoAcid'] == token].iloc[0]
                token = ""
                if row['Coil_p'] > row['Helix_p'] and row['Coil_p'] > row['Strand_p']:
                    token = "C"
                elif row['Strand_p'] > row['Helix_p']:
                    token = "E"
                else:
                    token = "H"
                prediction.append(token)
                
            else:
                elem_in_window = sequence[idx - self.window_size // 2: idx + 1 + self.window_size // 2]
                wind_probs = {'C': 0, 'E': 0, 'H': 0}
                
                # Iterate over window for each element
                for wind_pos, weight in enumerate(self.weights):
                    counts = self.dict_decision[elem_in_window[wind_pos]]
                    wind_probs['C'] +=  counts[0] * weight
                    wind_probs['E'] +=  counts[1] * weight
                    wind_probs['H'] +=  counts[2] * weight
                    
                prediction.append(max(wind_probs, key=wind_probs.get))
                #print("Elements in window: ", elem_in_window)
                #print("Dict: ", wind_probs)
                #print("Decision: ", max(wind_probs, key=wind_probs.get))
                
        return prediction  


# In[341]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=3)
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[343]:


predictor = window_predict(counts_train_df, window_size=5)
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[344]:


predictor = window_predict(counts_train_df, window_size=7)
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[345]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=3,weights = [0.25,0.5,0.25])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[349]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=3,weights = [0.4,0.1,0.4])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[350]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=3,weights = [0.1,0.8,0.1])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[351]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=3,weights = [0.2,0.6,0.2])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[352]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=5,weights = [0.05,0.15,0.5,0.15,0.05])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[348]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=5,weights = [0.1,0.2,0.5,0.2,0.1])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[375]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=5,weights = [0.12,0.18,0.4,0.18,0.12])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[374]:


sum([0.12,0.18,0.4,0.18,0.12])


# In[373]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=5,weights = [0.15,0.2,0.3,0.2,0.15])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[371]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=5,weights = [0.1,0.1,0.4,0.2,0.2])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[372]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=5,weights = [0.15,0.15,0.3,0.2,0.2])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# In[369]:


# Evaluate naive approach based on "testCounts.txt"
counts_train_df = parse_counts("trainCounts.txt", True)
predictor = window_predict(counts_train_df, window_size=7,weights = [0.5,0.1,0.2,0.4,0.2,0.1,0.5])
evalute_basic_predictor(predictor, ["trainSS.txt", "testSS.txt"])


# ## Exceptional attempt

# In[ ]:


#Can be found in the other file


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


#Archive


#Calculate the accuracy from the counts directly
def acc_naive_fromCounts(fileName):
    dataFrame = parse_counts(fileName)
    dataFrame = dataFrame.reset_index()  # make sure indexes pair with number of rows
    acc = 0
    
    for index, row in dataFrame.iterrows():
        token = ""
        if row['Coil_p'] > row['Helix_p'] and row['Coil_p'] > row['Strand_p']:
            acc += row['Coil_p']
        elif row['Strand_p'] > row['Helix_p']:
            acc += row['Strand_p']
        else:
            acc += row['Helix_p']
        
    acc = acc / len(dataFrame['Coil_p'])  
    return acc

# Here we compare with the accuracies from the counts directly to see that the implementation is correct
print("ACC with trainCounts:", acc_naive_fromCounts("trainCounts.txt"))
print("ACC with testCounts:", acc_naive_fromCounts("testCounts.txt"))

