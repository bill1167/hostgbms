"""
This script is used to process input X and Y
"""
import numpy as np
import os
import string
import re 
import pandas as pd

def read_data(x_file = 'dnn_batvirus_Science_pseudoinput.csv', y_file = 'dnn_batvirus_Science_labels.csv'):
    
    print(x_file)
    print(y_file)

    X = pd.read_csv(x_file)
    X = X.values # remove heading to become array

    Y = pd.read_csv(y_file)
    Y = Y.values

    np.save('virus_data/X_data.npy', X)
    np.save('virus_data/Y_data.npy', Y)
    
    return X, Y


if __name__ == '__main__':

    X, Y = read_data()
