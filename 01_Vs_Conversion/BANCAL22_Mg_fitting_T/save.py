import numpy as np
from datetime import datetime
import os 
import os.path 
import pandas as pd
import shutil

def create_output_directory(f = True):
    now = datetime.now()  #takes current time and passes it to the main program, and sets up an output directory
    now = now.strftime("%Y-%m-%d_%H-%M-%S")
    if f == True:
        os.makedirs('./output/' + now)
        shutil.copyfile('./options/data_selection.txt', './output/' + now + '/data_selection.txt')
        shutil.copyfile('./options/parameterisation_selection.txt', './output/' + now + '/parameterisation_selection.txt')
        shutil.copyfile('./options/algorithm_selection.txt', './output/' + now + '/algorithm_selection.txt')
        shutil.copytree('./data', './output/' + now + '/data', dirs_exist_ok=True)
    return now

def save_samples(samples, parameter_labels, now, n_burnin):
    df = pd.DataFrame(data = samples, columns = parameter_labels)
    df.to_csv('./output/' + now + '/samples.csv', index = False, sep = '\t')
    df2 = pd.DataFrame(data = samples, columns = parameter_labels).iloc[n_burnin:]
    df2.to_csv('./output/' + now + '/samples_postburnin.csv', index = False, sep = '\t')