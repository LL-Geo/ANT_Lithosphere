import numpy as np
from algorithm import *
from save import *
from data import *
from prior import get_starting_model
import time
import matplotlib.pyplot as plt

t_start = time.time()

now = create_output_directory(True) # get current time and produce an output directory using this timestamp
print("Beginning inversion at", now)

data, n_data = get_data() # retrieve inversion data
x0, m0, h0, priors, hyperpriors = get_starting_model() # generate a starting model

n_trials = 400000
n_burnin = int(0.5*n_trials)
n_static = 999
samples, RMS, track_posterior = run_test_algorithm(n_trials, n_burnin, n_static, x0, m0, h0, priors, hyperpriors, data, n_data)
stack = np.concatenate((track_posterior.T, samples.T, RMS.T), axis = 1)

t_end = time.time()
print("Inversion completed in", t_end - t_start, "seconds")

m_labels = np.loadtxt('./priors.txt', max_rows = 1, dtype = str) # Get param labels for saving 
h_labels = np.loadtxt('./hyperpriors.txt', max_rows = 1, dtype = str, comments = None)[1:] # Get hyperparam labels for saving
RMS_labels = []
for i in range(np.sum(n_data)):
    RMS_labels.append('RMS' + str(i + 1))
x_labels = np.concatenate((['Posterior'], m_labels, h_labels, RMS_labels)) # Join labels
save_samples(stack, x_labels, now, n_burnin) # Save model samples