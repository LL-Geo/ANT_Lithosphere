import numpy as np

const_pi = np.sqrt(2*np.pi)

def prior(prior_means,prior_sigmas,y):
    prior_y = np.zeros(len(y))
    for i in range(len(y)):
        mu=prior_means[i]
        sig=prior_sigmas[i]
        dist_val=y[i]
        prior_y[i] = -0.5*(((dist_val-mu)/sig)**2) - np.log(const_pi*sig)
    return np.sum(prior_y)

def get_starting_model():
    priors = np.loadtxt('./priors.txt', skiprows = 1)
    hyperpriors = np.loadtxt('./hyperpriors.txt', skiprows = 1)

    m0 = np.random.normal(loc = priors[0,:], scale = priors[1,:]) # Initialise anelasticity model parameters randomly
    h0 = np.full(len(hyperpriors[1,:]), 0) # Initialise data hyperparameters at h=1 (log_10(h) = 0)
    x0 = np.concatenate((m0, h0)) # Join model params and hyperparams
    return x0, m0, h0, priors, hyperpriors