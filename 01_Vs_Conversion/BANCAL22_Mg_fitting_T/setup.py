import numpy as np
import shutil

def select_data():

    data_type = ['plate', 'adiabat', 'attenuation', 'viscosity']
    data_name = ['plate.VseTz', 'adiabat.VseTz', 'attenuation.QeVsz', 'viscosity.neVsz']
    data_plate = []
    data_adiabat = []
    data_attenuation = []
    data_viscosity = []
    data_selection = np.loadtxt('./options/data_selection.txt', skiprows = 1)
    if data_selection[0] == 1:
        data_plate = [np.loadtxt('./data/' + data_type[0] + '/' + data_name[0], skiprows = 1).T]
    if data_selection[1] == 1:
        data_adiabat = [np.loadtxt('./data/' + data_type[1] + '/' + data_name[1], skiprows = 1).T]
    if data_selection[2] == 1:
        data_attenuation = [np.loadtxt('./data/' + data_type[2] + '/' + data_name[2], skiprows = 1).T]
    if data_selection[3] == 1:
        data_viscosity = [np.loadtxt('./data/' + data_type[3] + '/' + data_name[3], skiprows = 1).T]
    data_length = [len(data_plate), len(data_adiabat), len(data_attenuation), 0]
    n_data = int(np.sum(np.asarray(data_selection)*np.asarray(data_length)))

    hyperpriors = np.array([np.zeros(n_data), np.ones(n_data)])
    hyperprior_header = ''
    for i in range(n_data):
        hyperprior_header += 'h' + str(i+1) + ' '
    np.savetxt('./hyperpriors.txt', hyperpriors, header = hyperprior_header, fmt = '%.0f')

    return None

def select_param():

    param_selection = str(np.genfromtxt('./options/parameterisation_selection.txt', dtype = 'str'))
    shutil.copyfile('./anelasticity_parameterisation/'+param_selection+'/thermodynamic_conversions.py', './thermodynamic_conversions.py')
    shutil.copyfile('./anelasticity_parameterisation/'+param_selection+'/priors.txt', './priors.txt')

    return None

def select_algorithm():

    algorithm_selection = str(np.genfromtxt('./options/algorithm_selection.txt', dtype = 'str'))
    shutil.copyfile('./algorithm/'+algorithm_selection+'.py', './algorithm.py')

    return None

select_data()
select_param()
select_algorithm()