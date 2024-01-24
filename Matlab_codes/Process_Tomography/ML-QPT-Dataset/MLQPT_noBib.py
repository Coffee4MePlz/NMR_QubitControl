#from numpy import array, zeros, exp, random, dot, shape, reshape, meshgrid, linspace
import numpy as np
import utils
import NN_utils
import matplotlib.pyplot as plt # for plotting
import matplotlib
matplotlib.rcParams['figure.dpi']=300 # highres display

#from NN_utils import Weights, Biases, y_layer
# for subplots within subplots:
from matplotlib import gridspec

# for nice inset colorbars: (approach changed from lecture 1 'Visualization' notebook)
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

# for updating display 
# (very simple animation)
from IPython.display import clear_output
from time import sleep


if __name__ == "__main__":
    path1 = '/home/user/nmr/Matlab_codes/Process_Tomography/ML-QPT-Dataset/QPT_UMatrices.csv'
    path2 = '/home/user/nmr/Matlab_codes/Process_Tomography/ML-QPT-Dataset/QPT_Ulabels.csv'
    [labels, matrices] = utils.fetch_and_clean_dataset(path1,path2)
    batchsize = 30
    flattened_dataset = np.stack([matrix.flatten() for matrix in matrices], axis=1)

    #num_neurons = [32, 128, 1280, 5120, 5120, 1280, 1280,  128 ,9,9]
    #num_neurons = [16, 128, 1280, 5120, 5120, 1280, 1280,  128 ,9,9] # for unitary    num_neurons = [16, 128, 1280, 5120, 5120, 1280, 1280,  128 ,9,9] # for unitary
    #activations = ['sigmoid', 'reLU','reLU','reLU', 'reLU', 'reLU', 'reLU',  'reLU',  'custom']
    num_neurons = [16, 64 ,9,9] # for unitary
    activations = ['sigmoid', 'reLU', 'custom']  
    NN_utils.init_layer_variables(num_neurons, activations, weight_scale=0.1, bias_scale=0.1, yspread=1.0)
    '''
    #single run to test
    m_true = flattened_dataset[:,0]
    y_out  = NN_utils.apply_net(m_true)
    loss = NN_utils.ChoiLoss(y_out,m_true)
    '''
    # Flatten each matrix and stack them into a single array

    first_batch = flattened_dataset[:, :batchsize].T
    first_labels_batch = labels[:batchsize]
    steps = 2900
    loss = np.zeros(steps, dtype=float)
    for j in range(steps):
        batch = flattened_dataset[:, j*batchsize:(j+1)*batchsize].T
        labels_batch = labels[j*batchsize:(j+1)*batchsize]
        if j <= 200:
            loss[j] = NN_utils.train_net(batch,labels_batch,eta=0.5) #eta=0.000001
            if j%20 == 0:
                print(f"Loss[{j}] = {loss[j]}")
                y_out = NN_utils.apply_net(flattened_dataset[:,0])
                print(f" True label = {labels[0]} \n Y_out = {y_out} \n")        
        else:
            loss[j] = NN_utils.train_net(batch,labels_batch,eta=0.1) #eta=0.000000005
            if j%100 == 0:
                print(f"Loss[{j}] = {loss[j]}")
                y_out = NN_utils.apply_net(flattened_dataset[:,0])
                print(f" True label = {labels[0]} \n Y_out = {y_out} \n")
    print(loss[-1:-5])


