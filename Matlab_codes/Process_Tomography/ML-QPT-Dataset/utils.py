
import pandas as pd
import numpy as np
import qutip as qt
import re
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
#import tensorflow as tf
#from tensorflow import *
#from sklearn.model_selection import train_test_split
#from tensorflow.python.keras import backend as K
#from tensorflow.python.keras.layers import Conv2D


########## utils functions #################

def plotQPT(matrix, labels):
    xpos, ypos = np.meshgrid(np.arange(matrix.shape[0]), np.arange(matrix.shape[1]), indexing="ij")
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)

    dx = 0.8
    dy = 0.8
    zlim_real = (-1,1)
    zlim_imag = (-0.5,0.5)

    dz_real = matrix.real.flatten()
    dz_imag = matrix.imag.flatten()

    # Color mapping 
    norm_real = Normalize(vmin=zlim_real[0]+0.25, vmax=zlim_real[1])
    colors_real = cm.viridis(norm_real(dz_real))
    norm_imag = Normalize(vmin=zlim_imag[0]+0.25, vmax=zlim_imag[1])
    colors_imag = cm.viridis(norm_imag(dz_real))
    
    fig = plt.figure(figsize=(12, 6))

    # Plot for real part
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.bar3d(xpos, ypos, zpos, dx, dy, dz_real, shade=True, color=colors_real)
    ax1.set_xticks(np.arange(matrix.shape[0]))
    ax1.set_yticks(np.arange(matrix.shape[1]))
    ax1.set_xticklabels(labels)
    ax1.set_yticklabels(labels)
    ax1.set_title('Real Part')
    ax1.set_zlim(zlim_real)

    # Plot for imaginary part
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.bar3d(xpos, ypos, zpos, dx, dy, dz_imag, shade=True, color=colors_imag)
    ax2.set_xticks(np.arange(matrix.shape[0]))
    ax2.set_yticks(np.arange(matrix.shape[1]))
    ax2.set_xticklabels(labels)
    ax2.set_yticklabels(labels)
    ax2.set_title('Imaginary Part')
    ax2.set_zlim(zlim_imag)

    plt.tight_layout()
    plt.show()

########## data set functions #################

def fetch_and_clean_dataset(path_matrix,path_labels):

    df = pd.read_csv(path_matrix, header=None)
    # Regular expression pattern for matching numbers
    pattern = re.compile(r'([-+]?\d*\.\d+|[-+]?\d+)')

    # Convert the complex numbers to pairs of real numbers
    #df = df.applymap(lambda x: [float(i) for i in re.findall(pattern, x)] + [0] if len(re.findall(pattern, x)) == 1 else [float(i) for i in re.findall(pattern, x)])

    # Convert the DataFrame to a numpy array and reshape it
    #matrices = np.array(df.values.tolist()).reshape((-1, 4, 4, 2))
    

    #uncomment below to test only the real part of the matrix.
    matrices = np.array(df.values.tolist()).reshape((-1, 4, 4, 1))
    nth_slice = 0  # Change this to select a different slice
    selected_matrix = matrices[:, :, :, 0]
    matrices = selected_matrix


    labels_df = pd.read_csv(path_labels, header=None)
    labels = labels_df.values
    return [labels, matrices]

########## NN and Tensorflow assistant functions #################

def test_matrix(y_true,m_pred_np, Ntests):
    for i in range(0,Ntests):
        rho = qt.rand_dm(4, density=1.0, dims=[[2, 2], [2, 2]])


def custom_activation(z):
    K1 = K.clip(K.tanh(z[...,:4])*np.pi, -np.pi, np.pi)
    K2 = K.clip(K.relu(z[...,4:]), 0, 0.2)
    return K.concatenate([K1,K2], axis=-1)

def numpy_loss(m_true, y_pred):
    #m_true_np = m_true.numpy()
    #y_pred_np = y_pred.numpy()
    # gets the numpy matrices, transforms them into qutip matrices, finds the fidelity 
    # between randomly Ntest state matrices and returns the loss function.
    #loss_value = np.sum(np.abs(y_true_np - y_pred_np))
    with tf.GradientTape() as tape:
        m_pred_np = tf.numpy_function(generate_choi_matrix, [y_pred], tf.float32)
        #m_pred_np = generate_choi_matrix(y_pred_np)
        loss_value = tf.reduce_sum(tf.square(m_true - m_pred_np))
        #loss_value = np.sum(np.abs(m_pred_np-m_true_np))
    gradients = tape.gradient(loss_value, model.trainable_variables)

    ''' below is the test for multiple N states as inputs'''
    #m_true = qt.Qobj(m_true_np, dims=[[2,2],[2,2]])
    #Ntests = 4
    #loss_value = test_matrix(m_true,m_pred, Ntests)
    return loss_value #tf.convert_to_tensor(loss_value, dtype=tf.float32)

def custom_loss(y_true, y_pred):
    # Use tf.numpy_function to wrap the numpy function
    return tf.numpy_function(numpy_loss, [y_true, y_pred], tf.float32)

##########  useful QIP functions #################

def Rx(thet):
    return qt.Qobj([[np.cos(thet/2), -1j*np.sin(thet/2)],
                 [-1j*np.sin(thet/2), np.cos(thet/2)]])

def Ry(thet):
    return qt.Qobj([[np.cos(thet/2), -np.sin(thet/2)],
                 [np.sin(thet/2), np.cos(thet/2)]])

def Rz(thet):
    return qt.Qobj([[np.exp(-1j*thet), 0],
                 [0, np.exp(1j*thet)]])

def generate_choi_matrix(y_pred_np):
    true_par = y_pred_np
    #inputs the set of input parameters and finds the choi matrix.
        # define the basis to run, in this case is |0>, |1>, |x+>, |y->
    input_states = [qt.basis(2, 0) * qt.basis(2, 0).dag(),
                    qt.basis(2, 1) * qt.basis(2, 1).dag(),
                    0.5 * (qt.basis(2, 0) + qt.basis(2, 1)) * (qt.basis(2, 0) + qt.basis(2, 1)).dag(),
                    0.5 * (qt.basis(2, 0) + 1j * qt.basis(2, 1))*(qt.basis(2, 0) + 1j * qt.basis(2, 1)).dag()]
    
    [thet_global , thet1 , thet2, thet3] = true_par[0:4]

    unitary = np.exp(-1j * thet_global) * Rz(thet3) * Ry(thet2) * Rx(thet1)
    output_states = [unitary * state * unitary.dag() for state in input_states]

    p = [true_par[4], true_par[5], true_par[6], true_par[7]]
    q = true_par[8]

    # Assuming you have the Apply_composite_map function defined
    for i in range(len(input_states)):
        output_states[i] = Apply_composite_map(output_states[i], p, q)


    [Choi_matrix,Kraus,ll,Soma_Kraus] = process_tomography(input_states,output_states)
    #plotQPT(Choi_matrix, ['iI', 'X', 'Y', 'Z'])
    
    m_pred = Choi_matrix
    #m_pred = qt.Qobj(Choi_matrix)

    return m_pred

def Apply_composite_map(rho_in, p, q):
    rho = rho_in
    for j in range(4):
        KrausOps = KrausOps_gen(j+1, p[j], q)  # generate the Kraus operators for j-th map
        rho = Apply_Channel(rho, KrausOps)     # apply the j-th map
    return rho

def Apply_Channel(rho, krausOperators):
    output = qt.Qobj(np.zeros(rho.shape))
    for K in krausOperators:
        output += K * rho * K.dag()
    return output


def KrausOps_gen(channel, p, q):
    id = qt.qeye(2)
    s_k = [qt.sigmax(), qt.sigmay(), qt.sigmaz()]

    M = []
    if channel in [1, 2, 3]:  # Pauli channels
        M.append(np.sqrt(1 - p/2) * id)
        M.append(np.sqrt(p/2) * s_k[channel - 1])
    elif channel == 4:  # Generalized amplitude damping
        M.append(np.sqrt(q) * qt.Qobj([[1, 0], [0, np.sqrt(1 - p)]]))
        M.append(np.sqrt(q) * qt.Qobj([[0, np.sqrt(p)], [0, 0]]))
        M.append(np.sqrt(1 - q) * qt.Qobj([[np.sqrt(1 - p), 0], [0, 1]]))
        M.append(np.sqrt(1 - q) * qt.Qobj([[0, 0], [np.sqrt(p), 0]]))
    
    return M


def process_tomography(input_states, output_states):
    op = [1j * qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()]

    A_aux = []
    B_aux = []
    for k in range(4):
        r11, r12, r21, r22 = input_states[k][0,0], input_states[k][0,1], input_states[k][1,0], input_states[k][1,1]
        r_out11, r_out12, r_out21, r_out22 = output_states[k][0,0], output_states[k][0,1], output_states[k][1,0], output_states[k][1,1]
        
        A_k = np.array([
            [r11, r12*1j, -r12,  r11*1j, -r21*1j, r22,  r22*1j,  r21, -r21, -r22*1j,  r22, -r21*1j, -r11*1j,  r12,  r12*1j,  r11],
            [r12, r11*1j,  r11, -r12*1j, -r22*1j, r21, -r21*1j, -r22, -r22, -r21*1j, -r21,  r22*1j, -r12*1j,  r11, -r11*1j, -r12],
            [r21, r22*1j, -r22,  r21*1j, -r11*1j, r12,  r12*1j,  r11,  r11,  r12*1j, -r12,  r11*1j,  r21*1j, -r22, -r22*1j, -r21],
            [r22, r21*1j,  r21, -r22*1j, -r12*1j, r11, -r11*1j, -r12,  r12,  r11*1j,  r11, -r12*1j,  r22*1j, -r21,  r21*1j,  r22]
        ])
        
        B_k = np.array([r_out11, r_out12, r_out21, r_out22])
        
        A_aux.append(A_k)
        B_aux.append(B_k)

    A_conc = np.concatenate(A_aux)
    B_conc = np.concatenate(B_aux)

    Choi_coef = np.linalg.solve(A_conc, B_conc)

    Choi_matrix = np.round(np.array([
        [Choi_coef[0], Choi_coef[1], Choi_coef[2], Choi_coef[3]],
        [Choi_coef[4], Choi_coef[5], Choi_coef[6], Choi_coef[7]],
        [Choi_coef[8], Choi_coef[9], Choi_coef[10], Choi_coef[11]],
        [Choi_coef[12], Choi_coef[13], Choi_coef[14], Choi_coef[15]]
    ]), 4)

    eigenvalues, eigenvectors = np.linalg.eig(Choi_matrix)

    Soma_Kraus = qt.Qobj(np.zeros((2, 2)))
    Kraus = []
    ll = []
    for j in range(4):
        if eigenvalues[j] >= 0.0000001:
            ll.append(eigenvalues[j])
            kraus_j = np.sqrt(eigenvalues[j]) * (eigenvectors[:, j][0] * op[0] + eigenvectors[:, j][1] * op[1] + eigenvectors[:, j][2] * op[2] + eigenvectors[:, j][3] * op[3])
            Kraus.append(kraus_j)
            Soma_Kraus += kraus_j.dag() * kraus_j

    return Choi_matrix, Kraus, ll, Soma_Kraus

