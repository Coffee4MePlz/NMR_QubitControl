import pandas as pd
import numpy as np
import re
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Flatten, Input, Reshape
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Activation, Dropout
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras import backend as K
from tensorflow.keras.layers import Conv2D





def fetch_and_clean_dataset(path_matrix,path_labels):

    df = pd.read_csv(path_matrix, header=None)
    # Regular expression pattern for matching numbers
    pattern = re.compile(r'([-+]?\d*\.\d+|[-+]?\d+)')

    # Convert the complex numbers to pairs of real numbers
    df = df.applymap(lambda x: [float(i) for i in re.findall(pattern, x)] + [0] if len(re.findall(pattern, x)) == 1 else [float(i) for i in re.findall(pattern, x)])

    #matrices = np.array(df.values.tolist()).reshape((-1, 4, 4, 1))
    # Convert the DataFrame to a numpy array and reshape it
    matrices = np.array(df.values.tolist()).reshape((-1, 4, 4, 2))
    #nth_slice = 0  # Change this to select a different slice

    #    selected_matrix = matrices[:, :, :, nth_slice]
    #    matrices = selected_matrix
    
    labels_df = pd.read_csv(path_labels, header=None)
    labels = labels_df.values
    return [labels, matrices]

def custom_activation(z):
    K1 = K.clip(K.tanh(z[...,:4])*np.pi, -np.pi, np.pi)
    K2 = K.clip(K.relu(z[...,4:]), 0, 0.2)
    return K.concatenate([K1,K2], axis=-1)

if __name__ == "__main__":
    #load .csv files
    path1 = './QPT_Matrices.csv'
    path2 = './QPT_labels.csv'
    [labels, matrices] = fetch_and_clean_dataset(path1,path2)
    print(f" size of labels {labels.shape}")
    # Split the data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(matrices, labels, test_size=0.2, random_state=42)

    # Build the neural network
    model = Sequential([
        Flatten(input_shape=(4, 4, 2)),
        Dense(64, activation='relu'),
        Dense(320, activation='sigmoid'),
        Dropout(.2),
        Dense(1028, activation='sigmoid'),
        Dropout(.4),
        Dense(64, activation='relu'),                
        Dropout(.2),
        Dense(32, activation='tanh'),
        Dense(9, Activation(custom_activation), name='found_par')  # 4 or 9 outputs without activation function for regression problem
    ])
    '''        
    Dense(32, activation='tanh'),
    Dense(64, activation='relu'),
    Dense(64, activation='relu'),
    Dense(16, activation='relu'),        
    Reshape((4,4,1)) 
    '''

#    model.add(Dense(9, input_shape=(9,)))
#    model.add(Activation(custom_activation))

    # Compile the model
    model.compile(optimizer='adam', loss='mse', metrics=['mae'])

    # Train the model
#    model.fit(X_train, y_train, epochs=100, batch_size=128, validation_data=(X_test, y_test))
    model.fit(X_train, y_train, epochs=200, batch_size=16, validation_data=(X_test, y_test))
    print('\n training done')
    # Generate predictions for the first 5 elements of the test set
    predictions = model.predict(X_test[:2])

    # Print the predictions
    for prediction in predictions:
        print(f"Predicted values for test set element : {prediction}")

    # You may also want to print the true labels for comparison
    true_labels = y_test[:2]
    for label in true_labels:
        print(f"True label for test set element : {label}")
        
    '''    
    intermediate_layer_model = Model(inputs=model.input,
                                 outputs=model.get_layer('found_par').output)
    input_data = X_test[:2]
    intermediate_output = intermediate_layer_model.predict(input_data)
    print(f"true labels = \n {y_test[:2]}")
    print(f"predicted labels =\n  {intermediate_output}")
    '''    


