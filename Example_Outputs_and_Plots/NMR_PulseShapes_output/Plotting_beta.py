import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('Shape_(1Xsy)_NP-100_Opt2_Dec.RF', comments="#", usecols=(0, 1))



# separate the data into phase, amplitude, and timesteps
phase = data[:, 0]
amplitude = data[:, 1]

# find the indices where the phase value is close to 360
discontinuous = np.logical_or(np.isclose(phase, 0, atol=1), np.isclose(phase, 360, atol=5))
phase[discontinuous] = np.nan

# create plots
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# plot phase
ax1.plot(phase, color='blue')
ax1.set_ylabel('Phase')

# plot amplitude
ax2.plot(amplitude, color='red')
ax2.set_ylabel('Amplitude')

# set x-axis label
plt.xlabel('Timesteps')

# display the plot
plt.show()
