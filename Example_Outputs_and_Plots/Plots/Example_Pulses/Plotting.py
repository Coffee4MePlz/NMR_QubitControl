import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('Clorophorm-190us_Shape_(sxXsx)_NP-180_SA8SP14_Dec.RF', comments="#", usecols=(0, 1)) # for actual data
#data = np.loadtxt('Naive_Obs.RF', comments="#", usecols=(0, 1)) 

# extract columns into separate variables
phase = data[:, 0]
amplitude = data[:, 1]


max_diff = 180
D = np.diff(phase)

indices = np.where(abs(D) > max_diff)[0] + 1
subarrays = np.split(phase, indices)


for i in range(len(phase)):
    if (350 <= phase[i] < 360) or (0 <= phase[i] < 10):
        phase[i] = np.nan


# create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)

# plot the phase data
mask = np.abs(np.diff(phase)) > 300 # find indices where a jump of > 300 occurs
x1 = np.concatenate(([0], np.where(mask)[0]+1, [len(phase)-1]))
y1 = phase[x1]
ax1.plot(phase, color='tab:blue', lw=2.5)

last_index=0
k = 0
last_bound = 1
for index in iter(indices):
    if (D[index-1]>=0):
        sub_to_fill = subarrays[k]#[last_index:index]
        ax1.fill_between(np.arange(last_index,index),sub_to_fill , 15, interpolate=True, alpha=0.3, color='tab:blue')
        last_bound = 1
    if (D[index-1]<0):
        sub_to_fill = subarrays[k]#[last_index:index]
        ax1.fill_between(np.arange(last_index,index),sub_to_fill , 350, interpolate=True, alpha=0.3, color='tab:blue')    
        last_bound = -1
    last_index = index
    k +=1
#ax1.plot(x1, y1, 'o', color='tab:blue')

sub_to_fill = subarrays[k]#[last_index:len(phase)]
term = len(phase) - last_index
if last_bound<=0:
    ax1.fill_between(np.arange(last_index, len(phase)),sub_to_fill , 3, interpolate=True, alpha=0.3, color='tab:blue')    
else:
    ax1.fill_between(np.arange(last_index,len(phase)), sub_to_fill, 350, interpolate=True, alpha=0.3, color='tab:blue')    



# plot the amplitude data
ax2.plot(amplitude, color='tab:orange')
ax2.fill_between(np.arange(len(amplitude)), amplitude, where=amplitude>=amplitude[0], interpolate=True, alpha=0.3, color='tab:orange')

# set titles and labels
ax1.set_title('Phase Dec')
ax2.set_title('Amplitude Dec')
ax1.set_ylim([0,360])
ax2.set_ylim([0,1124])
ax2.set_xlabel('Pulse point number/timestep (us)')
ax1.set_ylabel('Degrees')
ax2.set_ylabel('Amplitude (power)')

# set the x-axis limits
ax2.set_xlim(0, len(phase)+1)

# display the plots
plt.show()
