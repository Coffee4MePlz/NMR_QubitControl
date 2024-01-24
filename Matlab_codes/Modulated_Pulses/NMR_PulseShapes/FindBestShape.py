# this script should be run in the same path as the output files

import numpy as np
import re
import os


pwd = os.getcwd()
pattern = re.compile("Chi_2 =") # expression to be found

gates = iter(['1Xsy']) #, 'sxXsy', 'syXsy']) # defining gates in filename
with open( pwd +'/BestShapes.txt', 'w') as Best_file:
    Best_file.write("# With The Following sum limits: s_a,s_p = 6,12 \n# \n")
    for gate in gates:
        Chi_2 = np.zeros((4,36	)) # must be of the size of the output
        for NP in range(1,5): #from 50 to 200 in +50 groups
            for j in range(1,37):    
                path = pwd + "/xShape_(" + gate + ")_NP-" + str(NP*50) + "_Opt" + str(j) +"_Dec.RF"
                for i, line in enumerate(open(path)):
                    for match in re.finditer(pattern, line):
                        line_split = line.split('=', 7)
                        line_split2 = line_split[1].split(',')
                        Chi_2[NP-1][j-1] = float(line_split2[0])
            Best_Chi = np.min(Chi_2[NP-1])
            ind = np.where(Chi_2[NP-1] == Best_Chi)[0][0] +1
            Title_toprint = "# Best Chi_2 for gate " + gate +" , NP=" + str(NP*50) + ", Opt"+ str(ind) + "\n"
            content_toprint = "\n   Chi_2 = " + str(Best_Chi) + " \n\n"
            Best_file.write(Title_toprint)
            Best_file.write(content_toprint)

Best_file.close()
