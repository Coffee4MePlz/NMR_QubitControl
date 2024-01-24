# this script should be run in the same path as the output files

import numpy as np
import re
import os

pwd = os.getcwd()
pattern = re.compile("Chi_2 =")  # expression to be found

# usual file names: Shape_(Roty(90)Xid)_PW-100NP-100_SA8SP14__Opt16_Obs.RF

gates = iter(['Roty(90)Xid'])  # , 'sxXsy', 'syXsy']) # defining gates in filename
# gates = iter(['1Xsy', 'sxXsy', 'syXsy']) # defining multiple gates in filename
with open(pwd + '/BestShapes.txt', 'w') as Best_file:
    Best_file.write("# With The Following sum limits: s_a,s_p = 8,14 \n# \n")
    for gate in gates:
        Chi_2 = np.zeros((3, 19))  # must be of the size of the output
        for PW in range(0, 3):  # from 100 to 200 in +50 groups
            for NP in range(0, 3):  # from 100 to 200 in +50 groups
                if NP <= PW:
                    for j in range(1, 20):
                        path = pwd + "/Shape_(" + gate + ")_PW-" + str(100 + 50 * PW) + "NP-" + str(100 + NP *50) + "_SA8SP14_" + "_Opt" + str(j) + "_Obs.RF"
                        for i, line in enumerate(open(path)):
                            for match in re.finditer(pattern, line):
                                line_split = line.split('=', 7)
                                line_split2 = line_split[1].split(',')
                                Chi_2[NP][j - 1] = float(line_split2[0])#Chi_2[NP - 1][j - 1] = float(line_split2[0])
                    Best_Chi = np.min(Chi_2[NP])
                    ind = np.where(Chi_2[NP] == Best_Chi)[0][0] + 1
                    Title_toprint = "# Best Chi_2 for gate " + gate + " , PW = " + str(100 + PW * 50) + " , NP=" + str(
                        100 + NP * 50) + ", Opt" + str(ind) + "\n"
                    content_toprint = "\n   Chi_2 = " + str(Best_Chi) + " \n\n"
                    Best_file.write(Title_toprint)
                    Best_file.write(content_toprint)

Best_file.close()


