import numpy as np

tempo = 0
change = False

pos = [[0, 0, 0],[0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1]]

pos2 = []

for i in range(12):
    temp = []
    for j in range(12):
        temp = np.append(temp, np.random.randint(0, 2))
    pos2 = np.append(pos2, [temp])

print(pos2)

razao = 2

def sum():
    s = 0
    for i in range(5):
        for j in range(np.size(pos[i])):
            s += pos[i][j]
    return s 

print(sum())

print(pos)

while((pos[0][0]*pos[0][1]*pos[0][2] == 0) and (sum() != 0)) :
    for i in range(5) :
        for j in range(np.size(pos[i])) : 
            change = False
            if((i == 0) and (pos[i][j] == 1)) :
                pos[i][j] = 0
                change = True
            if(((change == False) and ((j == 0) and (i != 0)) and ((pos[i][j] == 1) and (pos[i][j+1] == 0)))) :
                pos[i][j] = 0
                pos[i][j+1] = 1
                change = True
            if(((change == False) and ((j == (np.size(pos[i]) - 1)) and (i != 0)) and (pos[i][j] == 1))) :
                if((pos[i-1][j-2] == 0)):
                    pos[i][j] = 0
                    pos[i-1][j-2] = 1
                elif((pos[i][j-1] == 0)) :    
                    pos[i][j] = 0
                    pos[i][j-1] = 1
                change = True
            if((i > 0) and ((change == False) and ((j < (np.size(pos[i]) - 1)) and j > 0))):
                if((pos[i-1][j-1] == 0) and (pos[i][j] == 1)) :
                    pos[i][j] = 0
                    pos[i-1][j-1] = 1
    tempo += 1
    print(pos)

print(tempo)