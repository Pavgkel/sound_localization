import csv
import re

import matplotlib.pyplot as plt
aList=[]
plt.axes()

with open('Coordinates_3.csv', newline='') as File:
    reader = csv.reader(File,skipinitialspace=True)

    #print(*reader)
    for row in reader:
        #aList.append(row)
        #if row!=[]: plt.scatter(row[0],row[1])
        if row!=[]:
            #xy=row[0].split(',')
            #print(row[1])
            plt.scatter(float(row[0]),float(row[1]))
    plt.xlabel("Ось X")
    plt.ylabel("Ось Y")

    #plt.xlim(0, 99)
    plt.ylim(-1000, 1000)
    plt.show()