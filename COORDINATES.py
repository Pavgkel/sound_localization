import math
def cordinate(x1,x2,x3,x4,y1,y2,y3,y4,angle1,angle2,angle3,angle4):
    X=x1*math.tan(angle1)-x2*math.tan(angle2)-x3*math.tan(angle3)-x4*math.tan(angle4)-y1+y2-y3+y4
    Y=y1+(X-x1)*math.tan(angle1)
    return X,Y

angle = [[0, 0, 0, 0],
       [90, 90, 90, 0],
       [0, 0, 0, 0],
       [0, 90, 0, 0],
       [3.00000000e+01, -2.75335829e-31, 5.40591856e+00, -5.40591856e+00],
       [90, 0, 90, 0],
       [90, 90, 0, 90],
       [90, 90, 90, 0],
       [0, 0, 24.91500561, -24.91500561],
       [90, 0, 0, 0],
       [0, 90, 90, 90],
       [-9.00000000e+01, -8.12265493e-31, 1.27830426e+01, -1.27830426e+01],
       [0, 0, 0, 90],
       [90, 90, 0, 90],
       [0, 0, 90, 0],
       [90, 0, 0, 90],
       [90, 90, -34, 43935254, 34, 43935254],
       [90, 0, 0, 90],
       [1.28659837e-16, -1.45437344e-30, 3.96704309e+01, -3.96704309e+01],
       [30.00000053, -21.18664631, - 30.00000023, 0],
       [-8.75132114e-16, -9.00000000e+01, 6.57799796e-16, -7.65398297e-15],
       [90, 90, 0, 90],
       [0, 0, 90, 0],
       [30, 90, 10.51660265, - 10.51660265],
       [2.23261722e-15, 4.94415161e+00, -3.02390134e+01, -4.75205666e+01],
       [0, 90, 23.56696753, -23.56696753],
       [0, 0, 0, 0],
       [-2.64813263e-31, 2.13404798e+01, - 2.13404798e+01, - 8.59059559e-16],
       [0, 90, 90, 0],
       [0, 90, 0, 0],
       [0, 0, 0, 90],
       [0, 90, 0, 90],
       [90, 0, 90, 90],
       [90, -90, 62.54422823, -62.54422823],
       [6.79731171e-16, 9.00000000e+01, -1.72721838e-16, -9.00000000e+01],
       [90, 90, 0, 90],
       [90, 0, 0, 0],
       [0.00000000e+00, 9.00000000e+01, - 9.00000000e+01, 1.29555572e-15],
       [0, 0, 13.86970439, -13.86970439],
       [0, 0, 90, 0],
       [0, 90, 90, 90],
       [9.00000000e+01, -3.32956284e-38, 2.60431719e-25, 9.00000000e+01]]
x1=0
x2=50
x3=0
x4=0
y1=0
y2=50
y3=50
y4=0
import matplotlib.pyplot as plt
aList=[]
plt.axes()
for i in range(len(angle)):
    for j in range(4):
        if angle[i][j]<0:
            angle[i][j]=90+angle[i][j]
#print(angle)
#coordinates=cordinate(x1,x2,x3,x4,y1,y2,y3,y4,angle[9][0],angle[9][1],angle[9][2],angle[9][3])
xy=[]
for i in range(len(angle)):
    xy.append(cordinate(x1,x2,x3,x4,y1,y2,y3,y4,angle[i][0],angle[i][1],angle[i][2],angle[i][3]))
plt.scatter(xy[0],xy[1])
plt.show()