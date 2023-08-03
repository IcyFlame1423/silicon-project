import numpy as np

import matplotlib.pyplot as plt



siGFM = 0.028
siDensity = 2200
CMass = 0.2 #in grams
R = 8.31
T = 550 + 273
standardPressure = 1.01 * 10 ** 5
pressure = 0.5 * 1.01 * 10 ** 5 #pressure in the tube
poreDiameter = 10 * 10 ** -9 #initial diameter
poreVolume = 0.39
poreTolerance = 0.5 * 10 ** -9 #make sure this value is never zero
tubeLength = 0.08
tubeRadius = 0.01
flowRate = 250
silane = 0.005 #proportion of gas that is absorbed
particleRadius = 3 * 10 ** -6 # for carbon
particleVolume = (4/3) * np.pi * particleRadius ** 3 # also for carbon
SSA = 300
surfaceArea = SSA / (1/(siDensity *1000)) #surface area per volume
epsilon0 = poreVolume * 10 ** -3 / (poreVolume * 10 ** -3 + 1 / siDensity)
DTube = (0.2 * (T / 300) ** (3/2)) / (pressure / standardPressure) * 10 ** -4
D_i0 = 1 #Initial diffusivity
lambda0 = 68 * 10 ** -9 * ((T-273)/298) / (pressure / standardPressure)
velocity = flowRate/10 ** 6/60 * (standardPressure/pressure) * (T / 298) / (np.pi * tubeRadius)
pressureParticle = 0.8 * 1.01 * 10 ** 5 #pressure of the particle
c = silane * pressureParticle / R / T
k_o = 1 * 10 ** -7 /3600 * siDensity / siGFM / c / 3
k_i = 1 * 10 ** -7 /3600 * siDensity / siGFM / c / 3


outer_radii = particleRadius
inner_radii = 0
radii = outer_radii - inner_radii
radii_points = 20 #amount of points, which is equal to the number of sections which will be represented by a point
r = np.linspace(radii/2/radii_points,radii - radii/2/radii_points, radii_points) #list of points used 
dr = radii / radii_points


outer_length = tubeLength
inner_length = 0
length = outer_length - inner_length
length_points = 10
x = np.linspace(length/2/length_points, length - length/2/length_points, length_points)
dx = length / length_points

dt = 1
totalTime = 2 * 10 #2 hours
time_points = int(totalTime / dt)
criticalTime = 1800 / dt
time = np.linspace(dt, totalTime, time_points)

C = np.zeros((time_points, length_points, radii_points)) #concentration of gas
Jx = np.zeros(length_points - 1) #flux
t = np.zeros((length_points, radii_points)) #thickness
w = np.zeros((length_points, radii_points)) #pore diameter
D_i = np.ones((length_points, radii_points))  #Diffusivity
S = np.zeros((length_points, radii_points)) #surface area
inner_rou = np.zeros(length_points)
outer_rou = np.zeros(length_points)
rou = np.zeros(length_points)

P = CMass / 1000 / siDensity / (tubeRadius ** 2 * np.pi * tubeLength)  #Volume fraction of the carbon
P = np.ones(length_points) * P

y_min = 0
y_max = poreDiameter - 2 * poreTolerance #maximum thickness




# Initial Conditions
C_0 = 0.5 * 10 ** (-2) * pressure / R / T
C[0,:,radii_points-1] = C_0

T_i = np.zeros((length_points, radii_points)) #Thickness of inner
T_o = np.zeros((1, length_points)) #Thickness of surface of the particle


#phase one
for i in range(time_points):
    
    temp = w
    for a in range(length_points):
        for b in range(radii_points):
            S[a,b] = surfaceArea * (w[a,b] / poreDiameter) ** 2

    for c in range(length_points - 1):
        Jx[c] = -DTube * (C[i, c + 1, radii_points - 1] - C[i, c, radii_points - 1]) / dx + ((C[i, c + 1, radii_points - 1] + C[i, c, radii_points - 1]) / 2) * velocity
    
    for a in range(length_points):
        int_value = 0
        for b in range(radii_points):
            int_value += S[a,b] * C[i,a,b] * 4 * np.pi * r[b] ** 2 * dr

        inner_rou[a] = P[a] * k_i / particleVolume * int_value 

    for a in range(length_points):
        outer_rou[a] = 3 * P[a] * k_o * C[i, a, radii_points - 1] / particleRadius

    for a in range(length_points):
        rou[a] = inner_rou[a] + outer_rou[a]


    for l in range(length_points): #phase two
        matrixA = np.zeros((radii_points, radii_points)) #KiSc(x,r,t)/Di(r,t)
        matrixB = np.zeros((radii_points, radii_points)) #partialC * 2/r
        matrixC = np.zeros((radii_points, radii_points)) #double derivative of C
        matrixD = np.zeros((radii_points, radii_points)) #partialDi * partialC * 1/Di(r,t)
        matrixE = np.zeros((radii_points, radii_points)) #boundary matrix
        for b in range(radii_points):
            matrixA[b,b] = k_i * S[l,b] / D_i[l,b]
            matrixB[b,b] = 0
            matrixC[b,b] = -2 / (dr ** 2)
            matrixD[b,b] = 0

        for b in range(radii_points - 2):
            matrixA[b+1,b] = 0
            matrixB[b+1,b] = -1 / (r[b+1] * dr)
            matrixC[b+1,b] = 1 / (dr ** 2)
            matrixD[b+1,b] = -(D_i[l, b+2] - D_i[l,b]) / D_i[l,b+1]
        for b in range(radii_points - 2):
            matrixA[b+1,b+2] = 0
            matrixB[b+1,b+2] = 1 / (r[b+1] * dr)
            matrixC[b+1,b+2] = 1 / (dr ** 2)
            matrixD[b+1,b+2] = (D_i[l, b+2] - D_i[l,b]) / D_i[l, b+1]
        for b in range(radii_points):
            matrixA[0,b] = 0
            matrixA[radii_points - 1, b] = 0
            matrixB[0,b] = 0
            matrixB[radii_points - 1, b] = 0
            matrixC[0,b] = 0
            matrixC[radii_points - 1, b] = 0    
            matrixD[0,b] = 0
            matrixD[radii_points - 1, b] = 0




        matrixE[0,0] = 1
        matrixE[0,1] = -1
        matrixE[radii_points - 1, radii_points - 1] = 1

        b = np.zeros((radii_points, 1))
        matrixSum = matrixA + matrixB + matrixC + matrixD + matrixE
        C[i,l] = np.reshape(np.linalg.solve(matrixSum, b), (20,))


        #setting values for more matrixes
    for b in range(radii_points):
        for l in range(length_points):
            T_i[l,b] = k_i * C[i,l,b] * siGFM / siDensity
    for l in range(length_points):
        T_o[0,l] = k_o * C[i,l,radii_points - 1] * siGFM / siDensity
    for b in range(radii_points):
        for l in range(length_points):
            t[l,b] += T_i[l,b] * DTube
    for b in range(radii_points):
        for l in range(length_points):
            w[l,b] = poreDiameter - 2 * t[l,b]
    for b in range(radii_points):
        for l in range(length_points):
            porosity = epsilon0 * (w[l,b] / poreDiameter) ** 2
            D_i[l,b] = porosity ** (3/2) * D_i0 * (w[l,b] / poreDiameter) ** (1/2)
    
plt.plot(time,C[:, length_points - 1, radii_points - 1],label='')
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.show()


