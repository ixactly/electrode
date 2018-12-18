# coding utf-8
import trajectory1 as tra_one
import trajectory2 as tra_two
import trajectory3 as tra_thr
import matplotlib.pyplot as plt
import numpy as np
itera = 30

particles1 = tra_one.calc_trajectory(itera)
particles2 = tra_two.calc_trajectory(itera)
particles3 = tra_thr.calc_trajectory(itera)



plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['font.size'] = 12 #フォントの大きさ
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ


plt.figure(figsize=(11.1 ,3))
zform1 = np.array([50, 0, 0, 2, 2, 50, 50])

yform1 = np.array([20, 20, 5, 5, 18.5, 18.5, 20])

zform2 = np.array([72, 12, 12, 14, 14, 72, 72])
yform2 = np.array([13, 13, 5, 5, 11.5, 11.5, 13])
    
zform3 = np.array([82, 102, 102, 82, 82])
zform4 = np.array([107, 137, 137, 107, 107])
zform5 = zform4 + 35

yform3 = np.array([19.7, 19.7, 21.4, 21.4, 19.7])


plt.plot(zform1, yform1, color="k")
plt.plot(zform2, yform2, color="k")
plt.plot(zform1, yform1*(-1), color="k")
plt.plot(zform2, yform2*(-1), color="k")
plt.plot(zform3, yform3, color="k")
plt.plot(zform4, yform3, color="k")
plt.plot(zform5, yform3, color="k")
plt.plot(zform3, yform3*(-1), color="k")
plt.plot(zform4, yform3*(-1), color="k")
plt.plot(zform5, yform3*(-1), color="k")


for i in range(13):
    particles1[0][i].extend(particles2[0][i])
    particles1[0][i].extend(particles3[0][i])

    particles1[1][i].extend(particles2[1][i])
    particles1[1][i].extend(particles3[1][i])

for i in range(13):
    plt.plot(particles1[0][i], particles1[1][i], color = "r")

plt.xlim(-5, 180)
plt.ylim(-25, 25)
plt.xlabel("x[mm]")
plt.ylabel("y[mm]")

plt.show()
