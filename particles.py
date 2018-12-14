import trajectory1 as tra_one
import trajectory2 as tra_two
import trajectory3 as tra_thr
import matplotlib.pyplot as plt

itera = 2

data = tra_one(itera)
data = tra_two(itera, data)
data = tra_thr(itera, data)

for i in range(13):
    plt.plot(data[1][0][i], data[])