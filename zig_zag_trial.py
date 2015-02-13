import matplotlib.pyplot as plt
import scipy.optimize
from numpy import*

time = [0, 7, 9, 12, 14, 17, 19, 27, 30, 33, 35, 37, 40, 42, 43, 45, 46, 48, 50, 67, 69, 74, 75, 76, 77, 78, 80, 83, 86, 88, 99, 101, 102, 104, 106, 108, 110]
drift = [0, 2.5, 5, 7.5, 10.0, 12.5, 15.0, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 0, -2.5, -5.0, -7.5, -10.5, -12.5, -10.0, -7.5, -5.0, -2.5, 0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 0]
delta = [0, 10, 10, -10, -10, 10, 10, -10, -10]
delta_time = [0, 3.38, 13.46, 16.93, 45.46, 49.02, 79.69, 83.38, 108.27]

Px = array(time)
Py = array(drift)

init_parameter = [4, 0.1, 0.5]
                  
def fit_func(parameter, x, y):
        A = parameter[0]
        theta = parameter[1]
        phi = parameter[2]
        residual = y - (A * sin(x * theta - phi))
        return residual

result = scipy.optimize.leastsq(fit_func, init_parameter, args = (Px, Py))
print (result[0])
        
def theta_curve(x, p):
        return (p[0] * sin(x * p[1] - p[2]))

x = arange(0, 110, 1)
graph = theta_curve(x, result[0])
fig = plt.figure()

plt.ylabel('delta:rudder angle(degree) or theta:rate of turn(degree)')
plt.xlabel('time(sec)')

plt.title("Zig Zag trial(veering: 10degrees)")
plt.axhline(y=0, color="black") 
plt.grid()
plt.plot(x, graph, label = "theta", linewidth = 1)
plt.plot(delta_time, delta, "-", label = "delta", linewidth = 3)
plt.plot(time, drift, "o", label = "measured")
plt.legend()

ax = fig.gca()
ax.set_xticks(arange(0, 111, 5))
ax.set_yticks(arange(-20, 20, 1))

plt.show()


