import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


df = pd.read_csv (r'simulation_result.csv')
states = df.to_numpy()

time = states[:,0]
u = states[:,1]
v = states[:,2]
w = states[:,3]
p = states[:,4]
q = states[:,5]
r = states[:,6]
e0 = states[:,7]
e1 = states[:,8]
e2 = states[:,9]
e3 = states[:,10]
Xe = states[:,11]
Ye = states[:,12]
h = states[:,13]
throttle = states[:,14]
elevator = states[:,15]
aileron = states[:,16]
rudder = states[:,17]
Fx = states[:,18]
Fy = states[:,19]
Fz = states[:,20]
Mx = states[:,21]
My = states[:,22]
Mz = states[:,23]


plt.figure("Aircraft states", figsize=(16,9))
plt.subplot(4,3,1)
plt.plot(time, u)
plt.xlabel("Time (sec)")
plt.ylabel("u (m/s)")
plt.grid("true")
plt.subplot(4,3,2)
plt.plot(time, v)
plt.xlabel("Time (sec)")
plt.ylabel("v (m/s)")
plt.grid("true")
plt.subplot(4,3,3)
plt.plot(time, w)
plt.xlabel("Time (sec)")
plt.ylabel("w (m/s)")
plt.grid("true")
plt.subplot(4,3,4)
plt.plot(time, p)
plt.xlabel("Time (sec)")
plt.ylabel("p (rad/s)")
plt.grid("true")
plt.subplot(4,3,5)
plt.plot(time, q)
plt.xlabel("Time (sec)")
plt.ylabel("q (rad/s)")
plt.grid("true")
plt.subplot(4,3,6)
plt.plot(time, r)
plt.xlabel("Time (sec)")
plt.ylabel("r (rad/s)")
plt.grid("true")
plt.subplot(4,3,7)
plt.plot(time, e0)
plt.xlabel("Time (sec)")
plt.ylabel("e0")
plt.grid("true")
plt.subplot(4,3,8)
plt.plot(time, e1)
plt.xlabel("Time (sec)")
plt.ylabel("e1")
plt.grid("true")
plt.subplot(4,3,9)
plt.plot(time, e2)
plt.xlabel("Time (sec)")
plt.ylabel("e2")
plt.grid("true")
plt.subplot(4,3,10)
plt.plot(time, e3)
plt.xlabel("Time (sec)")
plt.ylabel("e3")
plt.grid("true")
plt.subplot(4,3,11)
plt.plot(time, Xe, time, Ye)
plt.xlabel("Time (sec)")
plt.ylabel("Xe and Ye (m)")
plt.legend(["Xe", "Ye"])
plt.grid("true")
plt.subplot(4,3,12)
plt.plot(time, h)
plt.xlabel("Time (sec)")
plt.ylabel("h (m)")
plt.ylim([0,np.max(h)+10])
plt.grid("true")
plt.show()

plt.figure("Control Inputs", figsize=(16,9))
plt.subplot(4,1,1)
plt.plot(time, throttle)
plt.xlabel("Time (sec)")
plt.ylabel("throttle")
plt.grid("true")
plt.subplot(4,1,2)
plt.plot(time, elevator)
plt.xlabel("Time (sec)")
plt.ylabel("elevator")
plt.grid("true")
plt.subplot(4,1,3)
plt.plot(time, aileron)
plt.xlabel("Time (sec)")
plt.ylabel("aileron")
plt.grid("true")
plt.subplot(4,1,4)
plt.plot(time, rudder)
plt.xlabel("Time (sec)")
plt.ylabel("rudder")
plt.grid("true")
plt.show()

plt.figure("Forces & Moments", figsize=(16,9))
plt.subplot(2,3,1)
plt.plot(time, Fx)
plt.xlabel("Time (sec)")
plt.ylabel("Fx")
plt.grid("true")
plt.subplot(2,3,2)
plt.plot(time, Fy)
plt.xlabel("Time (sec)")
plt.ylabel("Fy")
plt.grid("true")
plt.subplot(2,3,3)
plt.plot(time, Fz)
plt.xlabel("Time (sec)")
plt.ylabel("Fz")
plt.grid("true")
plt.subplot(2,3,4)
plt.plot(time, Mx)
plt.xlabel("Time (sec)")
plt.ylabel("Mx")
plt.grid("true")
plt.subplot(2,3,5)
plt.plot(time, My)
plt.xlabel("Time (sec)")
plt.ylabel("My")
plt.grid("true")
plt.subplot(2,3,6)
plt.plot(time, Mz)
plt.xlabel("Time (sec)")
plt.ylabel("Mz")
plt.grid("true")
plt.show()

fig = plt.figure('3D Flight Path', figsize=(16,9))
ax = fig.add_subplot(111, projection='3d')
ax.plot(Xe,Ye,h)
ax.set_xlabel('North, m')
ax.set_ylabel('East, m')
ax.set_zlabel('Altitude, m')
ax.set_title('3D Flight Path')
ax.set_xlim(min(min(Xe), min(Ye))-10, max(max(Xe), max(Ye))+10)
ax.set_ylim(min(min(Xe), min(Ye))-10, max(max(Xe), max(Ye))+10)
ax.set_zlim(0, max(h)+10)
plt.show()
