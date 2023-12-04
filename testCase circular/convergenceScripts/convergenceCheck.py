import readFile
import matplotlib.pyplot as plt
import numpy as np
import os, shutil

residualForFields=["volFieldValue.dat"]
labels=[r"$U_{max}$",r"$\kappa_{max}$",r"$\epsilon_{max}$",r"$M_{3}$"]

os.chdir("../postProcessing/cellMax/0")

x=readFile.readVelocity(residualForFields[0]).x()
y=readFile.readVelocity(residualForFields[0]).y()
k=readFile.readConvergence(residualForFields[0]).k()
epsilon=readFile.readConvergence(residualForFields[0]).epsilon()
M3=readFile.readConvergence(residualForFields[0]).M3()

print("\n"+"|| Convergence check ||")
print("Final velocity value (m/s): "+str(y[-1]))
print("Final k value (m/s): "+str(k[-1]))
print("Final epsilon value (m/s): "+str(epsilon[-1]))

print("Final velocity variation (m/s): "+str(abs(y[-1]-y[-2])))
print("Final k variation (m^2/s^2): "+str(abs(k[-1]-k[-2])))
print("Final epsilon variation (m^2/s^3): "+str(abs(epsilon[-1]-epsilon[-2])))
print("Final M3 variation (-): "+str(abs(M3[-1]-M3[-2])))

plt.figure()
plt.plot(x,y,label=labels[0])
plt.plot(x,k,label=labels[1])
plt.plot(x,epsilon,label=labels[2])
plt.plot(x,M3,label=labels[3])

# plotting OH- concentration
plt.xlabel(r'$time,\ s$', fontsize=20)
plt.title(r"$Convergence\ check$")
plt.legend(fontsize=15,handlelength=3)
plt.yscale("log")
plt.tight_layout()
plt.show()