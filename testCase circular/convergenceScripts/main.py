import readFile
import matplotlib.pyplot as plt
import numpy as np
import os, shutil

residualForFields=["Ux_0","Uy_0","Uz_0","k_0","epsilon_0","M3_0"]
labels=[r"$U_{rel,x}$",r"$U_{rel,y}$",r"$U_{rel,z}$",r"$\kappa$",r"$\epsilon$",r"$m_3$"]
os.chdir("../logs")

for i in range(0,np.size(residualForFields)):
    x=readFile.readFile(residualForFields[i]).x()
    y=readFile.readFile(residualForFields[i]).y()
    plt.plot(x,y,label=labels[i])
    if i==0:
        print("\n"+"|| Residuals check ||")
        print("Final time is: "+str(x[-1]))
        print("Final residuals are:")
    print(str(residualForFields[i])+": "+str(y[-1]))
os.chdir("../convergenceScripts")

# plotting OH- concentration
plt.xlabel(r'$time,\ s$', fontsize=20)
plt.title(r"$residuals$")
plt.legend(fontsize=15,handlelength=3)
plt.ylim([1e-10,1])
plt.yscale("log")
plt.tight_layout()
plt.show()
# import convergenceCheck
