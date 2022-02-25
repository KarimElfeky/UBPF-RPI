# Importing COM Interface and Creating the OpenDSS engine object.
import subprocess
from timeit import default_timer, timeit
#import win32com.client
import dss
from guppy import hpy

#h = hpy()
#print(h.heap())
#dssObj = win32com.client.Dispatch("OpenDSSEngine.DSS")
dssObj = dss.DSS
#The Text is used to enter command in the DSS script
dssText = dssObj.Text

#DSS circuit is the class of the active circuit in the DSS script 
dssCircuit = dssObj.ActiveCircuit

#DSS Solution is the solution property in the circuit class of the DSS.
dssSolution = dssCircuit.Solution

#DSS Elem is the active element class, this can be used to modify any element. Howver
#it is mostly useful when modifying elements that do not have COM Interface like PV.
dssElem = dssCircuit.ActiveCktElement

#Load Interface
dssLoad = dssCircuit.Loads
dssLoadshapes = dssCircuit.LoadShapes
# Setting the file path and the dss file path.
dssText.Command = "set datapath=/home/pi/Documents/Karim/Unbalanced_PowerFlow-main/OpenDSS_Files-main"
dssText.Command = "Redirect DGrid_OpenDSS_case1_double_v2.dss"
start = timeit,default_timer()
dssSolution.Solve()
stop = timeit,default_timer()
print('Total Time: ',stop[1]-start[1])
subprocess.call("/home/pi/Documents/Karim/Unbalanced_PowerFlow-main/memmon_v1")
dssText.Command = "set datapath=/home/pi/Documents/Karim/Unbalanced_PowerFlow-main/OpenDSS_Files-main/Results"
dssText.Command = "export Voltages"
dssText.Command = "export Currents"
