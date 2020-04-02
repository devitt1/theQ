# The Q
Quantum Computing simulation project, By S.J. Devitt and C. Ferrie at QSI@UTS
# Build
On Linux or Mac simply: ```gcc main.c Simulator/sim.c Simulator/norm.c -lm -o sim```

The following worked on a Windows 10 machine:
```
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\Common7\Tools\VsDevCmd.bat"
cl /W4 main.c Simulator/norm.c Simulator/sim.c /link /out:sim.exe

```
# Current main.c file

The current main.c file creates simulations cycling between 1 and 30 qubits, and between 1 and 100 gates.  Each gate is a  Y-rotation, on a Random qubit by some random angle.  It then simply clocks the simulation time from that circuit and outputs to file the number of qubits in the simulation, the number of gates and the total time. 

For each qubit number and number of gates, a simulator is freshly initialised and subsequently destroyed before the next run.
