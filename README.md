# devittsim
State vector simulator by Simon Devitt
# Build
On Linux simply: ```gcc main.c Simulator/sim.c Simulator/norm.c -lm -o sim```

The following worked on a Windows 10 machine:
```
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\Common7\Tools\VsDevCmd.bat"
cl /W4 main.c Simulator/norm.c Simulator/sim.c /link /out:sim.exe```
