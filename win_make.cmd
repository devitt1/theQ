call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\Common7\Tools\VsDevCmd.bat"
cl /W4 main.c Simulator/norm.c Simulator/sim.c /link /out:sim.exe
