# Quplexity

Quplexity is a blazingly fast and lightweight library that provides Quantum Computer (QC) simulators with their mathematical and essential logic. 
Quplexity is written in "Assembly" aka "Assembler" for insanely fast execution/computation times. Quplexity has exstensive and carefully crafted documentation, to help anyone no matter what their technological fluency is, to integrate Quplexity into their project or contribute to the project itself. Documentation and Examples can be found in the folder (./Documentation). Quplexity is currently in the process of being intergrated into Qrack (https://github.com/unitaryfund/qrack) to provide accelerated performance and better tailored hardware support. 

## Do I need to know Assembly to use Quplexity?
Not at all, Quplexity is designed to make it easy and convinient to get started with. Quplexity has exstensive and carefully crafted documentation, to help anyone no matter what their technological fluency is, intergrate Quplexity into their project or contribute to the project itself. Documentation and Examples can be found in the folder (./Documentation).

## Getting Started!

### Dependencies

Install the following to build/run Quplexity on your machine: 
* nasm (for intel ASM)
* as   (for ARM/ARM64 ASM)
* gcc & g++

### Compiling and Linking!

After you installed the dependencies and ensured everything is working your ready to start using Quplexity:
###### Compiling:
```bash
nasm -f elf64 assembly_file.asm -o assembly_object_file.o
```
###### Then link with your C/C++ file:
```bash
gcc -no-pie cpp_file.asm assembly_object_file.o -o test
```
###### To run the example above:
```bash
./test
```

#### For Extensive and indepth documentation please see the *.pdf files found in the folder/directory (./Documentation)

## Authors

Jacob Gill  

contact: jacobygill@outlook.com 

Discord: @mrgill0651
