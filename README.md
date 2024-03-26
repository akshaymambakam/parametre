# The parameTRE tool
This repository is for sharing the source code of the parameTRE tool. The tool is for parametric timed pattern matching (PTPM) using parametric timed regular expressions (PTRE). We deal with three different types of PTRE. They are Parametric Signal Regular Expressions (PSRE), Parametric Timed Regular Expressions (event-based) and Parametric Event-bounded Timed Regular Expressions (PE-TRE). In addition to this we also do parametric identification for PSRE.


## Installation
We give installation instructions for Ubuntu.
1. Install Parma Polyhedra Library (PPL).
```
sudo apt install ppl-dev
```
2. Install C++ compiler.
```
sudo apt install g++
```
3. Install parsing and lexing tools.
```
sudo apt install bison
sudo apt install flex
```
4. Install Python if needed.
```
sudo apt install python
```
5. Install gnuplot
```
sudo apt install gnuplot-x11
```
6. Install make
```
sudo apt install make
```
7. Navigate into the **doctre** folder which contains the source code.
8. Once inside the folder install and generate the binary named **ptre** by running make.
```
make
```
## Reproducing results from the paper.
The paper is HSCC23_paper_1519.pdf. Instructions to reproduce the results are given in the reproduce.txt file.

## Authors and acknowledgment
This tool was developed at the VERIMAG laboratory located in the Grenoble city of France. This work is built on top of the theory of timed pattern matching developed by [Dogan Ulus](https://www.cmpe.boun.edu.tr/tr/people/dogan.ulus). Check out his [github page](https://github.com/doganulus). Timed pattern matching has been implemented in [montre](https://github.com/doganulus/montre) and [timedrel](https://github.com/doganulus/timedrel).

## License
For open source projects, say how it is licensed.
