# Description of the host platform:
# RAM memory: 16 GB
# Processor: Intel® Core™ i7-8665U CPU @ 1.90GHz × 8
# Operating system: Ubuntu 20.04.5 LTS 64-bit

# Requirements are as follows:
# Operating system: Ubuntu 20.04
# Parma Polyhedra library
sudo apt install ppl-dev
# c++ compiler
sudo apt install g++
# gnuplot: Only to produce figures 5b and 5c. Otherwise, can skip.
sudo apt install gnuplot-x11

# To install and produce the tool binary:
# Run "make" in "doctre" folder where all the data and source code is present
make
# The binary "ptre" is produced

# Producing results:

# To produce the results in Table 1 for synthetic examples.
./ptre intersection.ptre testInter.csv 3 0
./ptre etre.ptre testEtre.csv 2 0
./ptre kleenePlus.ptre testKleene.csv 2 0
./ptre eloop.ptre testBlowup.csv 3 1

# To produce results for ECGs in Table 2.
./ptre qrs.ptre 205L.csv 3 0
./ptre qrs.ptre 221L.csv 3 0
./ptre qrs.ptre 123L.csv 3 0

# ECG parametric identification: Subsubsection 6.2.2
./ptre param5ecg.ptre debugTest.csv 5 0 ecg2.label
./ptre param6ecg.ptre debugTest.csv 6 0 ecg2.label
./ptre param4ecg.ptre debugTest.csv 4 0 ecg2.label

# STL+PTRE: Booleanization and Matching for ECGs: Subsubsection 6.2.3
./ptre ecgstl.ptre bflist205.txt 2 2 # This is on the whole ECG 205 signal
./ptre ecgstl.ptre bflist221.txt 2 2 # This is on the whole ECG 221 signal
./ptre ecgstl.ptre bflist123.txt 2 2 # This is on the whole ECG 123 signal

# Marine traffic rule violation: Subsubsection 6.2.4
cp ptre ./tsdata/. # Copy the tool binary into the folder (tsdata) containing data
cd tsdata # Go into the folder containing the data
python monall.py # Run the tool on the data

# Example to record time taken to execute command:
python timeRec.py "./ptre eloop.ptre testBlowup.csv 3 1"

# To produce image in Figure 5b
gnuplot ecgvis1.p

# To produce image in Figure 5c
gnuplot ecgvis2.p