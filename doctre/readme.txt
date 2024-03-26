# Description of the host platform:
# RAM memory: 16 GB
# Processor: Intel® Core™ i7-8665U CPU @ 1.90GHz × 8
# Operating system: Ubuntu 20.04.5 LTS 64-bit

# Requirements are as follows:
# Operating system: Ubuntu 20.04
# Python: Generally comes with Ubuntu
# Parma Polyhedra library
sudo apt install ppl-dev
# c++ compiler
sudo apt install g++
# gnuplot: Only to produce figures 5b and 5c. Otherwise, can skip.
sudo apt install gnuplot-x11

# Producing results:
# Start in the folder "doctre" which contains everything needed

# Example to record time taken to execute commands:
python timeRec.py "./ptre eloop.ptre testBlowup.csv 3 1"
# Will run the command and prints out the time taken in seconds at the end

# To install and produce the tool binary:
# Check if on line 13 of file "ptre.cpp" in "doctre" folder is: #define FILTER (1)
# Clean the installation:
make clean
# Run "make" in "doctre" folder where all the data and source code is present
make
# The binary "ptre" is produced


# To produce the results in Table 1 for synthetic examples.
./ptre intersection.ptre testInter.csv 3 0
./ptre etre.ptre testEtre.csv 2 0
./ptre kleenePlus.ptre testKleene.csv 2 0
./ptre eloop.ptre testBlowup.csv 3 1

# Now change the value of FILTER in "ptre.cpp" in "doctre" folder from (1) to (0)
# change line 13: #define FILTER (1)
# to: #define FILTER (0)
# then: run "make" again in "doctre" folder
# This disables filtering/removing polytopes that are already included in other polytopes
# This is not strictly needed but the "Matches" for "Booleanization and Matching for ECGs" will slightly change

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

# To produce image in Figure 5b
gnuplot ecgvis1.p
# To produce image in Figure 5c
gnuplot ecgvis2.p
