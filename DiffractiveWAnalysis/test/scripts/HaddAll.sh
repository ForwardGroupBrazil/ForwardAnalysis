#!/bin/tcsh

set JOBSOUTPUT="/storage1/dmf/CrabSubmit"
set MERGEOUTPUT="/storage1/dmf/Samples"
set DATE="10_march_2015"

mkdir $MERGEOUTPUT/DiffractiveW/$DATE
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE

mkdir $MERGEOUTPUT/DiffractiveW/$DATE/Electron2010A
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/Electron2010B
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/Muon_P1
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/Muon_P2
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/Muon_P3

mkdir $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoE_Minus
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoE_Plus
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoMu_Minus
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoMu_Plus
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/WtoENu
mkdir $MERGEOUTPUT/DiffractiveW/$DATE/WtoMuNu

mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/Electron2010A
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/Electron2010B
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/Muon_P1
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/Muon_P2
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/Muon_P3

mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoEE_Minus
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoEE_Plus
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoMuMu_Minus
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoMuMu_Plus
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/ZtoEE
mkdir $MERGEOUTPUT/DiffractiveZ/$DATE/ZtoMuMu



# Boson W Samples, muon and electron
./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWE_$DATE/Electron2010A_P1_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/Electron2010A -f ElectronA.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWE_$DATE/Electron2010B_P2_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/Electron2010B -f ElectronB.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWMu_$DATE/Muon2010_Castor_p1_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/Muon_P1 -f MuonP1.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWMu_$DATE/Muon2010_Castor_p2_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/Muon_P2 -f MuonP2.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWMu_$DATE/Muon2010_Castor_p3_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/Muon_P3 -f MuonP3.root -g 50



# Boson Z Samples, muon and electron
./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZEE_$DATE/Electron2010A_P1_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/Electron2010A -f ElectronA.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZEE_$DATE/Electron2010B_P2_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/Electron2010B -f ElectronB.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZMuMu_$DATE/Muon2010_Castor_p1_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/Muon_P1 -f MuonP1.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZMuMu_$DATE/Muon2010_Castor_p2_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/Muon_P2 -f MuonP2.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZMuMu_$DATE/Muon2010_Castor_p3_GoodCastor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/Muon_P3 -f MuonP3.root -g 50



# Boson W Monte Carlo Samples
./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWE_$DATE/SingleWtoE_minus_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoE_Minus -f SingleWtoEMinus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWE_$DATE/SingleWtoE_plus_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoE_Plus -f SingleWtoEPlus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWE_$DATE/WToE_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/WtoENu -f WtoE.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWMu_$DATE/SingleDiffractive_WMu_minus_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoMu_Minus -f SingleWtoMuMinus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWMu_$DATE/SingleDiffractive_WMu_plus_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/PomwigWtoMu_Plus -f SingleWtoMuPlus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveWMu_$DATE/WtoMu_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveW/$DATE/WtoMuNu -f WtoMu.root -g 50



# Boson Z Monte Carlo Samples
./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZEE_$DATE/SingleZtoEE_minus_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoEE_Minus -f SingleZtoEEMinus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZEE_$DATE/SingleZtoEE_plus_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoEE_Plus -f SingleZtoEEPlus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZEE_$DATE/DyToEE_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/ZtoEE -f ZtoEE.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZMuMu_$DATE/SingleDiffractive_ZMuMu_minus_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoMuMu_Minus -f SingleZtoMuMuMinus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZMuMu_$DATE/SingleDiffractive_ZMuMu_plus_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/PomwigZtoMuMu_Plus -f SingleZtoMuMuPlus.root -g 50

./mergeTTree.py -i $JOBSOUTPUT/DiffractiveZMuMu_$DATE/DyToMuMu_Castor_$DATE/res -o $MERGEOUTPUT/DiffractiveZ/$DATE/ZtoMuMu -f ZtoMuMu.root -g 50

