#!/bin/sh

# C R E A T I N G   F O L D E R S
#################################

mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014

mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010A
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010B
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P1
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P2
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P3
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Minus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Plus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Minus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Plus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoENu
mkdir /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoMuNu
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010A
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010B
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P1
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P2
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P3
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Minus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Plus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Minus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Plus
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoEE
mkdir /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoMuMu


# D A T A
#########

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWE/Electron2010A_P1_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010A -f ElectronA.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010A/ElectronWA.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010A/ElectronA*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010A/ElectronA*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWE/Electron2010B_P2_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010B -f ElectronB.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010B/ElectronWB.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010B/ElectronB*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010B/ElectronB*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWMu/Muon2010_Castor_p1_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P1 -f MuonP1.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P1/MuonWP1.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P1/MuonP1*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P1/MuonP1*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWMu/Muon2010_Castor_p2_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P2 -f MuonP2.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P2/MuonWP2.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P2/MuonP2*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P2/MuonP2*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWMu/Muon2010_Castor_p3_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P3 -f MuonP3.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P3/MuonWP3.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P3/MuonP3*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P3/MuonP3*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZEE/Electron2010A_P1_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010A -f ElectronA.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010A/ElectronZA.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010A/ElectronA*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010A/ElectronA*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZEE/Electron2010B_P2_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010B -f ElectronB.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010B/ElectronZB.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010B/ElectronB*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Electron2010B/ElectronB*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZMuMu/Muon2010_Castor_p1_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P1 -f MuonP1.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P1/MuonZP1.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P1/MuonP1*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P1/MuonP1*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZMuMu/Muon2010_Castor_p2_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P2 -f MuonP2.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P2/MuonZP2.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P2/MuonP2*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P2/MuonP2*.root

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZMuMu/Muon2010_Castor_p3_GoodCastor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P3 -f MuonP3.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P3/MuonZP3.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P3/MuonP3*.root
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/Muon_P3/MuonP3*.root


# M C   S A M P L E S
#####################

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWE/SingleWtoE_minus_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Minus -f SingleWtoEMinus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Minus/pomwigWtoEMinus.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Minus/SingleWtoEMinus*
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Minus/SingleWtoEMinus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWE/SingleWtoE_plus_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Plus -f SingleWtoEPlus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Plus/pomwigWtoEPlus.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Plus/SingleWtoEPlus*
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Plus/SingleWtoEPlus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWE/WToE_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoENu -f WtoE.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoENu/pythiaWtoE.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoENu/WtoE*
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoENu/WtoE*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWMu/SingleDiffractive_WMu_minus_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Minus -f SingleWtoMuMinus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Minus/pomwigWtoMuMinus.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Minus/SingleWtoMuMinus*
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Minus/SingleWtoMuMinus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWMu/SingleDiffractive_WMu_plus_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Plus -f SingleWtoMuPlus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Plus/pomwigWtoMuPlus.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Plus/SingleWtoMuPlus*
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Plus/SingleWtoMuPlus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveWMu/WtoMu_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoMuNu -f WtoMu.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoMuNu/pythiaWtoMu.root /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoMuNu/WtoMu*
rm /storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoMuNu/WtoMu*



./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZEE/SingleZtoEE_minus_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Minus -f SingleZtoEEMinus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Minus/pomwigZtoEEMinus.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Minus/SingleZtoEEMinus*
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Minus/SingleZtoEEMinus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZEE/SingleZtoEE_plus_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Plus -f SingleZtoEEPlus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Plus/pomwigZtoEEPlus.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Plus/SingleZtoEEPlus*
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoEE_Plus/SingleZtoEEPlus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZEE/DyToEE_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoEE -f ZtoEE.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoEE/pythiaZtoEE.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoEE/ZtoEE*
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoEE/ZtoEE*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZMuMu/SingleDiffractive_ZMuMu_minus_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Minus -f SingleZtoMuMuMinus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Minus/pomwigZtoMuMuMinus.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Minus/SingleZtoMuMuMinus*
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Minus/SingleZtoMuMuMinus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZMuMu/SingleDiffractive_ZMuMu_plus_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Plus -f SingleZtoMuMuPlus.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Plus/pomwigZtoMuMuPlus.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Plus/SingleZtoMuMuPlus*
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/PomwigZtoMuMu_Plus/SingleZtoMuMuPlus*

./mergeTTree.py -i /storage/dmf/uerj-1/CrabSubmit/DiffractiveZMuMu/DyToMuMu_Castor_4_november_2014/res -o /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoMuMu -f ZtoMuMu.root -g 50
hadd /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoMuMu/pythiaZtoMuMu.root /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoMuMu/ZtoMuMu*
rm /storage/dmf/uerj-1/Samples/DiffractiveZ/4_november_2014/ZtoMuMu/ZtoMuMu*

