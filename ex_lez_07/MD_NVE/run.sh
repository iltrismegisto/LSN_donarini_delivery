rm argon_ave*.dat
rm ave*.dat
rm output*.*

#generic
cp config.fcc config.0
cp input.generic input.dat
./MolDyn_NVE.exe

sed  -i'' -e '1s/.*/1/' input.dat #restart

./MolDyn_NVE.exe
./MolDyn_NVE.exe

sed  -i'' -e '2s/.*/1/' input.dat #restart
./MolDyn_NVE.exe

#solid
cp config.fcc config.0
cp input.solid input.dat
./MolDyn_NVE.exe

#liquid
cp config.fcc config.0
cp input.liquid input.dat
./MolDyn_NVE.exe

sed  -i'' -e '1s/.*/1/' input.dat # restart

rm -rf *epot_liquid*
rm -rf *ekin_liquid*
rm -rf *etot_liquid*
rm -rf *temp_liquid*
rm -rf frames/*.xyz

./MolDyn_NVE.exe

#gas
cp config.fcc config.0
cp input.liquid input.dat
./MolDyn_NVE.exe
cp config.fcc config.0
cp input.gas input.dat

sed  -i'' -e '1s/.*/1/' input.dat # restart

rm -rf *epot_gas*
rm -rf *ekin_gas*
rm -rf *etot_gas*
rm -rf *temp_gas*

./MolDyn_NVE.exe
