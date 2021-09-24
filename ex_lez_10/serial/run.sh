#----------------Remove Files-------------------
rm best*.gene
rm first*.gene
rm evolution*.dat

#----------------Circumference Distribution----------------------
sed -i'' -e '5s/.*/0/' input.dat #City Distribution
sed -i'' -e '2s/.*/1/' input.dat #City Distribution
./sim_ann.exe
sed -i'' -e '2s/.*/3/' input.dat #City Distribution
./sim_ann.exe
sed -i'' -e '2s/.*/10/' input.dat #City Distribution
./sim_ann.exe
#----------------Square Distribution----------------------
sed -i'' -e '5s/.*/1/' input.dat #City Distribution
sed -i'' -e '2s/.*/1/' input.dat #City Distribution
./sim_ann.exe
sed -i'' -e '2s/.*/3/' input.dat #City Distribution
./sim_ann.exe
sed -i'' -e '2s/.*/10/' input.dat #City Distribution
./sim_ann.exe
