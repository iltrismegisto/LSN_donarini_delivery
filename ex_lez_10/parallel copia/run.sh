#----------------Remove Files-------------------
rm best*
rm first*
rm evolution*

#----------------Square Distribution----------------------
sed -i'' -e '5s/.*/1/' input.dat #City Distribution
sed -i'' -e '2s/.*/1/' input.dat #temp
mpirun -np 4 ./sim_ann.exe
sed -i'' -e '2s/.*/3/' input.dat #temp
mpirun -np 4 ./sim_ann.exe
sed -i'' -e '2s/.*/10/' input.dat #temp
mpirun -np 4 ./sim_ann.exe
