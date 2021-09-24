#----------------Remove Files-------------------
rm best*.gene
rm first*.gene
rm evolution*.dat
#----------------Circumference Distribution----------------------
sed -i'' -e '8s/.*/0/' input.dat #City Distribution
./GA_TSP.exe
#----------------Square Distribution----------------------
sed -i'' -e '8s/.*/1/' input.dat #City Distribution
./GA_TSP.exe
