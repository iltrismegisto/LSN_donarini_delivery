#solid
sed -i'' -e '1s/.*/0.8/' input.dat #Temperature
sed -i'' -e '3s/.*/1.1/' input.dat #Density
sed -i'' -e '4s/.*/2.2/' input.dat #Cutoff Radius
sed -i'' -e '5s/.*/0.12/' input.dat #Delta
./Monte_Carlo_NVT.exe
#liquid
sed -i'' -e '1s/.*/1.1/' input.dat #Temperature
sed -i'' -e '3s/.*/0.8/' input.dat #Density
sed -i'' -e '4s/.*/2.5/' input.dat #Cutoff Radius
sed -i'' -e '5s/.*/0.2/' input.dat #Delta
./Monte_Carlo_NVT.exe
#gas
sed -i'' -e '1s/.*/1.2/' input.dat #Temperature
sed -i'' -e '3s/.*/0.05/' input.dat #Density
sed -i'' -e '4s/.*/5./' input.dat #Cutoff Radius
sed -i'' -e '5s/.*/3.8/' input.dat #Delta
./Monte_Carlo_NVT.exe
