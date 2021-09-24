
sed -i'' -e '2s/.*/temperature				1.25/' input.dat #Temperature
./qmc1d
mv potential.dat potential_t1.dat
mv probability.dat probability_t1.dat

sed -i'' -e '2s/.*/temperature				2.5/' input.dat #Temperature
./qmc1d
mv potential.dat potential_t2.dat
mv probability.dat probability_t2.dat

sed -i'' -e '2s/.*/temperature				5/' input.dat #Temperature
./qmc1d
mv potential.dat potential_t3.dat
mv probability.dat probability_t3.dat

sed -i'' -e '2s/.*/temperature				10/' input.dat #Temperature
./qmc1d
mv potential.dat potential_t4.dat
mv probability.dat probability_t4.dat
