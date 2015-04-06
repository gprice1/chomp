SETTINGS="-n 127 -m 40 -e 1e-12"

MAP3A="-c 2.7,-2.7,-2.7,2.7 ../maps/map3.txt"
MAP3B="-c -2.7,-2.7,2.7,2.7 ../maps/map3.txt"

gamma1="0.5"
alpha1="0.02"
gamma2="0.2"
alpha2="0.03"

../build/map2d_eval -A GLOBAL_CHOMP    -g $gamma1 $SETTINGS -a $alpha1 $MAP3B -o accel
../build/map2d_eval -A COVARIANT_CHOMP -g $gamma1 $SETTINGS -a $alpha1 $MAP3B -o accel
../build/map2d_eval -A GLOBAL_CHOMP    -g $gamma2 $SETTINGS -a $alpha2 $MAP3B -o vel 
../build/map2d_eval -A COVARIANT_CHOMP -g $gamma2 $SETTINGS -a $alpha2 $MAP3B -o vel



