SETTINGS="-p 5 -n 127 -m 40 -e 1e-12"

MAP3A="-c 2.7,-2.7,-2.7,2.7 maps/map3.txt"
MAP3B="-c -2.7,-2.7,2.7,2.7 maps/map3.txt"

alpha1="0.00"
gamma1="0.05"

alpha2="0.1"
gamma2="0.3"

../build/map2d_eval -l GLOBAL_CHOMP -g $gamma1 $SETTINGS -a $alpha1 -o accel $MAP3A
../build/map2d_eval -k -l GLOBAL_CHOMP -g $gamma1 $SETTINGS -a $alpha1 -o accel $MAP3A 
../build/map2d_eval -l GLOBAL_CHOMP -g $gamma2 $SETTINGS -a $alpha2 -o vel $MAP3A
../build/map2d_eval -k -l GLOBAL_CHOMP -g $gamma2 $SETTINGS -a $alpha2 -o vel $MAP3A

../build/map2d_eval -l TEST -g $gamma1 $SETTINGS -a $alpha1 -o accel $MAP3A
../build/map2d_eval -k -l TEST -g $gamma1 $SETTINGS -a $alpha1 -o accel $MAP3A 
../build/map2d_eval -l TEST -g $gamma2 $SETTINGS -a $alpha2 -o vel $MAP3A
../build/map2d_eval -k -l TEST -g $gamma2 $SETTINGS -a $alpha2 -o vel $MAP3A


