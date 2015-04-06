SETTINGS="-m -e 1e-12 -n 100 -p"

../build/circle_eval -A GLOBAL_CHOMP  $SETTINGS -o accel  
../build/circle_eval -A GLOBAL_CHOMP  $SETTINGS -o accel -C
../build/circle_eval -A GLOBAL_CHOMP  $SETTINGS -o vel   
../build/circle_eval -A GLOBAL_CHOMP  $SETTINGS -o vel   -C

../build/circle_eval -A TEST $SETTINGS -o accel  
../build/circle_eval -A TEST $SETTINGS -o accel -C
../build/circle_eval -A TEST $SETTINGS -o vel   
../build/circle_eval -A TEST $SETTINGS -o vel   -C

