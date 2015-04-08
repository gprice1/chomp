SETTINGS="-m -e 1e-12 -n 100 -p"

../build/circle_eval -l GLOBAL_CHOMP  $SETTINGS -o accel  
../build/circle_eval -l GLOBAL_CHOMP  $SETTINGS -o accel -k
../build/circle_eval -l GLOBAL_CHOMP  $SETTINGS -o vel   
../build/circle_eval -l GLOBAL_CHOMP  $SETTINGS -o vel   -k

../build/circle_eval -l TEST $SETTINGS -o accel  
../build/circle_eval -l TEST $SETTINGS -o accel -k
../build/circle_eval -l TEST $SETTINGS -o vel   
../build/circle_eval -l TEST $SETTINGS -o vel   -k

