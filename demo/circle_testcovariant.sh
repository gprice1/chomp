SETTINGS="-m -e 1e-12 -n 100 -p"

../build/circle_eval -A GLOBAL_CHOMP     $SETTINGS -o accel  
../build/circle_eval -A COVARIANT_CHOMP  $SETTINGS -o accel 
../build/circle_eval -A GLOBAL_CHOMP     $SETTINGS -o vel   
../build/circle_eval -A COVARIANT_CHOMP  $SETTINGS -o vel    

