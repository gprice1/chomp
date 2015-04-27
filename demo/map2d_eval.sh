SETTINGS="-p 1000 -n 127 -m 400 -e 1e-12 -b -k -d data.txt"

MAP3A="-c 2.7,-2.7,-2.7,2.7 ../demo/maps/map3.txt"
MAP3B="-c -2.7,-2.7,2.7,2.7 ../demo/maps/map3.txt"

ALGORITHMS=(
#    'CHOMP'
    'MMA'
    'CCSAQ'
    'LBFGS'
    'NEWTON'
    'TNEWTON_RESTART'
    'TNEWTON_PRECOND_RESTART'
    'VAR1'
    'VAR2'
    )

GAMMAS=('0.8' '0.4' '0.2' '0.1' '0.05' '0.025'
            '0.0125' '0.005' '0.0025' '0.00125' '0.0006') 

ALPHAS=('0.8' '0.4' '0.2' '0.1' '0.05' '0.025'
            '0.0125' '0.005' '0.0025' '0.00125' '0.0006') 


for alg in ${ALGORITHMS[@]}; do
    
    for gamma in ${GAMMAS[@]}; do
        
        if [ $alg == 'CHOMP'   ] ||
           [ $alg == 'TEST' ]; then
            
            for alpha in ${ALPHAS[@]}; do
                ../build/map2d_demo -l $alg -g $gamma $SETTINGS -a $alpha  -o accel    $MAP3A
                ../build/map2d_demo -l $alg -g $gamma $SETTINGS -a $alpha  -o vel      $MAP3A
                ../build/map2d_demo -l $alg -g $gamma $SETTINGS -a $alpha  -o accel -k $MAP3A
                ../build/map2d_demo -l $alg -g $gamma $SETTINGS -a $alpha  -o vel   -k $MAP3A

                echo "../build/map2d_eval -A $alg -g $gamma $SETTINGS -a $alpha $MAP3A -o vel"
            done
        else 
            ../build/map2d_demo -l $alg -g $gamma $SETTINGS -o accel   $MAP3A
            ../build/map2d_demo -l $alg -g $gamma $SETTINGS -o vel     $MAP3A
            ../build/map2d_demo -l $alg -g $gamma $SETTINGS -o accel  $MAP3A
            ../build/map2d_demo -l $alg -g $gamma $SETTINGS -o vel    $MAP3A

            echo "../build/map2d_eval -A $alg -g $gamma -C $SETTINGS -a 0.1 -o vel $MAP3A "
        fi 
    done 
     
done 

