SETTINGS="-n 127 -m 400 -e 1e-12"

MAP3A="-c 2.7,-2.7,-2.7,2.7 ../maps/map3.txt"
MAP3B="-c -2.7,-2.7,2.7,2.7 ../maps/map3.txt"

ALGORITHMS=(
    'GLOBAL_CHOMP'
    'TEST'
    'TNEWTON_PRECOND_RESTART'
    'MMA'
    'CCSAQ'
    'SLSQP'
    'LBFGS'
    'NEWTON'
    'VAR1'
    'VAR2'
    'TNEWTON_RESTART'
    )

GAMMAS=('0.8' '0.4' '0.2' '0.1' '0.05' '0.025'
            '0.0125' '0.005' '0.0025' '0.00125' '0.0006') 

ALPHAS=('0.8' '0.4' '0.2' '0.1' '0.05' '0.025'
            '0.0125' '0.005' '0.0025' '0.00125' '0.0006') 

CHOMP_ALGORITHMS=( 'LOCAL_CHOMP' 'GLOBAL_CHOMP' 'COVARIANT_CHOMP' )

for alg in ${ALGORITHMS[@]}; do
    
    for gamma in ${GAMMAS[@]}; do
        
        if [ $alg == 'LOCAL_CHOMP'    ] ||
           [ $alg == 'GLOBAL_CHOMP'   ] ||
           [ $alg == 'TEST' ]; then
            
            for alpha in ${ALPHAS[@]}; do
                ../build/map2d_eval -A $alg -g $gamma $SETTINGS -a $alpha  -o accel    $MAP3A
                ../build/map2d_eval -A $alg -g $gamma $SETTINGS -a $alpha  -o vel      $MAP3A
                ../build/map2d_eval -A $alg -g $gamma $SETTINGS -a $alpha  -o accel -C $MAP3A
                ../build/map2d_eval -A $alg -g $gamma $SETTINGS -a $alpha  -o vel   -C $MAP3A

                echo "../build/map2d_eval -A $alg -g $gamma $SETTINGS -a $alpha $MAP3A -o vel"
            done
        else 
            ../build/map2d_eval -A $alg -g $gamma $SETTINGS -o accel -a 0.1    $MAP3A
            ../build/map2d_eval -A $alg -g $gamma $SETTINGS -o vel   -a 0.1    $MAP3A
            ../build/map2d_eval -A $alg -g $gamma $SETTINGS -o accel -a 0.1 -C $MAP3A
            ../build/map2d_eval -A $alg -g $gamma $SETTINGS -o vel   -a 0.1 -C $MAP3A

            echo "../build/map2d_eval -A $alg -g $gamma -C $SETTINGS -a 0.1 -o vel $MAP3A "
        fi 
    done 
     
done 

