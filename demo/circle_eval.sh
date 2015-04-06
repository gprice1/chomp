SETTINGS="-m 400 -e 1e-12"


ALGORITHMS=(
    'LOCAL_CHOMP'
    'GLOBAL_CHOMP'
    'COVARIANT_CHOMP'
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

ALPHAS=('0.8' '0.4' '0.2' '0.1' '0.05' '0.025'
            '0.0125' '0.005' '0.0025') 

STARTING_N=('15' '31' '63' '127')

for alg in ${ALGORITHMS[@]}; do
    
    for n in ${STARTING_N[@]}; do
        
        if [ $alg == 'LOCAL_CHOMP'    ] ||
           [ $alg == 'GLOBAL_CHOMP'   ] ||
           [ $alg == 'COVARIANT_CHOMP' ]; then
            
            for alpha in ${ALPHAS[@]}; do
                ../build/circle_eval -A $alg -n $n $SETTINGS -a $alpha -o accel
                ../build/circle_eval -A $alg -n $n $SETTINGS -a $alpha -o vel
                echo "Algorithm: $alg n = $n, alpha = $alpha" 
            done
        else 
            ../build/circle_eval -A $alg -n $n $SETTINGS -o accel  
            ../build/circle_eval -A $alg -n $n $SETTINGS -o vel    
            ../build/circle_eval -A $alg -n $n $SETTINGS -o accel -C 
            ../build/circle_eval -A $alg -n $n $SETTINGS -o vel   -C 

            echo "Algorithm: $alg, n = $n" 
        fi 
    done 
     
done 

