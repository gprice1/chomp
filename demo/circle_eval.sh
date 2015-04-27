SETTINGS="-p 10 -e 1e-12"


ALGORITHMS=(
    'TEST'
    'CHOMP'
    'MMA'
    'CCSAQ'
    'SLSQP'
    'LBFGS'
    'NEWTON'
    'TNEWTON_RESTART'
    'TNEWTON_PRECOND_RESTART'
    'VAR1'
    'VAR2'
    )

ALPHAS=('0.8' '0.4' '0.2' '0.1' '0.05' '0.025'
            '0.0125' '0.005' '0.0025') 

STARTING_N=('7' '15' '31' '63' '127')

for alg in ${ALGORITHMS[@]}; do
    
    for n in ${STARTING_N[@]}; do
        
        if [ $alg == 'TEST'    ] ||
           [ $alg == 'CHOMP'   ]; then
            
            for alpha in ${ALPHAS[@]}; do
                ../build/circle_demo -l $alg -n $n $SETTINGS -a $alpha -o accel
                ../build/circle_demo -l $alg -n $n $SETTINGS -a $alpha -o vel
                ../build/circle_demo -l $alg -n $n $SETTINGS -a $alpha -o accel -k
                ../build/circle_demo -l $alg -n $n $SETTINGS -a $alpha -o vel   -k
                echo "Algorithm: $alg n = $n, alpha = $alpha" 
            done
        else 
            ../build/circle_demo -l $alg -n $n $SETTINGS -o accel  
            ../build/circle_demo -l $alg -n $n $SETTINGS -o vel    
            ../build/circle_demo -l $alg -n $n $SETTINGS -o accel -k 
            ../build/circle_demo -l $alg -n $n $SETTINGS -o vel   -k 

            echo "Algorithm: $alg, n = $n" 
        fi 
    done 
     
done 

