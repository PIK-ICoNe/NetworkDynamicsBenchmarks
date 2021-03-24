echo "begin networkdynamics" 
for system in `sqlite3 -list -newline " " $KURABENCH_DB "select distinct id from systems;"`; 
do
    echo "begin --------- networkdynamics system number $system"
    for tol in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12; 
    do
       echo "nd.jl system  $system, $tol,..."
       JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia ./src/julia/KuraND.jl $system 20 $tol $tol
   done
    echo "complete --------- networkdynamics system number $system"
done
echo "complete networkdynamics"