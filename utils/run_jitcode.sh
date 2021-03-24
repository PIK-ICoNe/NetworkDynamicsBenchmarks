echo "begin jitcode"
for system in `sqlite3 -list -newline " " $KURABENCH_DB "select distinct id from systems;"`;
do
    echo "begin --------- jitcode system number $system"
    for tol in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12;
    do
       echo "jitcode system $system, $tol"
       OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 ./src/python/jitcode_complete.py $system 20 $tol $tol
   done
    echo "complete --------- jitcode system number $system"
done
echo "complete jitcode"
