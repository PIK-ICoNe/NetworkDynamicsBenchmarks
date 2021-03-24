#for i in 1 2 3 4 5 6 7 8 9;
for system in `sqlite3 -list -newline " " $KURABENCH_DB "select distinct id from systems;"`; 
do
    julia ./src/julia/radau_reference.jl $system 20 1e-12 1e-10
done 
