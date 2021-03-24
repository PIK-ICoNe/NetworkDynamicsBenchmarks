import sys
import numpy as np
from numpy.random import uniform
from numpy.random import choice
import sqlite3
import networkx as nx
import os

db=os.environ['KURABENCH_DB']
conn = sqlite3.connect(db)

np.random.seed(112020)

c = 3.0
q = 0.2

system = 0 # For system names
for n in [10,100,1000]:
    edges = q*n*(n-1)
    coupleConstant = c / (n-1)
    node_degree = 4
    for k in [1,2,3,4,5,6,7,8,9,10]: # make ten of each size
        system = system + 1 # name of system
        # Generate edges first so we can populate the systems table
        # with the correct number
        graph = nx.watts_strogatz_graph(n, node_degree, q)
        A = nx.adjacency_matrix(graph)
        A = A.todense()
        sources,dests = np.where(A)
        #### System Table
        # name | # nodes | # edges | q | c | coupleConstant
        command = f'insert into systems\
                    (id,name,nodes,edges,q,c,coupling_constant) \
                    values \
                ({system},"{n}_{k}",{n},{len(sources)},{q},{c},{coupleConstant});'
        conn.execute(command)

        #### Edges Table
        for edge in range(0,len(sources)):
            command = f'insert into edges (system,id,source,dest)\
                        values \
                        ({system},{edge+1},{sources[edge]+1},{dests[edge]+1});'
            conn.execute(command)

        #### Nodes Table
        omega = uniform(-0.5,0.5,n)
        x0 = uniform(0,2*np.pi,n)
        omega.sort()
        for node in range(0,n):
            command = f'insert into nodes (system,id,omega,initialCondition)\
                        values\
                        ({system},{node+1},{omega[node]},{x0[node]});'
            conn.execute(command)

        conn.commit()
        print(f'Generated system {system}')

print("Done!")
conn.close()
