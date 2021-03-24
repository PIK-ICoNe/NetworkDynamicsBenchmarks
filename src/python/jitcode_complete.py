# example explained in Ansmann paper
import time
startTstartup= time.process_time()

import sys
import numpy as np
from scipy.sparse import csr_matrix 
from symengine import sin
from jitcode import y
import os
endTstartup= time.process_time()

# Parse command line args
systemid = int(sys.argv[1])
tend = float(sys.argv[2])
atol = float(sys.argv[3])
rtol = float(sys.argv[4])
integrator = "dopri5"
# Time loading data
startTload = time.process_time()
import sqlite3
db=os.environ['KURABENCH_DB']

systemid = sys.argv[1]

# Load data

## First read the system info, inc. # of edges for allocating arrays
conn = sqlite3.connect(db)

cursor = conn.execute(f'SELECT name, nodes, edges, coupling_constant \
                      from systems where id={systemid}')

name,nodes,edges,couplingConstant = cursor.fetchone()

## Now read the connectivity
sources = np.empty(edges,dtype=int)
dests   = np.empty(edges,dtype=int)
data    = np.ones(edges,dtype=bool)


cursor = conn.execute(f'SELECT id, source, dest\
                      from edges where system={systemid}')

for edgeId,s,d in cursor: #TODO: validate rows returned == number of edges
    sources[edgeId-1] = s-1 # NOTE: Python is zero-based indexing
    dests[edgeId-1] = d-1   # but the database stores in 1-based

## Store the node omegas 
cursor = conn.execute(f'SELECT id, omega,initialCondition\
                      from nodes where system={systemid}')

omegas = np.empty(nodes,dtype=np.float64) # TODO: Confirm number of rows == # nodes expected
x0 = np.empty(nodes,dtype=np.float64) # TODO: Confirm number of rows == # nodes expected
for rowId,omega,init in cursor:
    omegas[rowId-1] = omega
    x0[rowId-1] = float(init)


## Build the network matrix
networkMatrix = csr_matrix((np.ones(edges), (sources, dests)),\
                           shape = (nodes, nodes),\
                           dtype = bool)


# cleanup 
#del(sources)
#del(dests)
#del(cursor)
#del(data)
## Force GC
#import gc
#gc.collect()

endTload = time.process_time() # Do you want to ignore data load time in comparison?
#print(f'Loading {nodes} nodes and {edges} edges from SQLite takes {finLoad-startLoad}s')

# T load and jit
startTjit = time.process_time()

def kuramotos_f():
	for i in range(n):
		coupling_sum=sum(
			sin(y(j)-y(i))
			for j in range(n)
			if networkMatrix[j,i]
		)
		yield omegas[i]+couplingConstant*coupling_sum


from jitcode import jitcode
n = nodes
I=jitcode(kuramotos_f,n=nodes)
I.set_integrator(integrator, atol=atol, rtol=rtol)#, nsteps=200)

endTjit = time.process_time()
initial_state=x0
I.set_initial_value(initial_state, time=0.0)
dt = 2
times = np.arange(0,float(tend)+1,dt)

states = []
startTintegrate = time.process_time()
for t in times:
	states.append(I.integrate(t)% (2*np.pi))
endTintegrate = time.process_time()

tstartup = endTstartup - startTstartup
tload = endTload - startTload
tjit = endTjit - startTjit
tintegrate = endTintegrate - startTintegrate
ttotal = tstartup+tload+tjit+tintegrate

## Write the experiment information into the database
githash = os.popen('git rev-parse HEAD').read()
gitchanges = os.popen('git status -s -uno').read()
hostname = os.popen('hostname').read()
datetime = os.popen('date').read()
commandline = ' '.join(sys.argv)

conn.execute("insert into experiments \
             (runtime,datetime,githash,gitchanges,commandline,host,\
             system,solver,atol,rtol,tstart,tend,\
             tstartup,tload,tjit,tintegration,ttotal)\
             values"+ \
             f'(\"jitcode\",\
             \"{datetime}\",\"{githash}\",\"{gitchanges}\",\
             \"{commandline}\",\"{hostname}\",\
             {systemid},\"{integrator}\",{atol},{rtol},0,{tend},\
             {tstartup},{tload},{tjit},{tintegrate},{ttotal})')
### Get the newly created unique experimentid
cursor = conn.execute(f'select last_insert_rowid()')
experimentid = cursor.fetchone()[0]
## Write the state vector out to the database.
for i,row in enumerate(states):
    for k,state in enumerate(row):
        conn.execute(f'insert into states (experimentid, time,idx,node,state)\
                 values\
                 ({experimentid},{times[i]},{k+1},{k+1},{state})')
conn.commit()

## Finally, find a reference trajectory for this trajectory and compute the 
refid=-1
cursor = conn.execute(f'select experimentid from experiments \
                      where \
                        system={systemid} and \
                        tend={tend} and \
                        solver is "radau" and \
                        rtol=1e-10 and \
                        atol=1e-12;')

try:
    try:
        refeid = cursor.fetchone()[0]
        errors = np.empty(len(states[-1]))
        for idx,state in enumerate(states[-1]):
            cursor = conn.execute(f'select state from states where \
                                   experimentid={refeid} \
                                   and time={tend} \
                                   and idx={idx+1};')
            refstate = cursor.fetchone()[0]
            errors[idx] = abs(refstate - state)

        avgErr = np.mean(errors)
        conn.execute(f'update experiments \
                     set \
                        err_v_ref={avgErr}, \
                        refeid={refeid} \
                     where \
                        experimentid={experimentid};')
        conn.commit()
    except Exception as e:
        print("Comparison to reference trajectory failed with exception: ")
        print(e)
        pass

except Exception:
    print("No reference id found!")

print(f'Finished system {systemid}')
