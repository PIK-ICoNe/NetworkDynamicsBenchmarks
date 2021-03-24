systemid=parse(Int8,ARGS[1])
tend=parse(Float64,ARGS[2])
atol=parse(Float64,ARGS[3])
rtol=parse(Float64,ARGS[4])

using LinearAlgebra: BLAS
BLAS.set_num_threads(1) # Recording CPU time means we should
                        # stick to one thread

# Load data
load_data = :(begin
    using SQLite
    db = SQLite.DB(ENV["KURABENCH_DB"])

    # First read the system info
    results = DBInterface.execute(db,string(
            "select ",
                "name, nodes, edges, coupling_constant " ,
            "from ",
                "systems ",
            "where " ,
                "id=$systemid;"))
    name = ""
    for row in results
        global name = row.name
        global nodes = row.nodes
        global edges = row.edges
        global coupling_constant = row.coupling_constant
    end
    # Now read connectivity
    # Allocate
    sources = Array{Int32}(undef,edges)
    dests = Array{Int32}(undef,edges)
    data = ones(Int32,edges)

    results = DBInterface.execute(db,string(
            "select ",
                "id,source,dest ",
             "from " ,
                "edges " ,
             "where " ,
                "system=$systemid;"))

    for (i,row) in enumerate(results)
        sources[i] = row.source
        dests[i]   = row.dest
    end
    # Store the node omegas
    omegas = Array{Float64}(undef,nodes)
    x0     = Array{Float64}(undef,nodes)
    results = DBInterface.execute(db,string(
            "select ",
                "id,omega,initialCondition ",
            "from " ,
                "nodes " ,
            "where " ,
                "system=$systemid;"))
    for (i,row) in enumerate(results)
        omegas[i] = row.omega
        x0[i] = row.initialCondition
    end
    my_matrix = sparse(sources,dests,data)
end)

prepare_ode = :(begin
    g = SimpleGraph(my_matrix')
    integrator = DP5
    ### Network dynamics vertex and edge functions

	@inline Base.@propagate_inbounds function kuramoto_vertex!(dv, v, edges, p, t)
		dv[1] = p # parameter for vertices are eigenfrequencies omega
		@inbounds for e in edges
	            dv[1] += e[1]
	    end
	    nothing
	end

    @inline Base.@propagate_inbounds function kuramoto_edge!(e,v_s,v_d,p,t)
        e[1] = p * sin(v_s[1] - v_d[1])
        nothing
    end

    ### Constructing the network dynamics

    odevertex = ODEVertex(f! = kuramoto_vertex!, dim = 1)
    staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1, coupling = :antisymmetric)

    # generating random values between -0.5 and 0.5 for the parameter value ω_0 of the vertices
    e_pars = [coupling_constant for i in 1:g.ne]
    parameters = (omegas, e_pars)

    # setting up the  network_dynamics
    kuramoto_network! = network_dynamics(odevertex, staticedge, g)

    ### Simulation
    # constructing random initial conditions for nodes (variable θ)

    prob = ODEProblem(kuramoto_network!, x0, (0.,tend), parameters, abstol=atol, reltol=rtol)

end)

# Precompile eval
eval(:nothing)

using Pkg
Pkg.activate(@__DIR__)
println("Enviroment Activated, starting program...")
using CPUTime

# Times
tstartup = @CPUelapsed begin
	using NetworkDynamics: ODEVertex, StaticEdge, network_dynamics
	using LightGraphs: SimpleGraph
	using OrdinaryDiffEq: ODEProblem, solve, DP5
    using SparseArrays
end
times = 0.:1.:tend
println("time to startup is $tstartup")
tload    = @CPUelapsed(eval(load_data))
println("time to load is $tload")
tjit     = @CPUelapsed(eval(prepare_ode)) +
           @CPUelapsed(solve(prob, integrator(), saveat=times))

println("time to jit is $tjit")
tintegrate = @CPUelapsed(sol = solve(prob, integrator(), saveat=times))

println("time to solve is $tintegrate")

ttotal = tstartup + tload + tjit + tintegrate
# Now write the config to the experiments table and the states to the states table
#

## Config
tidy = cmd -> chomp(String(read(cmd)))

githash = tidy(`git rev-parse HEAD`)
hostname = tidy(`hostname`)
datetime = tidy(`date`)
gitchanges = String(read(`git status -s -uno`))
argsstring = join(ARGS," ")
intstring = string(integrator)
DBInterface.execute(db,string(
    "insert into experiments ",
    "(runtime,datetime,githash,gitchanges,commandline,host,",
    "system,solver,atol,rtol,tstart,tend,",
    "tstartup,tload,tjit,tintegration,ttotal)",
    "values ",
    "(\"networkdynamics.jl\",",
    "\"$datetime\",\"$githash\",\"$gitchanges\",\"$argsstring\",\"$hostname\",",
    "$systemid,\"$intstring\",$atol,$rtol,0,$tend,",
    "$tstartup,$tload,$tjit,$tintegrate,$ttotal);"
                              )
                   )

results = DBInterface.execute(db,"select last_insert_rowid()")
experimentid = 0
for row in results
    # TODO: How to get a single element from julia iterator?
    global experimentid = row[1]
end

for (i,t) in enumerate(sol.t)
    for (k,state) in enumerate(sol[i])
        state = mod(state,2*pi)
        DBInterface.execute(db,string(
            "insert into states",
            "(experimentid,time,idx,node,state)",
            "values",
            "($experimentid,$t,$k,$k,$state);"
                                     )
                           )
    end
end

# TODO: Get the reference trajectory for this request (if it exists)
# and update the expermients table (see the python_complete file for guide)
results =DBInterface.execute(db,string(
            "select experimentid from experiments where",
            " system=\"$systemid\"",
            " and tend=$tend",
            " and atol=1e-12",
            " and rtol=1e-10",
            " and solver is \"radau\";"
                                     )
                   )
refid = -1
for row in results
    global refid=row.experimentid
    break
end
errors = zeros(length(sol[end]))
for (idx,state) in enumerate(last(sol))
    local results = DBInterface.execute(db,string(
            "select state from states",
            " where experimentid=$refid",
            " and time=$tend",
            " and idx=$idx;"
                                           )
                                 )
    for row in results
        global errors[idx] = abs(row.state-mod(state,2*pi))
    end
end

avgErr = sum(errors) / length(errors)

DBInterface.execute(db,string(
        "update experiments",
        " set err_v_ref=$avgErr,",
        " refeid=$refid",
        " where experimentid=$experimentid;"
                             )
                   )
