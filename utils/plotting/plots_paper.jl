using Pkg
Pkg.activate(string(@__DIR__, "/../../src/julia/"))

using Plots
using SQLite
using DataFrames

db = SQLite.DB(string(@__DIR__,"/../../data/KuraBenchXeon.db"))
_df = DataFrame(DBInterface.execute(db, "select * from experiments"))

df = _df[_df.solver .!= "radau", :]

# The database contains experiments on 30 systems of coupled ODEs
# Each system is defined on a unique network topology
# Systems 1:10  have 10 nodes
# Systems 11:20 have 100 nodes
# Systems 21:30 have 1000 nodes
# Every system is simulated with 10 different solver tolerances [0.001, 0.0001, ..., 10e-12]
# In total we end up with 300 different experiments
# These are repeated for 5 runtimes [scipy, jl_sparse, jl_nd, jitcode, Fortran]
#
# The reported results are averages of all systems with the same solver tolerances


const ntols = 10 # number of different tolerance settings


"""
      function system_size_mean(df, system_ids, runtime, column)

      Selects all rows in df that match the specified `runtime` and one of the `system_ids` and computes their mean of the given `column` with respect to a fixed solver tolerance
"""
function system_size_mean(df, system_ids, runtime, column)
      nsystems = length(system_ids)
      rows = df[(df.runtime .== runtime) .& map(x -> x ∈ system_ids, df.system), :]
      return sum(reshape(rows[:,column], ntols, nsystems), dims = 2) ./ nsystems
end

# Simple Test. Average of tolerance should give back the same tolerance.
@assert system_size_mean(df, 1:10, "jitcode", "rtol") ≈ df[1:10,:].rtol


# We produce different plots for each system size
ids_per_size = [1:10, 11:20, 21:30]
i = 1
saving = true

## Error vs tolerance
begin
      plot(system_size_mean(df, ids_per_size[i], "jitcode", "rtol"),
           system_size_mean(df, ids_per_size[i], "jitcode", "err_v_ref"),
           xlabel = "Requested tolerance",
           ylabel = "Error",
           scale = :log10,
           label = "jitcode",
           color_palette = :seaborn_colorblind6,
           legend = :bottomright)

      for runtime in ["networkdynamics.jl", "DifferentialEquations.jl", "scipy", "Fortran"]
            Plots.display(
                  plot!(system_size_mean(df, ids_per_size[i], runtime, "rtol"),
                        system_size_mean(df, ids_per_size[i], runtime, "err_v_ref"),
                        label = runtime))
      end
end


## WPD


for i = 1:3
      begin
            plot(system_size_mean(df, ids_per_size[i], "jitcode", "err_v_ref"),
                 system_size_mean(df, ids_per_size[i], "jitcode","tintegration"),
                 xlabel = "Error",
                 ylabel = "Integration time (s)",
                 title = "N = $(10^i)",
                 scale = :log10,
                 marker = :diamond,
                 markersize = 5,
                 label = "JiTCODE",
                 color_palette =:seaborn_colorblind6,
                 legend = :topright)

            runtime = "networkdynamics.jl"
            plot!(system_size_mean(df, ids_per_size[i], runtime, "err_v_ref"),
                  system_size_mean(df, ids_per_size[i], runtime, "tintegration"),
                  marker = :star5,
                  markersize = 5,
                  label = "NetworkDynamics.jl")

            runtime = "DifferentialEquations.jl"
            plot!(system_size_mean(df, ids_per_size[i], runtime, "err_v_ref"),
                  system_size_mean(df, ids_per_size[i], runtime, "tintegration"),
                  marker = :utriangle,
                  markersize = 5,
                  label = "SparseArrays.jl")

            runtime = "Fortran"
            Plots.display(plot!(system_size_mean(df, ids_per_size[i], runtime, "err_v_ref"),
                        system_size_mean(df, ids_per_size[i], runtime, "tintegration"),
                        marker = :x,
                        markersize = 5,
                        label = runtime))

            if i == 1
                  runtime = "scipy"

                  Plots.display(plot!(system_size_mean(df, ids_per_size[i], runtime, "err_v_ref"),
                              system_size_mean(df, ids_per_size[i], runtime, "tintegration"),
                              marker = :o,
                              markersize = 5,
                              label = "SciPy"))
            end
            saving ? savefig(string(@__DIR__, "/WPD$(10^i).pdf")) : nothing
      end
end
## WPD jit

begin

      plot(system_size_mean(df, ids_per_size[i], "jitcode", "err_v_ref"),
            system_size_mean(df, ids_per_size[i], "jitcode","tjit"),
            xlabel = "Error",
            ylabel = "Preparation + integration time (s)",
            scale = :log10,
            label = "JiTCODE",
            color_palette = :seaborn_colorblind6,
            legend = :topright)


      runtime = "networkdynamics.jl"
      plot!(system_size_mean(df, ids_per_size[i], runtime, "err_v_ref"),
            system_size_mean(df, ids_per_size[i], runtime, "tjit"),
            label = "NetworkDynamics.jl")

      runtime = "DifferentialEquations.jl"
      plot!(system_size_mean(df, ids_per_size[i], runtime, "err_v_ref"),
            system_size_mean(df, ids_per_size[i], runtime, "tjit"),
            label = "SparseArrays.jl")
end



# Scaling of compute time vs system size
# For Jitcode the inital improvement is porbably due to SymEngine replacing SymPy (?),
# the later increase due to increased cost of symbolic preprocessing

begin

      plot([10, 100, 1000],
           [system_size_mean(df, sys_ids, "jitcode","ttotal")[4] for sys_ids in ids_per_size],
           xlabel = "Number of nodes",
           ylabel = "Preparation + integration time (s)",
           scale = :log10,
           color_palette = :seaborn_colorblind6,
           label = "JiTCODE",
           marker = :diamond,
           markersize = 5,
           legend = :topleft)

      runtime = "networkdynamics.jl"
      plot!([10, 100, 1000],
            [system_size_mean(df, sids, runtime,"ttotal")[4] for sids in ids_per_size],
            marker = :star5,
            markersize = 5,
            label = "NetworkDynamics.jl")

      runtime = "DifferentialEquations.jl"
      plot!([10, 100, 1000],
            [system_size_mean(df, sids, runtime,"ttotal")[4] for sids in ids_per_size],
            marker = :utriangle,
            markersize = 5,
            label = "SparseArrays.jl")

      Plots.display(
            plot!([10, 100, 1000],
                  [system_size_mean(df, sids, "scipy","ttotal")[4] for sids in ids_per_size],
                  marker = :circle,
                  markersize = 5,
                  label = "SciPy"))


      saving ?  savefig(string(@__DIR__, "/size_vs_jit.pdf")) : nothing
end
