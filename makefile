

database:
	@echo "Building database.."
	@sqlite3 ${KURABENCH_DB} < ./data/KuraBench.sql

systems:
	@echo "Generating systems.."
	@KURABENCH_DB=${KURABENCH_DB} python3 ./utils/generate_systems.py

reference:
	@echo "Computing reference trajectories"
	@KURABENCH_DB=${KURABENCH_DB} ./utils/gen_ref.sh

jitcode:
	@echo "Computing jitcode trajectories"
	@KURABENCH_DB=${KURABENCH_DB} ./utils/run_jitcode.sh

scipy:
	@echo "Computing scipy trajectories"
	@KURABENCH_DB=${KURABENCH_DB} ./utils/run_scipy.sh

networkdynamics:
	@echo "Computing networkdynamics trajectories"
	@KURABENCH_DB=${KURABENCH_DB} ./utils/run_nd.sh

diffeq:
	@echo "Computing pure diffeq.jl trajectories"
	@KURABENCH_DB=${KURABENCH_DB} ./utils/run_diffeq.sh

fortran:
	@echo "Computing fortran trajectories"
	@KURABENCH_DB=${KURABENCH_DB} ./utils/run_f.sh

clean:
	@rm ${KURABENCH_DB} || echo "Nothing to cleanup.."
