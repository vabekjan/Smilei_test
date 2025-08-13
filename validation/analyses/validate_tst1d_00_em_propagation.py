import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)


# FIELD DIAG
Validate("List of fields in Field0", S.fieldInfo(0)["fields"] )

timesteps = list(S.Field.Field0("Ez").getAvailableTimesteps())
Validate("List of timesteps in Field0", timesteps )

Ez = S.Field.Field0("Ez", timesteps=timesteps[-1]).getData()[0]
Validate("Last Ez profile in Field0", Ez, 1e-7 )

# SUBGRID OF FIELD DIAG
Ez1 = S.Field.Field1("Ez", timesteps=timesteps[-1]).getData()[0]
subgrid = S.namelist.DiagFields[1].subgrid
Validate("Field subgrid works", (Ez1==Ez[subgrid]).all())

# TIME-AVERAGED FIELD DIAG
D = S.Field.Field2("Ez")
timesteps = list(D.getAvailableTimesteps())
Ez = D.getData(timesteps[-1])[0]
Validate("Last Ez profile in Field2", Ez, 1e-7 )

# 0-D PROBE
Ez = S.Probe(0,"Ez").getData()
Validate("Ez vs time in Probe0", Ez, 1e-7 )

# 1-D PROBE
Validate("List of fields in Probe1", S.probeInfo(1)["fields"] )

timesteps = S.Probe(1,"Ez").getAvailableTimesteps()
Validate("List of timesteps in Probe1", timesteps )

Ez = S.Probe(1,"Ez", timesteps=timesteps[-1]).getData()[0]
Validate("Last Ez profile in Probe1", Ez, 1e-7 )

# UBAL SCALAR
max_ubal = np.max( np.abs(S.Scalar.Ubal().getData()) )
Validate("Max Ubal is below 3%", max_ubal<0.03 )

# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5", "r") as f:
	dt = f["data/0"].attrs["dt"]
	dx = f["data/0/Ex"].attrs["gridSpacing"]
	patchSize = f["data/0"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
