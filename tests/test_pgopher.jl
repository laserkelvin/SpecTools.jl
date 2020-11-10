
using Pkg

Pkg.activate("SpecTools")
using SpecTools
using BSON: @save

mcta = PGopher("mcta", "methylcyanotriacetylene.pgo", 8)

# calculate partition function
print(Q(mcta.levels, 300.))

@save "mcta.bson" mcta

