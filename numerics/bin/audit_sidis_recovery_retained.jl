using JSON
include("SIDISRecoveryBindings.jl")
const RB=SIDISRecoveryBindings
length(ARGS)==5 || error("FROZEN PREDECESSOR MANIFEST MANIFEST_SHA FRESH_AUDIT")
const RM=RB.load_manifest(ARGS[1:4]...)
include(joinpath(RM.frozen,"numerics/bin/HistoricalSIDISJobRecords.jl"))
include("SIDISRecoveryDispatch.jl")
include("SIDISRetainedPairReview.jl")
exit(SIDISRetainedPairReview.run(RM,ARGS[5],HistoricalSIDISJobRecords,RB,check_recovery_transport))
