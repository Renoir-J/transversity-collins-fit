# Specialize only the two owned run types. All pre-existing methods remain intact.
const RecoveryBindings=SIDISRecoveryBindings
function HistoricalSIDISJobRecords.point_folder(run::Union{RecoveryBindings.RecoveryRun,RecoveryBindings.LegacyRun},phase,index)
    RecoveryBindings.point_folder(run,phase,index)
end

function HistoricalSIDISJobRecords.point_dependencies(ctx,run::RecoveryBindings.RecoveryRun,phase,index)
    producer=RecoveryBindings.origin(run,phase,index)
    producer===run ? invoke(HistoricalSIDISJobRecords.point_dependencies,Tuple{Any,Any,Any,Any},ctx,run,phase,index) :
        HistoricalSIDISJobRecords.point_dependencies(ctx,producer,phase,index)
end

function check_recovery_transport(run,phase,index,packet;complete_invocation=false)
    expected=RecoveryBindings.transport(run,phase,index)
    folder=RecoveryBindings.point_folder(run,phase,index)
    if expected===nothing
        !haskey(packet,"recovery") && !haskey(packet,"continuation") || error("original packet cannot claim a later transport")
    else
        key=expected["scope"]==RecoveryBindings.SCOPE ? "recovery" : "continuation"
        get(packet,key,nothing)==expected || error("unbound actual producing transport")
        if complete_invocation
            for name in ("invocation.json","audit.json")
                get(JSON.parsefile(joinpath(folder,name)),key,nothing)==expected || error("changed complete producer/auditor provenance")
            end
        end
    end
    true
end

function check_current_recovery_prerequisites(ctx,run,phase,packet)
    if phase=="native4"
        current=HistoricalSIDISJobRecords.selection(ctx,run,"native-select3")
        current["status"]=="NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED" && current["selected_level"]===0 || error("native4 not requested by the current full132 decision")
    elseif phase in ("tailref123","reference4")
        current=HistoricalSIDISJobRecords.selection(ctx,run,"native-selection")
        current["selected_level"]==packet["native_level"] || error("retained job disagrees with the current common native level")
        if phase=="reference4"
            derivative=HistoricalSIDISJobRecords.selection(ctx,run,"derivative-select3")
            derivative["status"]=="DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED" && derivative["selected_reference_level"]===0 || error("reference4 not requested by the current full132 decision")
        end
    end
    true
end

function HistoricalSIDISJobRecords.audit_point(ctx,run::Union{RecoveryBindings.RecoveryRun,RecoveryBindings.LegacyRun},phase,index;complete_invocation=false)
    packet=invoke(HistoricalSIDISJobRecords.audit_point,Tuple{Any,Any,Any,Any},ctx,run,phase,index;complete_invocation)
    check_recovery_transport(run,phase,index,packet;complete_invocation)
    check_current_recovery_prerequisites(ctx,run,phase,packet)
    packet
end
