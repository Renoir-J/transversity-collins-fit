# Add a method only for our new path type. The frozen generic method is retained.
function HistoricalSIDISJobRecords.point_folder(run::SIDISContinuationPaths.ContinuationRun,phase,index)
    SIDISContinuationPaths.point_folder(run,phase,index)
end
