module SIDISRecoveryLocations
using JSON, SHA

"Only immutable path strings and source indices; no model results or physical queries."
struct ScopedLocations
    root::String
    rows::Tuple
    signature::String
end

function snapshot(run,paths)
    paths.valid(run)
    rows=Tuple(Tuple(row[key] for key in paths.FIELDS) for row in run.entries)
    ScopedLocations(run.root,rows,run.signature)
end

function retained(run::ScopedLocations,phase,index,paths)
    paths.point_index(phase,index)
    for row in run.rows
        row[1]==phase && row[2]==index && return Dict{String,Any}(key=>value for (key,value) in zip(paths.FIELDS,row))
    end
    nothing
end

function current_path(run::ScopedLocations,parts::AbstractString...)
    path=joinpath(run.root,parts...)
    normpath(path)==path && (path==run.root || startswith(path,run.root*"/")) || error("current-run path escapes its immutable root")
    path
end

function folder(run::ScopedLocations,phase,index,paths)
    if phase in paths.POINT_PHASES
        row=retained(run,phase,index,paths)
        return row===nothing ? current_path(run,phase,"p"*lpad(string(index),4,'0')) : row["folder"]
    end
    phase in paths.GLOBAL_PHASES && index isa Int && index==0 || error("declared fresh global phase/index0 required")
    current_path(run,phase)
end

function fresh_folder(run::ScopedLocations,phase,index,paths)
    path=folder(run,phase,index,paths)
    startswith(path,run.root*"/") && !ispath(path) && !islink(path) || error("never dispatch retained or occupied outputs")
    parent=dirname(path)
    while !ispath(parent) && !islink(parent);parent=dirname(parent);end
    realpath(parent)==parent && !islink(parent) || error("new output cannot traverse aliases")
    path
end

function verify_boundary(run::ScopedLocations,paths)
    entries=[Dict{String,Any}(key=>value for (key,value) in zip(paths.FIELDS,row)) for row in run.rows]
    restored=paths.RecoveryRunPaths(run.root,entries)
    restored.signature==run.signature || error("changed immutable location snapshot")
    restored
end
end
