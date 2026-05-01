struct FusionPath
    objects::Vector{Int}
    total::Int
    channels::Vector{Int}
end

struct FusionTreeBasis
    rules::FusionRule
    objects::Vector{Int}
    total::Int
    paths::Vector{FusionPath}
end

Base.:(==)(a::FusionPath, b::FusionPath) =
    a.objects == b.objects && a.total == b.total && a.channels == b.channels
Base.hash(p::FusionPath, h::UInt) = hash((p.objects, p.total, p.channels), h)

Base.length(b::FusionTreeBasis) = length(b.paths)
Base.iterate(b::FusionTreeBasis, st...) = iterate(b.paths, st...)
dim(b::FusionTreeBasis) = length(b)

function _fusion_space_object_index(fr::FusionRule, x)
    x isa Integer || error("fusion-space object labels must be integers in 1:$(fr.rank); got $(repr(x))")
    idx = Int(x)
    _check_object(fr, idx)
    return idx
end

function fusion_paths(rules, objects::AbstractVector, total)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    objs = [_fusion_space_object_index(fr, x) for x in objects]
    tgt = _fusion_space_object_index(fr, total)
    isempty(objs) && return FusionPath[]
    length(objs) == 1 && return objs[1] == tgt ? [FusionPath(objs, tgt, [tgt])] : FusionPath[]
    paths = FusionPath[]
    function extend(pos::Int, current::Int, channels::Vector{Int})
        if pos > length(objs)
            current == tgt && push!(paths, FusionPath(copy(objs), tgt, copy(channels)))
            return
        end
        for c in fusion_channels(fr, current, objs[pos])
            extend(pos + 1, c, [channels; c])
        end
    end
    for c in fusion_channels(fr, objs[1], objs[2])
        extend(3, c, [c])
    end
    sort!(paths; by = p -> Tuple(p.channels))
    return paths
end

function fusion_basis(rules, objects::AbstractVector, total)
    fr = _fusion_rule(rules)
    return FusionTreeBasis(fr, [_fusion_space_object_index(fr, x) for x in objects],
                           _fusion_space_object_index(fr, total),
                           fusion_paths(fr, objects, total))
end
