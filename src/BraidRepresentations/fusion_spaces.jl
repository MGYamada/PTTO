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

_object_label(x::Integer) = Int(x)
_object_label(x::Symbol) = x === :one || x === Symbol("1") ? 1 :
    x === :τ || x === :tau || x === :Tau ? 2 :
    x === :s ? 2 :
    x === :σ || x === :sigma ? 2 :
    x === :ψ || x === :psi ? 3 :
    error("unsupported symbolic simple-object label $x; pass integer labels or one of :one, :τ/:tau, :s, :σ/:sigma, :ψ/:psi")
_object_labels(xs) = [_object_label(x) for x in xs]

function fusion_paths(rules, objects::AbstractVector, total)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    objs = _object_labels(objects)
    tgt = _object_label(total)
    foreach(a -> _check_object(fr, a), objs)
    _check_object(fr, tgt)
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

fusion_basis(rules, objects::AbstractVector, total) =
    FusionTreeBasis(_fusion_rule(rules), _object_labels(objects), _object_label(total),
                    fusion_paths(rules, objects, total))
