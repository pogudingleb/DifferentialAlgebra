using Singular

include("DifferentialPolynomials.jl")

# Should not be here?
AbstractAlgebra.needs_parentheses(x::Singular.n_unknown{Nemo.fmpq}) = false
AbstractAlgebra.needs_parentheses(x::Singular.n_unknown{AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}}) = true

#------------------------------------------------------------------------------

function lex_ranking(varnames::Array{String, 1})
    return (x, y) -> begin
        vx, ordx = decompose_derivative(x)
        vy, ordy = decompose_derivative(y)
        indx = findfirst(a -> a == vx, varnames)
        indy = findfirst(a -> a == vy, varnames)
        if indx != indy
            return indx < indy
        end
        return ordx < ordy
    end
end

function lex_ranking(R::DiffPolyRing)
    return lex_ranking(R.varnames)
end

#------------------------------------------------------------------------------

function get_singular_polys(diff_polys::Array{DiffPoly, 1}, ordering=:lex, ranking=nothing)
    R = parent(diff_polys[1])
    if ranking === nothing
        ranking = lex_ranking(R)
    end

    orders = Array{Integer, 1}()
    for v in gens(R)
        push!(orders, max([order(p, v) for p in diff_polys]...))
    end
    newvarnames = Array{String, 1}()
    for (i, name) in enumerate(R.varnames)
        for ord in 0:orders[i]
            push!(newvarnames, form_derivative(name, ord))
        end
    end
    sort!(newvarnames, lt=ranking)
    bring = base_ring(R.poly_ring)
    Rsing, varsing = Singular.PolynomialRing(bring, newvarnames, ordering=ordering)
    return Rsing, [parent_ring_change(algdata(p), Rsing) for p in diff_polys]
end

function from_singular(p, R::DiffPolyRing)
    return DiffPoly(R, parent_ring_change(p, R.poly_ring))
end

#------------------------------------------------------------------------------
