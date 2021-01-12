import Base
import AbstractAlgebra, Nemo
#import Oscar

include("utils.jl")
include("AbstractTypes.jl")

#------------------------------------------------------------------------------

function form_derivative(varname::String, order::Integer)
    return "$(varname)^($order)"
end

function decompose_derivative(dername::String)
    regex = r"^(.*)\^\((\d+)\)$"
    m = match(regex, dername)
    if m === nothing
        throw(DomainError("$dername is not a name of a derivative"))
    end
    return (String(m.captures[1]), parse(Int, m.captures[2]))
end

#------------------------------------------------------------------------------

mutable struct DiffPolyRing <: DiffRing
    base_ring::AbstractAlgebra.Ring
    poly_ring::MPolyRing
    max_ord::Integer
    varnames::Array{String, 1}
    derivation::Dict{Any, Any}

    function DiffPolyRing(R::AbstractAlgebra.Ring, varnames::Array{String, 1}, max_ord=20)
        all_varnames = [form_derivative(v, ord) for v in varnames for ord in 0:max_ord]
        poly_ring, _ = AbstractAlgebra.PolynomialRing(R, all_varnames)
        derivation = Dict()
        for v in varnames
            for ord in 0:(max_ord - 1)
                derivation[str_to_var(form_derivative(v, ord), poly_ring)] = 
                    str_to_var(form_derivative(v, ord + 1), poly_ring)
            end
        end
        return new(R, poly_ring, max_ord, varnames, derivation)
    end
end

function DifferentialPolynomialRing(R::AbstractAlgebra.Ring, varnames::Array{String, 1}, max_ord=20)
    R = DiffPolyRing(R, varnames, max_ord)
    return R, Tuple(gens(R))
end

#------------------------------------------------------------------------------

mutable struct DiffPoly <: DiffRingElem
    parent::DiffPolyRing
    algdata::AbstractAlgebra.MPolyElem

    function DiffPoly(R::DiffPolyRing, alg_poly::AbstractAlgebra.MPolyElem)
        return new(R, parent_ring_change(alg_poly, R.poly_ring))
    end
end

function Base.parent(a::DiffPoly)
    return a.parent
end

function algdata(a::DiffPoly)
    return a.algdata
end

AbstractAlgebra.elem_type(::Type{DiffPolyRing}) = DiffRingElem

AbstractAlgebra.parent_type(::Type{DiffRingElem}) = DiffPolyRing

AbstractAlgebra.parent_type(::Type{DiffPoly}) = DiffPolyRing

function get_der(R::DiffPolyRing, vname::String, ord::Integer)
    if ord > R.max_ord
        throw(DomainError("Requested order $ord exceed the maximal order of the ring $(R.max_ord)"))
    end
    return DiffPoly(R, str_to_var(form_derivative(vname, ord), R.poly_ring))
end

#------------------------------------------------------------------------------

function AbstractAlgebra.gens(R::DiffPolyRing)
    return [DiffPoly(R, str_to_var(form_derivative(v, 0), R.poly_ring)) for v in R.varnames]
end

#------------------------------------------------------------------------------

mutable struct DiffRational <: DiffFieldElem
    parent::DiffPolyRing
    algdata::Generic.Frac{<: AbstractAlgebra.MPolyElem}

    function DiffRational(R::DiffPolyRing, alg_rational::Generic.Frac{<: AbstractAlgebra.MPolyElem})
        p, q = map(x -> parent_ring_change(x, R.poly_ring), [numerator(alg_rational), denominator(alg_rational)])
        return new(R, p // q)
    end
end

function Base.parent(a::DiffRational)
    return a.parent
end

function algdata(a::DiffRational)
    return a.algdata
end

function Base.numerator(a::DiffRational)
    return DiffPoly(parent(a), numerator(algdata(a)))
end

function Base.denominator(a::DiffRational)
    return DiffPoly(parent(a), denominator(algdata(a)))
end

AbstractAlgebra.parent_type(::Type{DiffRational}) = DiffPolyRing

#------------------------------------------------------------------------------

function (R::DiffPolyRing)(b)
    return DiffPoly(R, R.poly_ring(b))
end

function (R::DiffPolyRing)(b::Generic.Frac)
    return DiffRational(R, b)
end

function (R::DiffPolyRing)(b::DiffPoly)
    return DiffPoly(R, algdata(b))
end

function (R::DiffPolyRing)(b::DiffRational)
    return DiffRational(R, algdata(b))
end

function (R::DiffPolyRing)()
    return zero(R)
end

function Base.zero(R::DiffPolyRing)
    return R(0)
end

function Base.one(R::DiffPolyRing)
    return R(1)
end

#------------------------------------------------------------------------------

function Base.show(io::IO, R::DiffPolyRing)
    print(io, "Differential polynomial ring over $(R.base_ring) in " * join(R.varnames, ", "))
end

function Base.show(io::IO, p::DiffRingElem)
    show(io, algdata(p))
end

#------------------------------------------------------------------------------

function check_parent(a::DiffRingElem, b::DiffRingElem)
    if parent(a) != parent(b)
        throw(DomainError("$a and $b are from different rings: $(parent(a)) and $(parent(b)), respectively"))
    end
end

#------------------------------------------------------------------------------

function Base.iszero(a::DiffRingElem)
    return algdata(a) == 0
end

#------------------------------------------------------------------------------

function Base.:+(a::DiffRingElem, b::DiffRingElem)
    check_parent(a, b)
    return parent(a)(algdata(a) + algdata(b))
end

function Base.:+(a::DiffRingElem, b)
    return parent(a)(algdata(a) + b)
end

function Base.:+(a, b::DiffRingElem)
    return parent(b)(algdata(b) + a)
end

function Base.:+(a::DiffRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational, Nemo.fmpq})
    return parent(a)(algdata(a) + b)
end

function AbstractAlgebra.addeq!(a::DiffRingElem, b::DiffRingElem)
    a = a + b
end

#------------------------------------------------------------------------------

function Base.:*(a::DiffRingElem, b::DiffRingElem)
    check_parent(a, b)
    return parent(a)(algdata(a) * algdata(b))
end

function Base.:*(a::RingElem, b::DiffRingElem)
    if typeof(a) <: DiffRingElem
        return parent(a)(algdata(a) * algdata(b))
    end
    return parent(b)(a * algdata(b))
end

#function Base.:*(a, b::DiffRingElem)
#    return parent(b)(algdata(b) * a)
#end

function Base.:*(a::Union{AbstractFloat, Integer, Rational}, b::DiffRingElem)
    return parent(b)(a * algdata(b))
end

function AbstractAlgebra.mul!(a::DiffRingElem, b::DiffRingElem, c::DiffRingElem)
    a = b * c
end

#------------------------------------------------------------------------------

function Base.:-(a::DiffRingElem, b::DiffRingElem)
    check_parent(a, b)
    return parent(a)(algdata(a) - algdata(b))
end

function Base.:-(a::DiffRingElem, b)
    return parent(a)(algdata(a) - b)
end

function Base.:-(a, b::DiffRingElem)
    return parent(b)(-algdata(b) + a)
end

function Base.:-(a::DiffRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational, Nemo.fmpq})
    return parent(a)(algdata(a) - b)
end

#------------------------------------------------------------------------------

function Base.:-(a::DiffRingElem)
    return parent(a)(-algdata(a))
end

#------------------------------------------------------------------------------

function Base.:^(a::DiffRingElem, i::Integer)
    return parent(a)(algdata(a)^i)
end

#------------------------------------------------------------------------------

function Base.:(==)(a::DiffRingElem, b::DiffRingElem)
    return algdata(a) == algdata(b)
end

#------------------------------------------------------------------------------

function Base.:(//)(a::DiffRingElem, b::DiffRingElem)
    check_parent(a, b)
    return parent(a)(algdata(a) // algdata(b))
end

function Base.:(//)(a::DiffRingElem, b)
    return parent(a)(algdata(a) // b)
end

function Base.:(//)(a, b::DiffRingElem)
    return parent(b)(a // algdata(b))
end

# Fixing ambiguity
function Base.:(//)(a::Union{Integer, Rational}, b::DiffRingElem)
    return parent(b)(a // algdata(b))
end

#------------------------------------------------------------------------------

function Base.hash(a::DiffRingElem)
    return hash(algdata(a))
end

#------------------------------------------------------------------------------

function d_aux(p::MPolyElem, der::Dict{Any, Any})
    result = zero(parent(p))
    for v in vars(p)
        if !(v in keys(der))
            throw(DomainError("No derivative defined for $v. Most likely you have exceeded the maximal order."))
        end
        result += der[v] * derivative(p, v)
    end
    return result
end

function d(a::DiffPoly)
    return DiffPoly(parent(a), d_aux(algdata(a), parent(a).derivation))
end

function d(a::DiffRational)
    p, q = numerator(algdata(a)), denominator(algdata(a))
    dp, dq = map(x -> d_aux(x, parent(a).derivation), [p, q])
    return DiffRational(parent(a), (dp * q - p * dq) // q^2)
end

function d(a::DiffRingElem, ord::Integer)
    if ord == 0
        return a
    end
    return d(d(a), ord - 1)
end

#------------------------------------------------------------------------------

function AbstractAlgebra.degree(p::DiffPoly, v::DiffPoly)
    return degree(algdata(p), algdata(v))
end

function AbstractAlgebra.total_degree(p::DiffPoly)
    return total_degree(algdata(p))
end

function AbstractAlgebra.gcd(a::DiffPoly, b::DiffPoly)
    return DiffPoly(parent(a), gcd(algdata(a), algdata(b)))
end

function AbstractAlgebra.coeffs(a::DiffPoly)
    return coeffs(algdata(a))
end

function AbstractAlgebra.divexact(p::DiffPoly, q::DiffPoly)
    return DiffPoly(parent(p), divexact(algdata(p), algdata(q)))
end

#------------------------------------------------------------------------------

function order(p::DiffPoly, v::DiffPoly)
    R = parent(p)
    varname = decompose_derivative(var_to_str(algdata(v)))[1]
    orders = []
    for v in vars(algdata(p))
        vname, ord = decompose_derivative(var_to_str(v))
        if vname == varname
            push!(orders, ord)
        end
    end
    if length(orders) == 0
        return -1
    end
    return max(orders...)
end

function order(r::DiffRational, v::DiffPoly)
    return max(order(numerator(r), v), order(denominator(r), v))
end

#------------------------------------------------------------------------------

function leader(p::DiffPoly, v::DiffPoly)
    """
    Finds the highest derivative of variable v appearing in p
    """
    R = parent(p)
    h = order(p, v)
    if h == -1
        throw(DomainError("Variable $v does not apperar in $p"))
    end
    return get_der(R, decompose_derivative(var_to_str(algdata(v)))[1], h)
end

function leader(p::DiffPoly)
    if length(parent(p).varnames) != 1
        throw(DomainError("Expected a unvariate differential polynomial but got $p"))
    end
    return leader(p, gens(parent(p))[1])
end

#------------------------------------------------------------------------------

function initial(p::DiffPoly, v::DiffPoly)
    """
    Finds the initail of p with respect to variable v
    """
    ld = leader(p, v)
    d = degree(p, ld)
    initial = extract_coefficients(algdata(p), [algdata(ld)])[ [d] ]
    return DiffPoly(parent(p), initial)
end

function initial(p::DiffPoly)
    if length(parent(p).varnames) != 1
        throw(DomainError("Expected a unvariate differential polynomial but got $p"))
    end
    return initial(p, gens(parent(p))[1])
end


#------------------------------------------------------------------------------

function separant(p::DiffPoly, v::DiffPoly)
    """
    Finds the separant of p with respect to variable v
    """
    ld = leader(p, v)
    return DiffPoly(parent(p), derivative(algdata(p), algdata(ld)))
end

function separant(p::DiffPoly)
    if length(parent(p).varnames) != 1
        throw(DomainError("Expected a unvariate differential polynomial but got $p"))
    end
    return separant(p, gens(parent(p))[1])
end

#------------------------------------------------------------------------------

function diffred(g::DiffPoly, f::DiffPoly, v::DiffPoly)
    """
    Performs a differential reduction of g with respect to f considered as polynomials v
    Returns a triple (q, a, b), where
      - q is the remainder
      - a is the power of the separant
      - b is the power of the initial
    """
    check_parent(f, g)
    result = g
    a, b = 0, 0

    h = order(f, v)
    ld = leader(f, v)
    deg = degree(f, ld)
    sep = separant(f, v)
    init = initial(f, v)

    while order(result, v) > h
        H = order(result, v)
        ld_res = leader(result, v)
        D = degree(result, ld_res)
        result = result * sep - initial(result, v) * ld_res^(D - 1) * d(f, H - h)
        a += 1
    end

    while degree(result, ld) > deg
        D = degree(result, ld)
        result = result * init - initial(result, v) * ld^(D - deg) * f
        b += 1
    end

    return (result, a, b)
end

function diffred(g::DiffPoly, f::DiffPoly)
    check_parent(f, g)
    if length(parent(g).varnames) != 1
        throw(DomainError("Expected a unvariate differential polynomial but got $p"))
    end
    return diffred(g, f, gens(parent(g))[1])
end

#------------------------------------------------------------------------------

function evaluate_alg(a::DiffPoly, eval_point::Array{<: Any, 1})
    return evaluate(algdata(a), eval_point)
end

function evaluate_alg(a::DiffRational, eval_point::Array{<: Any, 1})
    return evaluate(numerator(algdata(a)), eval_point) // evaluate(denominator(algdata(a)), eval_point)
end

function Generic.evaluate(a::DiffRingElem, eval_point::Array{<: Any, 1})
    R = parent(a)
    eval_dict = Dict()
    for (i, v) in enumerate(gens(R))
        eval_dict[v] = eval_point[i]
        for ord in 1:order(a, v)
            prev = eval_dict[v]
            v = d(v)
            eval_dict[v] = d(prev)
        end
    end
    eval_point_alg = [get(eval_dict, R(v), zero(R)) for v in gens(R.poly_ring)]
    return evaluate_alg(a, eval_point_alg)
end

#------------------------------------------------------------------------------

# orders must be sorted in the nondecreasing order
function wronskian(polys::Array{DiffPoly, 1}, orders::Array{<: Integer, 1})
    R = parent(polys[1])
    n = length(polys)
    S = AbstractAlgebra.MatrixSpace(R, n, n)
    W = zero(S)
    for (i, p) in enumerate(polys)
        W[i, 1] = d(p, orders[1])
        for j in 2:n
            W[i, j] = d(W[i, j - 1], orders[j] - orders[j - 1])
        end
    end
    return W
end

function wronskian(polys::Array{DiffPoly, 1})
    return wronskian(polys, [i - 1 for i in 1:length(polys)])
end

function wronskian(polys::Array{DiffPoly, 1}, ind_to_omit::Integer)
    return wronskian(polys, [i for i in 0:length(polys) if i != ind_to_omit])
end
