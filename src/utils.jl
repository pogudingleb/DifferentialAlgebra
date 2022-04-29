# Function taken from the elimination project, to insert a reference
#------------------------------------------------------------------------------

function dict_to_poly(dict_monom::Dict{Array{Int, 1}, <: RingElem}, poly_ring::MPolyRing)
    builder = MPolyBuildCtx(poly_ring)
    for (monom, coef) in pairs(dict_monom)
        push_term!(builder, poly_ring.base_ring(coef), monom)
    end
    return finish(builder)
end

#------------------------------------------------------------------------------

function parent_ring_change(poly::MPolyElem, new_ring::MPolyRing)
    """
    Converts a polynomial to a different polynomial ring
    Input
      - poly - a polynomial to be converted
      - new_ring - a polynomial ring such that every variable name
          appearing in poly appears among the generators
    Output: a polynomial in new_ring "equal" to poly
    """
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any, 1}()
    for u in symbols(old_ring)
        push!(
            var_mapping,
            findfirst(v -> (string(u) == string(v)), symbols(new_ring))
        )
    end
    builder = MPolyBuildCtx(new_ring)
    one_new = one(new_ring.base_ring)
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        # detecting Singular.n_unknown
        if (typeof(coef) <: Singular.n_unknown) && !(typeof(one_new) <: Singular.n_unknown)
            coef = Singular.libSingular.julia(Singular.libSingular.cast_number_to_void(coef.ptr))
        end
        if (typeof(coef) <: Singular.FieldElemWrapper) 
            coef = coef.data
        end
        new_exp = [0 for _ in gens(new_ring)]
        for i in 1:length(exp)
            if exp[i] != 0
                if var_mapping[i] == nothing
                    throw(Base.ArgumentError("The polynomial contains a variable not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
                end
            end
        end
        push_term!(builder, new_ring.base_ring(coef), new_exp)
    end
    return finish(builder)
end

#------------------------------------------------------------------------------

function str_to_var(s::String, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind == nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

#------------------------------------------------------------------------------

function var_to_str(v::MPolyElem)
    ind = findfirst(vv -> vv == v, gens(parent(v)))
    return string(symbols(parent(v))[ind])
end

#------------------------------------------------------------------------------

function switch_ring(v::MPolyElem, ring::MPolyRing)
    """
    For a variable v, returns a variable in ring with the same name
    """
    ind = findfirst(vv -> vv == v, gens(parent(v)))
    return str_to_var(string(symbols(parent(v))[ind]), ring)
end

#------------------------------------------------------------------------------

function extract_coefficients(poly::P, variables::Array{P, 1}) where P <: MPolyElem
    """
    Intput:
        poly - multivariate polynomial
        variables - a list of variables from the generators of the ring of p
    Output:
        dictionary with keys being tuples of length len(variables) and values being 
        polynomials in the variables other than variables which are the coefficients
        at the corresponding monomials (in a smaller polynomial ring)
    """
    var_to_ind = Dict([(v, findfirst(e -> (e == v), gens(parent(poly)))) for v in variables])
    indices = [var_to_ind[v] for v in variables]

    coeff_vars = filter(v -> !(var_to_str(v) in map(var_to_str, variables)), gens(parent(poly)))
    new_ring, new_vars = AbstractAlgebra.PolynomialRing(base_ring(parent(poly)), map(var_to_str, coeff_vars))
    coeff_var_to_ind = Dict([(v, findfirst(e -> (e == v), gens(parent(poly)))) for v in coeff_vars])
    FieldType = typeof(one(base_ring(new_ring)))

    result = Dict{Array{Int, 1}, Dict{Array{Int, 1}, FieldType}}()

    for (monom, coef) in zip(exponent_vectors(poly), coeffs(poly))
        var_slice = [monom[i] for i in indices]
        if !haskey(result, var_slice)
            result[var_slice] = Dict{Array{Int, 1}, FieldType}()
        end
        new_monom = [0 for _ in 1:length(coeff_vars)]
        for i in 1:length(new_monom)
            new_monom[i] = monom[coeff_var_to_ind[coeff_vars[i]]]
        end
        result[var_slice][new_monom] = coef
    end

    return Dict(k => dict_to_poly(v, new_ring) for (k, v) in result)
end


