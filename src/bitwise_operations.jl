module BitWise_Operations

using MLStyle

const Mask = UInt # bit representation of the state in integer for state indexing
const COMPLEX_ONE = ComplexF64(1)

"bit representation for the occupied/empty two-state vertex. Here `i` is the _1-based_ linear index of the vertex in the graph"
@inline bitmask_of_site(i::Int)::Mask = @fastmath Mask(1) << (i - 1)
"occupy vertex-i for mask `m`"
@inline occupy_site_for_mask(m::Mask, i::Int)::Mask = @fastmath m | bitmask_of_site(i)
"empty vertex-i for mask `m`"
@inline empty_site_for_mask(m::Mask, i::Int)::Mask = @fastmath m & ~bitmask_of_site(i)
"is vertex-i occupied in mask `m`?"
@inline is_site_occupied(m::Mask, i::Int)::Bool = @fastmath (m & bitmask_of_site(i)) != 0
"is vertex-i empty in mask `m`?"
@inline is_site_empty(m::Mask, i::Int)::Bool = @fastmath (m & bitmask_of_site(i)) == 0
"number of occupied vertices in mask `m`"
@inline n_occupied_for_mask(m::Mask) = Base.count_ones(m)


"""
Helper function to get an iterator of occupied site indices from a bit mask `m`
---
without explicit construction of the configuration vector to avoid allocation overhead.
- Args:
    - `m::Mask`: the input occupation basis state in bitmask form, where the i-th bit (1-based) indicates whether vertex-i is occupied or not
    - `n_site::Int`: the total number of sites in the lattice, which determines how many bits we need to check in the mask
- Returns:
    - `Base.Generator`: a generator that yields the indices of the occupied sites in the mask `m`. For example, if `m = 0b10110` (in binary) and `n_site = 5`, then the generator will yield `2, 3, 5` since the 2nd, 3rd, and 5th bits are set in `m`.
"""
function filled_site_iter_for_mask(m::Mask, n_site::Int)::Base.Generator
    return (i for i in 1:n_site if is_site_occupied(m, i))
end


"""
Helper function to get an iterator of empty site indices from a bit mask `m`
---
without explicit construction of the configuration vector to avoid allocation overhead.
- Args:
    - `m::Mask`: the input occupation basis state in bitmask form, where the i-th bit (1-based) indicates whether vertex-i is occupied or not
    - `n_site::Int`: the total number of sites in the lattice, which determines how many bits we need to check in the mask
- Returns:
    - `Base.Generator`: a generator that yields the indices of the empty sites in the mask `m`. For example, if `m = 0b10110` (in binary) and `n_site = 5`, then the generator will yield `1, 4` since the 1st and 4th bits are not set in `m`.
"""
function empty_site_iter_for_mask(m::Mask, n_site::Int)::Base.Generator
    return (i for i in 1:n_site if is_site_empty(m, i))
end





end # module `BitWise_Operations`