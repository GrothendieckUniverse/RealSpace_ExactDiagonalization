function ensure_parent(path::AbstractString)
    mkpath(dirname(path))
    return path
end

function write_key_values(path::AbstractString, pairs)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "key,value")
        for (key, value) in pairs
            println(io, key, ',', value)
        end
    end
    return path
end

function read_simple_csv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return NamedTuple[]
    header = Symbol.(split(lines[1], ','))
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = split(line, ','; keepempty=true)
        length(values) == length(header) || error("Malformed CSV row in $path: $line")
        push!(rows, NamedTuple{Tuple(header)}(Tuple(values)))
    end
    return rows
end

csv_float(value) = parse(Float64, value)
csv_int(value) = parse(Int, value)

function spectrum_table(ed_data)
    rows = NamedTuple[]
    for (irrep_idx, (values, _)) in ed_data.ed_scan_res
        sector = Tuple(Int.(ed_data.irrep_list[irrep_idx].label))
        for (level, energy) in enumerate(values)
            push!(rows, (energy=Float64(energy), sector=sector,
                irrep_idx=irrep_idx, level=level))
        end
    end
    isempty(rows) && error("ED result contains no eigenvalues.")
    sort!(rows; by=row -> row.energy)
    e0 = rows[1].energy
    return [(; row..., shifted=row.energy - e0) for row in rows]
end

function lowest_unique_sectors(table; count::Int=3)
    sectors = Tuple{Int,Int}[]
    for row in table
        row.sector in sectors && continue
        push!(sectors, row.sector)
        length(sectors) == count && break
    end
    return sectors
end

"Return the globally lowest `count` eigenstates, retaining both sector and in-sector level."
function lowest_manifold_states(table; count::Int=3)
    count > 0 || error("Manifold size must be positive; got $count.")
    length(table) >= count || error(
        "Requested a $count-state manifold, but the spectrum contains only $(length(table)) states.")
    return table[1:count]
end

function all_sector_labels(ed_data)
    return [Tuple(Int.(irrep.label)) for irrep in ed_data.irrep_list]
end

function write_spectrum_csv(path, table, sample, x)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "L1,L2,tpp_numerator,tpp_actual,rank,energy,energy_minus_E0,k1,k2,irrep_index,level_in_sector")
        for (rank, row) in enumerate(table)
            @printf(io, "%d,%d,%.16g,%.16g,%d,%.16g,%.16g,%d,%d,%d,%d\n",
                sample[1], sample[2], x, actual_tpp(x), rank, row.energy, row.shifted,
                row.sector[1], row.sector[2], row.irrep_idx, row.level)
        end
    end
    return path
end

function write_allowed_sf_csv(path, qx, qy, values)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "qx,qy,S_q")
        for i in eachindex(values)
            @printf(io, "%.16g,%.16g,%.16g\n", qx[i], qy[i], values[i])
        end
    end
    return path
end

function write_dense_sf_csv(path, kx, ky, values)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "kx,ky,S_q")
        for iy in eachindex(ky), ix in eachindex(kx)
            @printf(io, "%.16g,%.16g,%.16g\n", kx[ix], ky[iy], values[ix, iy])
        end
    end
    return path
end

function sf_metrics(qx, qy, values)
    abs_values = abs.(values)
    peak_idx = argmax(abs_values)
    mean_s = mean(values)
    mean_abs_s = mean(abs_values)
    requested_ratio = abs(mean_s) <= eps(Float64) ? NaN : abs_values[peak_idx] / abs(mean_s)
    robust_ratio = mean_abs_s <= eps(Float64) ? NaN : abs_values[peak_idx] / mean_abs_s
    return (
        peak_idx=peak_idx,
        peak_q=(qx[peak_idx], qy[peak_idx]),
        max_abs_S=abs_values[peak_idx],
        mean_S=mean_s,
        mean_abs_S=mean_abs_s,
        max_abs_S_over_mean_S=requested_ratio,
        max_abs_S_over_mean_abs_S=robust_ratio,
    )
end

function write_sf_metrics_csv(path, metrics, sample, x, ground_sector)
    ensure_parent(path)
    open(path, "w") do io
        println(io, "L1,L2,n_sites,n_particles,tpp_numerator,tpp_actual,ground_k1,ground_k2,peak_qx,peak_qy,max_abs_S,mean_S,mean_abs_S,max_abs_S_over_mean_S,max_abs_S_over_mean_abs_S")
        @printf(io, "%d,%d,%d,%d,%.16g,%.16g,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n",
            sample[1], sample[2], 2 * prod(sample), default_particle_number(sample),
            x, actual_tpp(x), ground_sector[1], ground_sector[2],
            metrics.peak_q[1], metrics.peak_q[2], metrics.max_abs_S,
            metrics.mean_S, metrics.mean_abs_S, metrics.max_abs_S_over_mean_S,
            metrics.max_abs_S_over_mean_abs_S)
    end
    return path
end

function result_point_dir(kind::AbstractString, sample::Tuple{Int,Int}, x::Real)
    return joinpath(RESULT_ROOT, kind, geometry_tag(sample), "x_$(tpp_tag(x))")
end
