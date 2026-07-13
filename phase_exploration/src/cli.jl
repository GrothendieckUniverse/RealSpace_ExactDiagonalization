"Minimal dependency-free parser for `--key value` and boolean `--flag` options."
function parse_cli(args::Vector{String})
    opts = Dict{String,String}()
    positional = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--")
            key = arg[3:end]
            isempty(key) && error("Empty CLI option.")
            if occursin('=', key)
                k, v = split(key, '='; limit=2)
                opts[k] = v
            elseif i < length(args) && !startswith(args[i + 1], "--")
                opts[key] = args[i + 1]
                i += 1
            else
                opts[key] = "true"
            end
        else
            push!(positional, arg)
        end
        i += 1
    end
    return opts, positional
end

function parse_geometry(text::AbstractString)
    normalized = replace(strip(text), '[' => "", ']' => "", ',' => "x", 'X' => "x")
    fields = split(normalized, 'x')
    length(fields) == 2 || error("Geometry must look like 3x5 or [3,5]; got `$text`.")
    return (parse(Int, strip(fields[1])), parse(Int, strip(fields[2])))
end

function parse_bool(text::AbstractString)
    value = lowercase(strip(text))
    value in ("true", "yes", "1", "on") && return true
    value in ("false", "no", "0", "off") && return false
    error("Expected boolean, got `$text`.")
end

parse_float_list(text::AbstractString) = Float64[parse(Float64, strip(x)) for x in split(text, ',') if !isempty(strip(x))]

function parse_mode(text::AbstractString)
    mode = Symbol(lowercase(strip(text)))
    mode in (:matrix, :matrixfree) || error("Mode must be matrix or matrixfree; got `$text`.")
    return mode
end
