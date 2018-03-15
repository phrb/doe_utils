if length(workers()) > 1
    rmprocs(workers())
end

addprocs()

import JuliaDB

@everywhere begin
    using JuliaDB

    mutable struct Measurement
        f1_value::Int64
        f2_value::Float64
        f3_value::Bool
        f4_value::String
        response::Float64
        complete::Bool
        id::UInt64

        Measurement(f1, f2, f3, f4, r) = new(f1, f2, f3, f4, r, true, hash("$f1$f2$f3$f4"))
        Measurement(f1, f2, f3, f4) = new(f1, f2, f3, f4, -1, false, hash("$f1$f2$f3$f4"))
    end

    function random_measurement()
        Measurement(rand(1:512), randexp(), bitrand()[1], rand(["--merge", "--no-merge"]), randexp())
    end

    function init_dummy_table()
        n = rand(4:10)

        parameters = Array{Measurement, 1}(n)

        # Generating meaningless, random parameters
        for i = 1:n
            parameters[i] = random_measurement()
        end

        names = fieldnames(Measurement)

        return table([[getfield(p, n) for p in parameters] for n in names]...,
                    names = names, pkey = :id)
    end

    function contains(t::NextTable, entry::Measurement)
        return entry.id in select(t, :id)
    end

    function add_complete_unique_entry!(t::NextTable, entry::Measurement)
        if entry.complete && !contains(t, entry)
            row_expr::Expr = :@NT

            for name in fieldnames(Measurement)
                push!(row_expr.args, Expr(Symbol("="), name, getfield(entry, name)))
            end

            push!(rows(t), eval(row_expr))

            # Sorting 't' by its primary key
            t = table(t, pkey = t.pkey, copy = false)
        end

        t
    end

    function measure(channel::RemoteChannel)
        put!(channel, init_dummy_table())
    end
end

using FileIO

results = RemoteChannel(() -> Channel(32))

for p in workers()
    @async remote_do(measure, p, results)
end

t = take!(results)

for i = 1:length(workers()) - 1
    t = merge(t, take!(results))
end

t

FileIO.save("./test.csv", t)
