if length(workers()) > 1
    rmprocs.(workers())
end

addprocs(2)

using IterTools, DataFrames, StatsBase, StatPlots

function build_linear_formula(factors::Int)
    linear_formula = Expr(:call)
    linear_formula.args = vcat(:+, [Symbol("x", i) for i = 1:factors])
    return Formula(:Y, linear_formula)
end

function get_model_variables(formula::DataFrames.Formula)
    variables = Array{Any, 1}()
    push!(variables, DataFrames.Terms(formula).terms...)
    return variables
end

function scale_orthogonal!(design::Array{Float64, 2},
                           factors::Array{T, 1}) where T<: Any
    design_range = (max(design[:,]...) + min(design[:,]...)) / 2

    for i = 1:size(design, 2)
        factor_range = (max(factors[i]...) - min(factors[i]...)) / 2
        design[:, i] = (design[:, i] - design_range) / factor_range
    end

    return design
end

function generate_model_matrix(formula::DataFrames.Formula,
                               design::Array{Float64, 2},
                               factors::Array{T, 1}) where T <: Any
    variables  = get_model_variables(formula)

    design     = DataFrame(scale_orthogonal!(design, factors))
    new_design = DataFrame(I = ones(size(design, 1)))

    for variable in variables
        if typeof(variable) == Expr && variable.args[1] == :&
            interaction             = Symbol(variable.args[2:end]...)
            new_design[interaction] = ones(size(design, 1))

            for s in variable.args[2:end]
                new_design[interaction] .*= design[s]
            end
        else
            new_design[variable] = float(design[variable])
        end
    end

    return new_design
end


function get_prediction_variances(model_matrix::Array{Float64, 2})
    information_matrix = model_matrix' * model_matrix

    if det(information_matrix) != 0.0
        dispersion_matrix = inv(information_matrix)
        rows              = size(dispersion_matrix, 1)

        prediction_variances = [dispersion_matrix[i, :]' * dispersion_matrix *
                                dispersion_matrix[i, :] for i = 1:rows]

        return prediction_variances
    else
        return 0.0
    end
end

function d_optimality(model_matrix::Array{Float64, 2})
    det_information_matrix = det(model_matrix' * model_matrix)

    if det_information_matrix < 0.0
        return 0.0
    else
        return ^(det_information_matrix, 1 / size(model_matrix, 2))
    end
end

function a_optimality(model_matrix::Array{Float64, 2})
    information_matrix = model_matrix' * model_matrix

    if det(information_matrix) != 0.0
        return trace(inv(information_matrix)) / size(model_matrix, 2)
    else
        return 0.0
    end

end

function v_optimality(model_matrix::Array{Float64, 2})
    prediction_variances = get_prediction_variances(model_matrix)
    rows                 = size(model_matrix, 1)

    return sum(prediction_variances) / rows
end

function g_optimality(model_matrix::Array{Float64, 2})
    prediction_variances = get_prediction_variances(model_matrix)
    return max(prediction_variances...)
end

function g_efficiency(model_matrix::Array{Float64, 2})
    prediction_variances = get_prediction_variances(model_matrix)
    max_variance         = max(prediction_variances...)

    g_e = size(model_matrix, 2) / max_variance

    if g_e == Inf
        return 0.0
    else
        return g_e
    end
end

function d_efficiency_lower_bound(model_matrix::Array{Float64, 2})
    g_e   = g_efficiency(model_matrix)
    d_elb = exp(1 - (1 / g_e))

    if d_elb == Inf
        return 0.0
    else
        return d_elb
    end
end

function condition_number(model_matrix::Array{Float64, 2})
    condition_number = cond(model_matrix)

    if condition_number == Inf
        return 0.0
    else
        return condition_number
    end
end

@everywhere begin
    function sample_full_factorial(factors::Array{T, 1}) where T <: Any
        return Array{Float64, 1}([rand(i) for i in factors])
    end

    function check_repeated_row(subset::SharedArray{Float64, 2}, row::Array{Float64, 1})
        for j = 1:size(subset, 1)
            if subset[j, :] == row
                return true
            end
        end

        return false
    end

    function full_factorial_subset(factors::Array{T, 1}, experiments::Int) where T <: Any
        subset = fill!(SharedArray{Float64, 2}(experiments, size(factors, 1)), 0.0)

        @sync @parallel for i = 1:experiments
            sample_row = sample_full_factorial(factors)

            while check_repeated_row(subset, sample_row)
                sample_row = sample_full_factorial(factors)
            end

            subset[i, :] = sample_row
        end

        return subset
    end
end

function generate_designs(factors::Array{T, 1},
                          formula::DataFrames.Formula,
                          sample_range::UnitRange{Int},
                          designs::Int;
                          check_bounds::Bool = true) where T <: Any
    println("> Factors: ", factors)

    full_factorial_size = prod(length, factors)
    full_factorial_subsets = 2.0 ^ full_factorial_size

    println("> Full Factorial Size: ", full_factorial_size)
    println("> Total Subsets: ", full_factorial_subsets)
    println("> Range of Design Sizes: ", sample_range)
    println("> Number of Design to Sample: ", designs)

    if check_bounds
        if designs > full_factorial_subsets
            println("> Requested too many designs, using ",
                    full_factorial_subsets, " instead")

            designs = full_factorial_subsets
        end

        if sample_range.stop > full_factorial_size
            println("> Requested too many maximum experiments, using ",
                    full_factorial_size, " instead")
            sample_range = sample_range.start:full_factorial_size
        end

        if sample_range.start == sample_range.stop
            println("> Total Subsets for Fixed Size ", sample_range.start, ": ",
                    factorial(float(full_factorial_size)) /
                    (factorial(float(full_factorial_size - sample_range.start)) *
                    factorial(float(sample_range.start))))
        end
    else
        println("> WARNING: Skipping bounds check!")
    end

    evaluation = DataFrame(Length = [],
                           D      = [],
                           log10D = [],
                           DELB   = [],
                           A      = [],
                           V      = [],
                           G      = [],
                           CN     = [],
                           log2CN = [],
                           GE     = [])


    for i in 1:designs
        samples      = rand(sample_range)
        subset       = full_factorial_subset(factors, samples)
        model_matrix = generate_model_matrix(formula, Array{Float64, 2}(subset), factors)
        candidate    = Array(model_matrix)

        d_opt = d_optimality(candidate)
        c_n   = condition_number(candidate)

        push!(evaluation, [size(candidate, 1),
                           d_opt,
                           log(10, abs(d_opt)),
                           d_efficiency_lower_bound(candidate),
                           a_optimality(candidate),
                           v_optimality(candidate),
                           g_optimality(candidate),
                           c_n,
                           log(2, abs(c_n)),
                           g_efficiency(candidate)])
    end

    return evaluation
end

function sample_subset(factors, sample_range, designs; check_bounds = true)
    formula = build_linear_formula(length(factors))
    #formula = @formula(y ~ x1 + x2 + x3)

    run_time = @elapsed sampling_subset = generate_designs(factors,
                                                           formula,
                                                           sample_range,
                                                           designs,
                                                           check_bounds = check_bounds)
    println("> Elapsed Time: ", run_time, " seconds")

    sort!(sampling_subset, cols = :D, rev = true)

    return sampling_subset
end

function sample_subsets(factors::Array,
                        ranges::Array{UnitRange{Int}, 1},
                        designs::Int;
                        check_bounds::Bool = true)
    sampled_subsets = []

    for subset = 1:length(ranges)
        label = " "

        if ranges[subset].start == ranges[subset].stop
            label = string(ranges[subset].start, " Experiments")
        else
            label = string(ranges[subset].start, " to ",
                          ranges[subset].stop, " Experiments")
        end

        label = string(label, ", ", length(factors[subset]), " Factors")

        sampled_subset = sample_subset(factors[subset],
                                       ranges[subset],
                                       designs,
                                       check_bounds = check_bounds)

        push!(sampled_subsets,
              (sampled_subset,
              max(sampled_subset[:Length]..., sampled_subset[:D]...),
              label))
    end

    return sampled_subsets
end

function plot_subsets(sampled_subsets)
    upscale = 2
    small_font = Plots.font("sans-serif", 10.0 * upscale)
    large_font = Plots.font("sans-serif", 14.0 * upscale)
    default(titlefont  = large_font,
            guidefont  = large_font,
            tickfont   = small_font,
            legendfont = small_font)
    default(size = (896 * upscale, 504 * upscale))
    default(dpi = 300)

    plotly()

    subplots = []

    for subset in sampled_subsets
        push!(subplots,
              histogram(Array(subset[1][:D]), labels = "D",
                        title = subset[3]),
              plot(Array(subset[1][[:D, :DELB, :Length]]),
                  ylims = (-1, subset[2] + 2),
                  labels = ["D" "DELB" "Length"],
                  title = subset[3],
                  linestyle = :solid,
                  linealpha=1.0,
                  linewidth=1.5 * upscale)
        )
    end

    plot(subplots...,
        layout = (length(sampled_subsets), 2)
    )
end
