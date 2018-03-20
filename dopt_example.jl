ranges = [4:4, 5:5, 6:6, 12:12, 24:24, 48:48, 4:3 ^ 4]
#ranges = [4:4, 5:5, 6:6]
factors = [[Array{Float64, 1}(1:3) for i = 1:4] for k in 1:length(ranges)]
designs = 8000

sampled_subsets = sample_subsets(factors, ranges, designs, scale = scale_box_encoding!)

plot_subsets(sampled_subsets)
