# From the repository root:
# julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
# julia --project=docs docs/generate_readme_figure.jl
using KomaMRI

root = dirname(@__DIR__)
readme = read(joinpath(root, "README.md"), String)
example = first(eachmatch(r"```julia\n(.*?)```"s, readme)).captures[1]

sequence = mktempdir() do directory
    cd(directory) do
        include_string(Main, example * "\nsequence", "README.md")
    end
end

assets = joinpath(@__DIR__, "assets")
mkpath(assets)
sequence_plot = plot_seq(sequence)
sequence_plot.layout[:updatemenus] = []
for trace in sequence_plot.data
    trace[:mode] = "lines"
end
kspace_plot = plot_kspace(sequence; view_2d=true)
kspace_range = [-150, 150]
kspace_plot.layout[:xaxis][:range] = kspace_range
kspace_plot.layout[:yaxis][:range] = kspace_range
kspace_plot.layout[:xaxis][:constrain] = "domain"
kspace_plot.layout[:yaxis][:constrain] = "domain"
savefig(sequence_plot, joinpath(assets, "spin_echo_sequence.png"); width=1200, height=420, scale=2)
savefig(kspace_plot, joinpath(assets, "spin_echo_kspace.png"); width=700, height=500, scale=2)
