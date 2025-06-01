module CairoMakieExt

using TransmissionChannelAnalysis
using CairoMakie: scatterlines!, barplot!, xlims!, PolyElement, LineElement,
    MarkerElement, Legend, Figure, Axis, GridPosition
using Makie: wong_colors

"""
    TCAPlot

Convenience struct for CairoMakie plotting functions in the TCA package.
Stores a `Figure` and a `Axis`, and implements custom behavior
for display and unpacking.

## Fields

- `fig::Figure`
- `ax::Axis`

## Notes

The `TCAPlot` struct allows:

- Displaying the figure when called directly in the REPL or script:
  ```julia
  tcaplot  # displays the figure
  tcaplot; # does not display figure (semicolon suppression)
  ```

- Flexible unpacking of components:
  ```julia
  fig, = tcaplot       # returns only the figure
  fig, ax = tcaplot    # returns figure and axis
  ```

The figure is displayed automatically in interactive contexts unless the
call is suppressed by a semicolon.
"""
struct TCAPlot
    fig::Figure
    ax::Axis
end
function Base.iterate(tcaplot::TCAPlot, state::Int=1)
    state == 1 && return (tcaplot.fig, 2)
    state == 2 && return (tcaplot.ax, 3)
    return nothing
end
Base.show(::IO, tcaplot::TCAPlot) = display(tcaplot.fig)


"""
    plot_decomposition(
        idx_outcome::Int,
        irfs::AbstractArray{<:Number,3},
        teffects::Vector{<:AbstractArray{<:Number,3}},
        channel_names::AbstractVector{<:String};
        title::String="",
        legend::Bool=false,
        colors=wong_colors()
    ) -> TCAPlot

Creates a decomposition plot of the total effects (`irfs`) and the
transmission effects (`teffects`) for a given outcome variable (idx_outcome).

The function returns a `TCAPlot` object, which displays the figure
immediately in interactive environments (unless suppressed with a
semicolon). You can also extract the separate `Figure` and `Axis`:

```julia
fig, = plot_decomposition(...)
fig, ax = plot_decomposition(...)
```

# Arguments

- `idx_outcome::Int`: Index of the outcome variable.
- `irfs::AbstractArray{<:Number,3}`: Total effects, with shape `[N, 1, H]`
  where `N` is the number of endogenous variables, and `H` is the number of
  horizons (including zero).
- `teffects::Vector{<:AbstractArray{<:Number,3}}`: Transmission effects for
  each transmission channel, same shape as `irfs`.
- `channel_names::AbstractVector{<:String}`: Names of the transmission
  channels.

# Keyword Arguments

- `title::String`: Optional figure title.
- `legend::Bool`: Whether to plot a legend.
- `colors`: Colors for the decomposition bars (default: Wong colors).

# Returns

- A `TCAPlot` that includes the figure and axis for further styling.
"""
function TransmissionChannelAnalysis.plot_decomposition(
    idx_outcome::Int,
    irfs::AbstractArray{<:Number,3},
    teffects::Vector{<:AbstractArray{<:Number,3}},
    channel_names::AbstractVector{<:String};
    title::String="",
    legend::Bool=false,
    colors=wong_colors()
)

    fig = Figure()
    ax = Axis(fig[1, 1]; title=title)

    ax = plot_decomposition!(ax, idx_outcome, irfs, teffects; colors=colors)

    if legend
        add_decomposition_legend!(fig[2, :], channel_names; colors=colors)
    end

    return TCAPlot(fig, ax)
end

"""
    plot_decomposition!(
        ax::Axis,
        idx_outcome::Int,
        irfs::AbstractArray{<:Number,3},
        teffects::Vector{<:AbstractArray{<:Number,3}};
        colors=wong_colors()
    ) -> Axis

Plots a decomposition plot of the total effects (`irfs`) and the
transmission effects (`teffects`) for a given outcome variable (idx_outcome)
into the axis `ax`.

# Arguments

- `idx_outcome::Int`: Index of the outcome variable.
- `irfs::AbstractArray{<:Number,3}`: Total effects, with shape `[N, 1, H]`
  where `N` is the number of endogenous variables, and `H` is the number of
  horizons (including zero).
- `teffects::Vector{<:AbstractArray{<:Number,3}}`: Transmission effects for
  each transmission channel, same shape as `irfs`.

# Keyword Arguments

- `colors`: Colors for the decomposition bars (default: Wong colors).

"""
function TransmissionChannelAnalysis.plot_decomposition!(
    ax::Axis,
    idx_outcome::Int,
    irfs::AbstractArray{<:Number,3},
    teffects::Vector{<:AbstractArray{<:Number,3}};
    colors=wong_colors()
)

    max_horizon = size(irfs, 3) - 1
    horizons = 0:max_horizon
    tbl = (;
        x=horizons,
        total=irfs[idx_outcome, 1, :],
        x_bar=repeat(horizons, length(teffects)),
        decomposition=reduce(vcat, te[idx_outcome, 1, :] for te in teffects),
        grp=repeat(1:length(teffects); inner=length(horizons))
    )

    barplot!(ax, tbl.x_bar, tbl.decomposition; stack=tbl.grp, color=colors[tbl.grp])
    scatterlines!(ax, tbl.x, tbl.total; color=:black)
    xlims!(ax, -1, max_horizon + 1)

    return ax
end

"""
    add_decomposition_legend!(
        gp::GridPosition, 
        channel_names::AbstractVector{<:String}; 
        colors=wong_colors()
    ) -> Legend

Add a legend to a decomposition plot. 

The decomposition plot can be created using `plot_decomposition` or 
`plot_decomposition!`. 

# Arguments

- `gp::GridPosition`: Position to plot the legend into. 
- `channel_names::AbstractVector{<:String}`: Names of the transmission
  channels.

# Keyword Arguments

- `colors`: Colors for the decomposition bars (default: Wong colors).

# Returns 

- A `Legend` object that can be further manipulated. 
"""
function TransmissionChannelAnalysis.add_decomposition_legend!(
    gp::GridPosition,
    channel_names::AbstractVector{<:String};
    colors=wong_colors()
)

    elements = [PolyElement(polycolor=colors[i]) for i = 1:length(channel_names)]
    elements = vcat([[LineElement(color=:black), MarkerElement(marker=:circle, color=:black)]], elements)
    lgd = Legend(
        gp,
        elements,
        vcat(["Total"], channel_names),
        "Effect";
        orientation=:horizontal,
        framevisible=false
    )

    return lgd
end

"""
    plot_decomposition_comparison(
        idx_outcome::Int,
        irfs::AbstractArray{<:Number,3},
        channel_names::AbstractVector{<:String},
        decomposition_names::AbstractVector{<:String},
        teffects::Vector{<:AbstractArray{<:Number,3}}...;
        title::String="",
        legend::Bool=true,
        colors=wong_colors()
    ) -> TCAPlot

Creates a decomposition comparison plot of the total effects (`irfs`) and the
transmission effects (`teffects`), of various transmission matrices, for a given 
outcome variable (idx_outcome).

The function returns a `TCAPlot` object, which displays the figure
immediately in interactive environments (unless suppressed with a
semicolon). You can also extract the separate `Figure` and `Axis`:

```julia
fig, = plot_decomposition_comparison(...)
fig, ax = plot_decomposition_comparison(...)
```

# Arguments

- `idx_outcome::Int`: Index of the outcome variable.
- `irfs::AbstractArray{<:Number,3}`: Total effects, with shape `[N, 1, H]`
  where `N` is the number of endogenous variables, and `H` is the number of
  horizons (including zero).
- `channel_names::AbstractVector{<:String}`: Names of the transmission
- `decomposition_names::AbstractVector{<:String}`: Names of the decompositions
- `teffects::Vector{<:AbstractArray{<:Number,3}}...`: Transmission effects for
  each transmission channel and each transmission matrix; Each transmission effect 
  must have the same shape as `irfs`.

# Keyword Arguments

- `title::String`: Optional figure title.
- `legend::Bool`: Whether to plot a legend.
- `colors`: Colors for the decomposition bars (default: Wong colors).

# Returns

- A `TCAPlot` that includes the figure and axis for further styling.
"""
function TransmissionChannelAnalysis.plot_decomposition_comparison(
    idx_outcome::Int,
    irfs::AbstractArray{<:Number,3},
    channel_names::AbstractVector{<:String},
    decomposition_names::AbstractVector{<:String},
    teffects::Vector{<:AbstractArray{<:Number,3}}...;
    title::String="",
    legend::Bool=true,
    colors=wong_colors()
)
    fig = Figure()
    ax = Axis(fig[1, 1]; title=title)

    ax = plot_decomposition_comparison!(ax, idx_outcome, irfs, teffects...; colors=colors)

    if legend
        add_decompare_legend!(fig[1, 2], channel_names, decomposition_names; colors=colors)
    end

    return TCAPlot(fig, ax)
end

"""
    plot_decomposition_comparison!(
        ax::Axis,
        idx_outcome::Int,
        irfs::AbstractArray{<:Number,3},
        teffects::Vector{<:AbstractArray{<:Number,3}}...;
        colors=wong_colors()
    ) -> Axis

Plots a decomposition comparison plot of the total effects (`irfs`) and the
transmission effects (`teffects`), of various transmission matrices, for a given 
outcome variable (idx_outcome) into the axis `ax`.

# Arguments

- `idx_outcome::Int`: Index of the outcome variable.
- `irfs::AbstractArray{<:Number,3}`: Total effects, with shape `[N, 1, H]`
  where `N` is the number of endogenous variables, and `H` is the number of
  horizons (including zero).
- `teffects::Vector{<:AbstractArray{<:Number,3}}...`: Transmission effects for
  each transmission channel and each transmission matrix; Each transmission effect 
  must have the same shape as `irfs`.

# Keyword Arguments

- `colors`: Colors for the decomposition bars (default: Wong colors).

"""
function TransmissionChannelAnalysis.plot_decomposition_comparison!(
    ax::Axis,
    idx_outcome::Int,
    irfs::AbstractArray{<:Number,3},
    teffects::Vector{<:AbstractArray{<:Number,3}}...;
    colors=wong_colors()
)

    n_decomps = length(teffects)
    n_channels = length(teffects[1])

    max_horizon = size(irfs, 3) - 1
    horizons = 0:max_horizon
    tbl = (;
        x=horizons,
        total=irfs[idx_outcome, 1, :],
        x_bar=repeat(horizons, n_channels * n_decomps),
        y_bar=reduce(vcat, reduce(vcat, te[idx_outcome, 1, :] for te in tes) for tes in teffects),
        grp_stack=repeat(1:n_channels; inner=length(horizons), outer=n_decomps),
        grp_dodge=repeat(1:n_decomps; inner=length(horizons)*n_channels),
        grp=repeat(1:(n_channels*n_decomps); inner=length(horizons))
    )

    barplot!(ax, tbl.x_bar, tbl.y_bar; stack=tbl.grp_stack, dodge=tbl.grp_dodge, color=colors[tbl.grp], dodge_gap=0.12)
    scatterlines!(ax, tbl.x, tbl.total; color=:black)
    xlims!(ax, -1, max_horizon + 1)

    return ax
end

"""
    add_decompare_legend!(
        gp::GridPosition, 
        channel_names::AbstractVector{<:String}, 
        decomposition_names::AbstractVector{<:String}; 
        colors=wong_colors()
    ) -> Legend

Add a legend to a decomposition comparison plot. 

The decomposition comparison plot can be created using 
`plot_decomposition_comparison` or `plot_decomposition_comparison!`. 

# Arguments

- `gp::GridPosition`: Position to plot the legend into. 
- `channel_names::AbstractVector{<:String}`: Names of the transmission
  channels.
- `decomposition_names::AbstractVector{<:String}`: Names for the decompositions.

# Keyword Arguments

- `colors`: Colors for the decomposition bars (default: Wong colors).

# Returns 

- A `Legend` object that can be further manipulated. 
"""
function TransmissionChannelAnalysis.add_decompare_legend!(
    gp::GridPosition, 
    channel_names::AbstractVector{String}, 
    decomposition_names::AbstractVector{String}; 
    colors=wong_colors()
)

    n_channels = length(channel_names)
    n_decomps = length(decomposition_names)

    elements = [[[LineElement(color=:black), MarkerElement(marker=:circle, color=:black)]]]
    for i = 1:n_decomps
        elements = vcat(elements, [[PolyElement(polycolor=colors[j]) for j = ((i-1)*n_channels+1):(i*n_channels)]])
    end

    lgd = Legend(
        gp[1,1],     
        elements, 
        vcat([["Total"]], [channel_names for _ in 1:n_decomps]), 
        vcat([""], decomposition_names); 
        framevisible=false, 
        orientation=:vertical, 
        nbanks=1, 
        titleposition=:top, 
    )

    return lgd
end

end
