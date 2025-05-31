module CairoMakieExt

using TransmissionChannelAnalysis
using CairoMakie: scatterlines!, barplot!, xlims!, PolyElement, LineElement,
    MarkerElement, Legend, Figure, Axis
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

    ax = plot_decomposition!(idx_outcome, irfs, teffects, channel_names; colors=colors)

    if legend
        elements = [PolyElement(polycolor=colors[i]) for i = 1:length(teffects)]
        elements = vcat([[LineElement(color=:black), MarkerElement(marker=:circle, color=:black)]], elements)
        Legend(
            fig[2, :],
            elements,
            vcat(["Total"], channel_names),
            "Effect";
            orientation=:horizontal,
            framevisible=false
        )
    end

    return TCAPlot(fig, ax)
end

"""
    plot_decomposition!(
        ax::Axis,
        idx_outcome::Int,
        irfs::AbstractArray{<:Number,3},
        teffects::Vector{<:AbstractArray{<:Number,3}},
        channel_names::AbstractVector{<:String};
        colors=wong_colors()
    ) -> TCAPlot

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
- `channel_names::AbstractVector{<:String}`: Names of the transmission
  channels.

# Keyword Arguments

- `colors`: Colors for the decomposition bars (default: Wong colors).

"""
function TransmissionChannelAnalysis.plot_decomposition!(
    ax::Axis,
    idx_outcome::Int,
    irfs::AbstractArray{<:Number,3},
    teffects::Vector{<:AbstractArray{<:Number,3}},
    channel_names::AbstractVector{<:String};
    colors=wong_colors()
)

    max_horizon = size(irfs, 3) - 1
    tbl = (;
        x=horizons,
        total=total[idx_outcome, 1, :],
        x_bar=repeat(horizons, length(teffects)),
        decomposition=reduce(vcat, te[idx_outcome, 1, :] for te in teffects),
        grp=repeat(1:length(teffects); inner=length(horizons))
    )

    barplot!(ax, tbl.x_bar, tbl.decomposition; stack=tbl.grp, color=colors[tbl.grp])
    scatterlines!(ax, tbl.x, tbl.total; color=:black)
    xlims!(ax, -1, max_horizon + 1)

    return ax
end
end
