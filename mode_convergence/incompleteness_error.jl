using HDF5
using Makie
using WGLMakie
using LsqFit

function get_path(Nm, Nc, ΩR, fname="out.h5")
    fpath  = "Nm$(Nm)_Nc$(Nc)_R$(replace(string(ΩR), '.'=>'p'))_a10_0_Em2_0_sx60/$fname"
    return joinpath(@__DIR__, fpath)
end

function plot_propagation()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"d", title = L"N_M = 5000", xticks = 0:100:500, yticks = [100, 300, 500])
    xlims!(ax1, 0, 500)
    ylims!(ax1, -0.1, 1)
    ax2 = Axis(fig[2,1], xlabel="Time (fs)", ylabel=L"d", title = L"N_M = 10000", xticks = 0:100:500, yticks = [100, 300, 500])
    xlims!(ax2, 0, 500)
    ylims!(ax2, -0.1, 1)
    ax3 = Axis(fig[3,1], xlabel="Time (fs)", ylabel=L"d", title = L"N_M = 20000", xticks = 0:100:500, yticks = [100, 300, 500])
    xlims!(ax3, 0, 500)
    ylims!(ax3, -0.1, 1)
    axs = [ax1, ax2, ax3]


    # Get time values
    (ti,δt,tf) = h5read(get_path(5000, 800, 0.1), "time_range")
    t = collect(ti:δt:tf) .* 1000

    for (i,Nm) in enumerate([5000,10000,20000])

        # Get data for Nc = 1601V
        ref = h5read(get_path(Nm, 800, 0.1), "mean_square_disp")

        for Nc in [50, 100, 200, 400, 500]
            x2 = h5read(get_path(Nm, Nc, 0.1), "mean_square_disp")

            err = abs.(x2 .- ref) ./ ref

            #d = (sqrt.(x2) .- sqrt(x2[1])) ./ 10
            lines!(axs[i], t, err, label = L"N_c = %$(2*Nc+1)")
        end
    end

    fig[1:3,2] = Legend(fig, ax1)
    fig
end

function plot_grid()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"d", title = L"N_M = 5000", xticks = 0:100:500, yticks=0:0.1:1)
    xlims!(ax1, 0, 500)
    ylims!(ax1, -0.1, 1)

    ax2 = Axis(fig[1,2], xlabel="Time (fs)", ylabel=L"d", title = L"N_M = 10000", xticks = 0:100:500, yticks=0:0.1:1)
    xlims!(ax2, 0, 500)
    ylims!(ax2, -0.1, 1)

    ax3 = Axis(fig[2,1], xlabel="Time (fs)", ylabel=L"d", title = L"N_M = 15000", xticks = 0:100:500, yticks=0:0.1:1)
    xlims!(ax3, 0, 500)
    ylims!(ax3, -0.1, 1)

    ax4 = Axis(fig[2,2], xlabel="Time (fs)", ylabel=L"d", title = L"N_M = 20000", xticks = 0:100:500, yticks=0:0.1:1)
    xlims!(ax4, 0, 500)
    ylims!(ax4, -0.1, 1)

    axs = [ax1, ax2, ax3, ax4]


    # Get time values
    (ti,δt,tf) = h5read(get_path(5000, 800, 0.1), "time_range")
    t = collect(ti:δt:tf) .* 1000

    for (i,Nm) in enumerate([5000,10000,15000,20000])

        # Get data for Nc = 1601V
        ref = h5read(get_path(Nm, 800, 0.1), "mean_square_disp")

        for Nc in [50, 100, 200, 400, 500]
            x2 = h5read(get_path(Nm, Nc, 0.1), "mean_square_disp")

            err = abs.(x2 .- ref) ./ ref

            #d = (sqrt.(x2) .- sqrt(x2[1])) ./ 10
            lines!(axs[i], t, err, label = L"N_c = %$(2*Nc+1)")
        end
    end

    #fig[1:3,2] = Legend(fig, ax1)
    fig
end

function plot_error()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    ylab = "Error w.r.t Nc=1601"
    ax1 = Axis(fig[1,1], xlabel=L"\mathrm{Number}\;\mathrm{of}\;\mathrm{cavity}\;\mathrm{modes}\;(N_c)", ylabel=ylab)
    ax2 = Axis(fig[2,1], xlabel=L"\mathrm{Cavity}\;\mathrm{energy}\;\mathrm{cutoff}\;(\mathrm{eV})", ylabel=ylab)
    
    xlims!(ax2, (2,3.5))
    ylims!(ax2, (-0.05,1))

    for Nm in [5000,8000,10000,15000,20000]

        # Get data for Nc = 1601V
        ref = h5read(get_path(Nm, 800, 0.1), "mean_square_disp")
        sumref2 = sum(ref.^2)
        err_vals = zeros(5)
        cav_e = zeros(5)
        Ncvals = [50, 100, 200, 400, 500]
        for (j,Nc) in enumerate(Ncvals)
            x2 = h5read(get_path(Nm, Nc, 0.1), "mean_square_disp")

            
            open(get_path(Nm, Nc, 0.1, "log.out"), "r") do io
                m = match(r"Maximum Cavity energy:\s+?(\d+?\.\d+)", read(io, String))
                if isnothing(m)
                    println("Couldnt match $(get_path(Nm, Nc, 0.1, "log.out"))")
                else
                    cav_e[j] = parse(Float64, m.captures[1])
                end
            end

            
            err_vals[j] = sum((x2 .- ref).^2) / sumref2

            #d = (sqrt.(x2) .- sqrt(x2[1])) ./ 10
        end
        scatter!(ax1, Ncvals, err_vals, label=L"N_M = %$(Nm)")
        lines!(ax1, Ncvals, err_vals)

        scatter!(ax2, cav_e, err_vals)#, label=L"N_M = %$(Nm)")
        #lines!(ax2, cav_e, err_vals, label=L"N_M = %$(Nm)")
        m(e,p) = exp.(p[1] * (e.-2))
        p0 = [1.0]
        fit = curve_fit(m, cav_e, err_vals, p0)
        println(fit.param)
        #lines!(ax2, cav_e, [exp(fit.param[1]*(e-2)) for e in cav_e])
    end

    lines!(ax2, 2:0.05:3.5, [exp(-10*(e-2)) for e in 2:0.05:3.5], color=:black, linestyle=:dash, label=L"\exp[-\alpha(E_\mathrm{cutoff}-E_\mathrm{min})]")

    axislegend(ax1)
    axislegend(ax2)

    fig
end
