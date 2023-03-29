using Makie
using WGLMakie
using HDF5
using DataFrames
using GLM
using LaTeXStrings


function propagation(;ΩR=0.1, σx=120, d2=false, showstd=false, σMvals=[0.005, 0.01, 0.02, 0.04])

    ΩRstr = replace(string(ΩR), '.'=>'p')

	fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)
	
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"d", ylabelsize=30, xticks=0:500:2500)
    xlims!(0,2500)

	Rstr = replace(string(ΩR), "." => "p") 

    for σM in σMvals
	    sm = replace(string(σM), "." => "p") 

        d = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "$(Int(σx))_avg_d") .- σx/10
	    std_d = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "$(Int(σx))_std_d")
	    tvals = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "d_tvals") .* 1000

	    if d2
	    	d = 50 .* d .^2 ./ tvals
	    	lines!(ax, tvals, d, linewidth=3.0)
	    	ax.ylabel=L"\left\langle x \right\rangle^2 /2t \;\mathrm{(nm)}"
	    else
            lsty = σM ≥ ΩR ? :dash : :solid  
        	lines!(ax, tvals, d, linewidth=3.0, linestyle=lsty, label=L"%$(round(σM/ΩR, digits=3))")
	    	if showstd
        		band!(ax, tvals, d .- std_d, d .+ std_d, transparency=true)
	    	end
	    end
	    #xlims!(ax, x_min, x_max)
	    #ylims!(ax, y_min, y_max)
    end
    fig[2,1] = Legend(fig, ax, L"\sigma_\mathrm{M} / \Omega_\mathrm{R}", orientation=:horizontal)
    fig
end

function plot_fit(σx=120)

    Rvals = [0.05, 0.1, 0.2, 0.3]
    σMvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    axs = [Axis(fig[i,j]) for i = 1:2, j = 1:2]

    for (ax, ΩR) in zip(axs, Rvals)
        ax.xticklabelsvisible = false
        ax.title = L"\Omega_R = %$(ΩR)\;\mathrm{eV}"
        xlims!(ax, 0, 300)
        ylims!(ax, -5, 200)
    end

    for (ax, ΩR) in zip(axs, Rvals)

	    Rstr = replace(string(ΩR), "." => "p") 
        for σM in σMvals

	        sm = replace(string(σM), "." => "p") 
        
            rg = 1:61

            d = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "$(Int(σx))_avg_d")[rg] .- σx/10
	        tvals = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "d_tvals")[rg] .* 1000

            df = DataFrame(X = tvals, Y = d) 
            ols = lm(@formula(Y ~ X), df)

            fit = predict(ols)

            # y = mx + b
            b = ols.model.pp.beta0[1]
            m = ols.model.pp.beta0[2]

            scatter!(ax, tvals[rg], d, linewidth=3.0, label = L"%$(σM)")
            lines!(ax, tvals[rg], fit, linewidth=3.0, label = L"%$(σM)")
        end
    end

    axs[2].xticklabelsvisible = true
    axs[4].xticklabelsvisible = true
    Label(fig[1:end, 0], L"d", rotation = pi/2, justification=:right, textsize = 25)
    Label(fig[3, 1:2], "Time (fs)", justification=:right, textsize = 20)
    #fig[4,1:2] = Legend(fig, axs[1], L"\delta \;(\mathrm{eV})", orientation=:horizontal)

    fig
end

function plot_slopes(σx=120)

    Rvals = [0.05, 0.1, 0.2, 0.3]
    σMvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]

    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel=L"\sigma_\mathrm{M}\;(eV)")
    ax2 = Axis(fig[2,1], xlabel=L"{\sigma_\mathrm{M}}/{\Omega_\mathrm{R}}", xticks=[0, 1, 3, 5, 7, 9])
    xlims!(ax1, -0.01, 0.5)
    xlims!(ax2, -0.1, 10)

    slopes = zeros(4,10)

    for (j,ΩR) in enumerate(Rvals)
	    Rstr = replace(string(ΩR), "." => "p") 
        for (k,σM) in enumerate(σMvals)

	        sm = replace(string(σM), "." => "p") 
            rg = 1:61

            d = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "$(Int(σx))_avg_d")[rg] * 10
	        tvals = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "d_tvals")[rg] .* 1000

            df = DataFrame(X = tvals, Y = d) 
            ols = lm(@formula(Y ~ X), df)
        
            # y = mx + b
            b = ols.model.pp.beta0[1]
            m = ols.model.pp.beta0[2]
            slopes[j,k] = m
        end
    end


    for i = axes(slopes,1)
        m1 = σMvals .< Rvals[i]
        m2 = m1 .== false
        scatter!(ax1, σMvals, slopes[i,:], label = L"%$(Rvals[i])")
        lines!(ax1, σMvals[m1], slopes[i,m1], color=Cycled(i))
        lines!(ax1, σMvals, slopes[i,:], linestyle=:dash, color=Cycled(i))

        scatter!(ax2, σMvals ./Rvals[i], slopes[i,:])
        lines!(ax2, σMvals[m1] ./Rvals[i], slopes[i,m1], color=Cycled(i))
        lines!(ax2, σMvals ./Rvals[i], slopes[i,:], linestyle=:dash, color=Cycled(i))
    end

    fig[1:2,2] = Legend(fig, ax1, L"\Omega_\mathrm{R}\;(\mathrm{eV})")#, orientation=:horizontal)
    Label(fig[1:end, 0], L"v_0\;(\mathrm{nm}\;\mathrm{fs}^{-1})", rotation = pi/2, justification=:right, textsize = 25)

    fig
end

function plot_slopes_vs_σx(ΩR=0.1)

	Rstr = replace(string(ΩR), "." => "p") 

    σxvals = [60, 120, 240, 360, 480]
    σMvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]

    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax = Axis(fig[1,1], ylabel=L"v_0\;(\mathrm{nm}\;\mathrm{fs}^{-1})", xlabel=L"{\sigma_\mathrm{M}}/{\Omega_\mathrm{R}}", 
    xticks=0:0.2:1.2, ylabelsize=30, title=L"\Omega_R = %$(ΩR)\; \mathrm{eV}")

    xlims!(ax, 0,1.2)

    slopes = zeros(5,10)

    for (j,σx) in enumerate(σxvals)
        for (k,σM) in enumerate(σMvals)

	        sm = replace(string(σM), "." => "p") 
            rg = 1:61

            d = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "$(Int(σx))_avg_d")[rg] * 10
	        tvals = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "d_tvals")[rg] .* 1000

            df = DataFrame(X = tvals, Y = d) 
            ols = lm(@formula(Y ~ X), df)
        
            # y = mx + b
            b = ols.model.pp.beta0[1]
            m = ols.model.pp.beta0[2]
            slopes[j,k] = m
        end
    end

    for i = axes(slopes,1)
        scatter!(ax, σMvals ./ ΩR, slopes[i,:], label = L"%$(σxvals[i])")
        lines!(ax, σMvals ./ ΩR , slopes[i,:], )
    end

    fig[2,1] = Legend(fig, ax, L"\sigma_x\;(\mathrm{nm})", orientation=:horizontal)

    fig
end

function animate(ΩR, σMvals, σx)

	Rstr = replace(string(ΩR), "." => "p") 
    fig = Figure()
    positions = [i*0.5 - 25  for i = 1:100]
    time = Observable(1)

    #### FIRST
    axQ1 = Axis(fig[1,1], xticks=-10:2:10, xlabelsize=25,
    title= L"\sigma_\mathrm{M} / \Omega_\mathrm{R}\; = \; %$(round(σMvals[1]/ΩR, digits=2))", titlesize=20)
    xlims!(axQ1, -8, 8)
    axQ1.xticklabelsvisible = false


	smQ1 = replace(string(σMvals[1]), "." => "p") 
	tvalsQ1 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ1/out.h5"), "wvp_tvals")
    wvpQ1 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ1/out.h5"), "$(Int(σx))_avg_wvp")
    wvp_stdQ1 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ1/out.h5"), "$(Int(σx))_std_wvp")

    wvQ1 = @lift begin
        wvpQ1[:, $time]
    end

    σwvQ1 = @lift begin
        wvp_stdQ1[:,$time]
    end

	barplot!(axQ1, positions, wvQ1)
	errorbars!(axQ1, positions, wvQ1, σwvQ1, color = :red4, whiskerwidth = 10)

    ### Second
    axQ2 = Axis(fig[2,1], xticks=-10:2:10, xlabelsize=25, 
    title= L"\sigma_\mathrm{M} / \Omega_\mathrm{R}\; = \; %$(round(σMvals[2]/ΩR, digits=2))", titlesize=20)
    xlims!(axQ2, -8, 8)
    axQ2.xticklabelsvisible = false


	smQ2 = replace(string(σMvals[2]), "." => "p") 
	tvalsQ2 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ2/out.h5"), "wvp_tvals")
    wvpQ2 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ2/out.h5"), "$(Int(σx))_avg_wvp")
    wvp_stdQ2 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ2/out.h5"), "$(Int(σx))_std_wvp")

    wvQ2 = @lift begin
        wvpQ2[:, $time]
    end

    σwvQ2 = @lift begin
        wvp_stdQ2[:,$time]
    end

	barplot!(axQ2, positions, wvQ2)
	errorbars!(axQ2, positions, wvQ2, σwvQ2, color = :red4, whiskerwidth = 10)

    #### THIRD
    axQ3 = Axis(fig[3,1], xticks=-10:2:10, xlabelsize=25, 
    title= L"\sigma_\mathrm{M} / \Omega_\mathrm{R}\; = \; %$(round(σMvals[3]/ΩR, digits=2))", titlesize=20)
    xlims!(axQ3, -8, 8)
    axQ3.xticklabelsvisible = false


	smQ3 = replace(string(σMvals[3]), "." => "p") 
	tvalsQ3 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ3/out.h5"), "wvp_tvals")
    wvpQ3 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ3/out.h5"), "$(Int(σx))_avg_wvp")
    wvp_stdQ3 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ3/out.h5"), "$(Int(σx))_std_wvp")

    wvQ3 = @lift begin
        wvpQ3[:, $time]
    end

    σwvQ3 = @lift begin
        wvp_stdQ3[:,$time]
    end

	barplot!(axQ3, positions, wvQ3)
	errorbars!(axQ3, positions, wvQ3, σwvQ3, color = :red4, whiskerwidth = 10)

    #### FOURTH
    axQ4 = Axis(fig[4,1], xlabel="Distance (μm)", xticks=-10:2:10, xlabelsize=25, 
    title= L"\sigma_\mathrm{M} / \Omega_\mathrm{R}\; = \; %$(round(σMvals[4]/ΩR, digits=2))", titlesize=20)
    xlims!(axQ4, -8, 8)


	smQ4 = replace(string(σMvals[4]), "." => "p") 
	tvalsQ4 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ4/out.h5"), "wvp_tvals")
    wvpQ4 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ4/out.h5"), "$(Int(σx))_avg_wvp")
    wvp_stdQ4 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ4/out.h5"), "$(Int(σx))_std_wvp")

    wvQ4 = @lift begin
        wvpQ4[:, $time]
    end

    σwvQ4 = @lift begin
        wvp_stdQ4[:,$time]
    end

	barplot!(axQ4, positions, wvQ4)
	errorbars!(axQ4, positions, wvQ4, σwvQ4, color = :red4, whiskerwidth = 10)

    Label(fig[1:end, 0], L"|\Psi|^2", textsize = 25, justification=:right, rotation = pi/2)
    Label(fig[1:end, 2], @lift("t = $(5*($time - 1)) fs"), textsize = 25, justification=:right, rotation = -pi/2)

    framerate = 20
    timestamps = 1:101

    record(fig, "wvp.gif", timestamps;
            framerate = framerate) do t
        time[] = t
    end
end

function animate_dif(ΩR, σMvals, σx)

	Rstr = replace(string(ΩR), "." => "p") 
    fig = Figure()
    positions = [i*0.5 - 25  for i = 1:100]
    time = Observable(1)

    #### FIRST
    axQ1 = Axis(fig[1,1], xticks=-10:2:10, xlabelsize=25, ylabel=L"|\Psi|^2", ylabelsize=25,
    title= L"\sigma_\mathrm{M} / \Omega_\mathrm{R}\; = \; %$(round(σMvals[1]/ΩR, digits=2))", titlesize=20)
    xlims!(axQ1, -8, 8)

	smQ1 = replace(string(σMvals[1]), "." => "p") 
	tvalsQ1 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ1/out.h5"), "wvp_tvals")
    wvpQ1 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ1/out.h5"), "$(Int(σx))_avg_wvp")
    wvp_stdQ1 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ1/out.h5"), "$(Int(σx))_std_wvp")

	smQ2 = replace(string(σMvals[2]), "." => "p") 
	tvalsQ2 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ2/out.h5"), "wvp_tvals")
    wvpQ2 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ2/out.h5"), "$(Int(σx))_avg_wvp")
    wvp_stdQ2 = h5read(joinpath(@__DIR__,"Nm500_Nc500_a10_1/R$Rstr/sm$smQ2/out.h5"), "$(Int(σx))_std_wvp")

    wvQ1 = @lift begin
        wvpQ1[:, $time] .- wvpQ2[:, $time]
    end

	barplot!(axQ1, positions, wvQ1)

    #Label(fig[2, 1], @lift("t = $(5*($time - 1)) fs"), textsize = 25, justification=:right)

    framerate = 20
    timestamps = 1:101

    record(fig, "wvp_dif.gif", timestamps;
            framerate = framerate) do t
        time[] = t
    end
end

