#=====================================================================================================#
# This is Bayesian inference by MultiNest based on Magglow 
# MultiNest: https://github.com/JohannesBuchner/MultiNest
# This part is edited by Kaori Obayashi
# Detailed information will be written in Kusafuka, Obayashi, and Asano (2025)
#=====================================================================================================#
include("../src/Magglow.jl")
using .Magglow
using PyCall
using PyPlot
pygui(false) 
using DelimitedFiles
using Printf
using Statistics
using JSON
pymultinest = pyimport("pymultinest.solve")

# Getting arguments
P_line = 20 # Define from number of P line in MagglowNest_sample.sh
P_column = 5 # Define from number of P colum in MagglowNest_sample.sh
version = ARGS[1]
comment = ARGS[2]
args_after5 = length(ARGS[5:end])
expept_P = args_after5 - P_line*P_column
num_element = expept_P/3
freq_str = ARGS[5:Int(5 + num_element - 1)]
data_file = ARGS[Int(5 + num_element):Int(5 + 2 * num_element - 1)]
color = ARGS[Int(5 + 2 * num_element):Int(5 + 3 * num_element- 1)]
P_str = ARGS[Int(5 + 3 * num_element):end]

freq = [parse(Float64, f) for f in freq_str]
P = zeros(Float64, P_line, P_column)
for (i, row) in enumerate(eachrow(P))
    row .= [parse(Float64, P_str[(i - 1) * 5 + j]) for j in 1:5]
end

# Getting observed data
function extract_data(filename::String)
    time_list = Float64[]
    flux_list = Float64[]
    time_lower_err_list = Float64[]
    time_upper_err_list = Float64[]
    flux_lower_err_list = Float64[]
    flux_upper_err_list = Float64[]

    open(filename, "r") do file
        for line in eachline(file)
            values = split(line)
            time = parse(Float64, values[1])
            flux = parse(Float64, values[2])
            time_lower_err = parse(Float64, values[3])
            time_upper_err = parse(Float64, values[4])
            flux_lower_err = parse(Float64, values[5])
            flux_upper_err = parse(Float64, values[6])

            push!(time_list, time)
            push!(flux_list, flux)
            push!(time_lower_err_list, time_lower_err)
            push!(time_upper_err_list, time_upper_err)
            push!(flux_lower_err_list, flux_lower_err)
            push!(flux_upper_err_list, flux_upper_err)
        end
    end

    return time_list, flux_list, time_lower_err_list, time_upper_err_list, flux_lower_err_list, flux_upper_err_list
end

time = []
flux = []
time_lower_err = []
time_upper_err = []
flux_lower_err = []
flux_upper_err = []

time_log = []
flux_log = []
time_lower_err_log = []
time_upper_err_log = []
flux_lower_err_log = []
flux_upper_err_log = []
flux_err_log = []

for (index, file) in enumerate(data_file)
    time_list, flux_list, time_lower_err_list, time_upper_err_list, flux_lower_err_list, flux_upper_err_list = extract_data(file)
    println(index, file)

    push!(time, time_list)
    push!(time_lower_err, time_lower_err_list)
    push!(time_upper_err, time_upper_err_list)
    push!(flux, flux_list)
    push!(flux_lower_err, flux_lower_err_list)
    push!(flux_upper_err, flux_upper_err_list)

    t_log = log10.(time_list)
    f_log = log10.(flux_list)
    t_pos_err_log = log10.(time_list .+ time_upper_err_list) .- t_log
    f_pos_err_log = log10.(flux_list .+ flux_upper_err_list) .- f_log
    t_neg_err_log = t_log .- log10.(abs.(time_list .- time_lower_err_list))
    f_neg_err_log = f_log .- log10.(abs.(flux_list .- flux_lower_err_list))

    push!(time_log, t_log)
    push!(time_lower_err_log, t_neg_err_log)
    push!(time_upper_err_log, t_pos_err_log)
    push!(flux_log, f_log)
    push!(flux_lower_err_log, f_neg_err_log)
    push!(flux_upper_err_log, f_pos_err_log)
    push!(flux_err_log, (f_neg_err_log .+ f_pos_err_log) ./ 2)
    
end
flat_time_list = sort(unique(item for sublist in time for item in sublist))

# Organize all/search/minimu/maximum parameters
all_param_names = ["log(E0)", "log(G0)", "log(S0)", "log(D0)", "log(n0)", "k", "log(eps_e)", "log(eps_B)", "p", "log(f_e)", "theta_j", "beta", "log(eps_e(RS))", "log(eps_B(RS))", "p(RS)", "log(f_e(RS))", "eps_p", "log(f_p)", "eps_p(RS)", "log(f_p(RS))"]
parameters = []
bounds_lower = []
bounds_upper = []
for i in eachindex(all_param_names)
    if P[i, 1] == 0
        push!(parameters, all_param_names[i])
        push!(bounds_lower, P[i, 3])
        push!(bounds_upper, P[i, 4])
        println(all_param_names[i] * "=[" * string(P[i, 3]) * "," * string(P[i, 4]) * "]" )
    else
        println(all_param_names[i] * " = " * string(P[i, 2]))
    end
end

const z=parse(Float64,ARGS[3]) # redshift 
const DL=parse(Float64,ARGS[4]) # luminosity distance [cm]
const t_in = flat_time_list
const nu_in=freq

function Lightcurve!(Syn,Synp,SSC,SynRS,SynpRS,SSCRS, i)
    rcParams = PyDict(matplotlib["rcParams"])
    rcParams["font.family"] ="Times New Roman"
    rcParams["mathtext.fontset"] = "stix" 
    rcParams["xtick.direction"] = "in"     
    rcParams["ytick.direction"] = "in"     
    rcParams["xtick.minor.visible"] = true 
    rcParams["ytick.minor.visible"] = true 
    rcParams["xtick.major.width"] = 1.5   
    rcParams["ytick.major.width"] = 1.5   
    rcParams["xtick.minor.width"] = 1.0   
    rcParams["ytick.minor.width"] = 1.0    
    rcParams["xtick.major.size"] = 6       
    rcParams["ytick.major.size"] = 6      
    rcParams["xtick.minor.size"] = 3      
    rcParams["font.size"] = 15            
    rcParams["axes.linewidth"] = 1.0     

    fig = plt.figure(figsize=(12,9))
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    ax = fig.add_subplot(111)

    ax.set_xscale("log")
    ax.set_yscale("log")

    coef = [1e0, 1e0, 1e0]
    colo = ["red", "blue", "green"]
    h1 = []
    for m in eachindex(nu_in)
        model = (Syn[:,m].+SynRS[:,m].+SSC[:,m].+SSCRS[:,m]).*(10^10)
        model = model .*10^19
        ax.errorbar(time[m], coef[m].*flux[m], 
            xerr=(time_lower_err[m], time_upper_err[m]), 
            yerr=(flux_lower_err[m], flux_upper_err[m]), 
            fmt=".", color=colo[m])
        a,=ax.plot(t_in, coef[m].*model,
            label=@sprintf("%g Hz (× %g)", nu_in[m], coef[m]),
            color=colo[m])
        push!(h1, a) 
    end

    fig.text(0.04, 0.5, raw"$F_ν$ [erg $cm^{-2}$ $s^{-1}$ $Hz^{-1}$]", va="center", rotation="vertical", fontsize=20) 
    fig.text(0.5, 0.04, raw"Observed Time $t_{\rm obs}$ [s]", ha="center", fontsize=20) 
        
    lab1=[q.get_label() for q in h1]
    fig.legend(h1,lab1,bbox_to_anchor=(0.5, 0.98), loc="upper center", ncol=4, fontsize=13)
    fig.savefig(version*"/lightcurve.pdf",dpi=300)
end

function Fluxvalue!(m, target_time, Syn,Synp,SSC,SynRS,SynpRS,SSCRS)
    Fluxlist = Syn[:,m].+SynRS[:,m].+SSC[:,m].+SSCRS[:,m]
    match_index = [argmin(abs.(t_in .- t)) for t in target_time]
    return Fluxlist[match_index]
end

function myprior!(cube)
    for n in eachindex(parameters)
        cube[n] = cube[n]*(bounds_upper[n] - bounds_lower[n]) + bounds_lower[n]
    end
    return cube
end

function myloglike!(cube)
    err=Float16[]
    n = 1
    for i in eachindex(all_param_names)
        if P[i, 1] == 0
            P[i, 5] = cube[n]
            n += 1
        else
            P[i, 5] = P[i, 2]
        end
    end
    
    # Input parameters
    Input=Float64[
        10^P[1,5]       ,    # 1 initial energy log(E0 [erg])
        10^P[2,5]       ,    # 2 initial Lorentz factor log(G0 at R0)
        10^P[3,5]       ,    # 3 initial magnetization log(S0)
        3e10*10^P[4,5]  ,    # 4 initial radial width log(D0 [cm])
        10^P[5,5]       ,    # 5 CSM number density log(n0 [cm^(k-3)])
        P[6,5]          ,    # 6 density slope k (ISM:0 --- 2:wind) 
        10^P[7,5]       ,    # 7 energy fraction of accelerated electrons log(epsilon_e)
        10^P[8,5]       ,    # 8 energy fraction of Weibel induced magnetic field log(epsilon_B)
        P[9,5]          ,    # 9 PD spectral index 2 < p_FS < 3, p = 2.2
        10^P[10,5]      ,    # 10 number fraction of accelerated electrons log(f_e)
        P[11,5]         ,    # 11 opening angle theta_j = 1.0
        P[12,5]         ,    # 12 viewing angle theta_o = 0.0
        10^P[13, 5]     ,    # 13 energy fraction of accelerated electrons log(epsiron_e,RS)
        10^P[14, 5]     ,    # 14 energy fraction of Weibel induced magnetic field log(epsiron_B,RS)
        P[15, 5]        ,    # 15 PD spectral index 2 < p_RS < 3
        10^P[16, 5]     ,    # 16 number fraction of accelerated electrons log(f_e,RS)
        10^P[17, 5]     ,    # 17 energy fraction of accelerated protons log(epsiron_p)
        10^P[18, 5]     ,    # 18 number fraction of accelerated protons log(f_p)
        10^P[19, 5]     ,    # 19 energy fraction of accelerated protons log(epsiron_p,RS)
        10^P[20, 5]          # 20 number fraction of accelerated protons log(f_p,RS)
    ]

    # Output e/p synchrotron & SSC flux 
    Output=[
        zeros(Float64,length(t_in),length(nu_in)),    # 1 e-synchrotron    FS
        zeros(Float64,length(t_in),length(nu_in)),    # 2 p-synchrotron    FS
        zeros(Float64,length(t_in),length(nu_in)),    # 3 e-SSC            FS
        zeros(Float64,length(t_in),length(nu_in)),    # 4 e-synchrotron    RS
        zeros(Float64,length(t_in),length(nu_in)),    # 5 p-synchrotron    RS
        zeros(Float64,length(t_in),length(nu_in)),    # 6 e-SSC            RS
        zeros(Float64,length(t_in),length(nu_in)),    # 7 pp e-neutrino    FS
        zeros(Float64,length(t_in),length(nu_in)),    # 8 pp mu-neutrino   FS
        zeros(Float64,length(t_in),length(nu_in)),    # 9 pp pi0 gamma     FS
        zeros(Float64,length(t_in),length(nu_in)),    # 10 pg e-neutrino   FS
        zeros(Float64,length(t_in),length(nu_in)),    # 11 pg mu-neutrino  FS
        zeros(Float64,length(t_in),length(nu_in)),    # 12 pg pi0 gamma    FS
        zeros(Float64,length(t_in),length(nu_in)),    # 13 pp e-neutrino   RS
        zeros(Float64,length(t_in),length(nu_in)),    # 14 pp mu-neutrino  RS
        zeros(Float64,length(t_in),length(nu_in)),    # 15 pp pi0 gamma    RS
        zeros(Float64,length(t_in),length(nu_in)),    # 16 pg e-neutrino   RS
        zeros(Float64,length(t_in),length(nu_in)),    # 17 pg mu-neutrino  RS
        zeros(Float64,length(t_in),length(nu_in))     # 18 pg pi0 gamma    RS
    ]

    # Which is calculated?
    is_calc=[
        true  , # 1 e-syn
        false , # 2 p-syn
        false , # 3 e-SSC
        false , # 4 pp
        false   # 5 pg
    ]

    # Model calculation
    @time Magnetic_Bullet_Afterglow!(z,DL,t_in,nu_in,Input,Output,is_calc)

    # Likelihood
    for m in eachindex(nu_in)
        model = Fluxvalue!(m, time[m], Output[1],Output[2],Output[3],Output[4],Output[5],Output[6])
        model = model .*1e26 #mJy=[10^(-26) erg cm^-2 s^-1 Hz^-1]
        log_model = log10.(model)
        sigma = flux_err_log[m].^2
        xerr = -0.5*sum((flux_log[m] .- log_model).^2 ./ sigma .+ log.(sigma))
        push!(err,xerr)
    end
    return sum(err)
end 

function Record_txt(T, result, InputP, filename)
    open(filename, "w") do file
        # comment
        println(file, "---  Comment  ---")
        println(file, comment)
        println(file)

        # Time of estimation
        println(file, "---Elapsed time---")
        println(file, "$T seconds")
        println(file)

        # data file
        println(file, "--- Data file ---")
        for m in eachindex(data_file)
            println(file, data_file[m])
        end
        println(file)

        # setting
        println(file, "---  Setting  ---")
        println(file, "redshift         :$z")
        println(file, "luminosity[cm]   :$DL")
        for i in eachindex(all_param_names)
            if P[i, 1] == 0
                println(file, all_param_names[i], ": bound = (", @sprintf("%g", P[i, 3]), ", ", @sprintf("%g", P[i, 4]), ")")
            else
                println(file, all_param_names[i], ": fixed = " * string(InputP[i, 2]))
            end
        end
        println(file)

        # evidence 
        println(file, "---  Results  ---")
        println(file, "evidence: $(round(result["logZ"], digits=1)) +- $(round(result["logZerr"], digits=1))")
        println(file)

        # params
        for (name, col) in zip(parameters, eachcol(result["samples"]))
            println(file, @sprintf("%15s : %.3f +- %.3f", name, mean(col), std(col)))
        end
    end
end

function Results!(T, result)
    n = 1
    for i in eachindex(all_param_names)
        if P[i, 1] == 0
            P[i, 5] = mean(result["samples"][:, n])
            n += 1
        else
            P[i, 5] = P[i, 2]
        end
    end

    # Input parameters
    Input=Float64[
        10^P[1,5]       ,    # 1 initial energy log(E0 [erg])
        10^P[2,5]       ,    # 2 initial Lorentz factor log(G0 at R0)
        10^P[3,5]       ,    # 3 initial magnetization log(S0)
        3e10*10^P[4,5]  ,    # 4 initial radial width log(D0 [cm])
        10^P[5,5]       ,    # 5 CSM number density log(n0 [cm^(k-3)])
        P[6,5]          ,    # 6 density slope k (ISM:0 --- 2:wind) 
        10^P[7,5]       ,    # 7 energy fraction of accelerated electrons log(epsilon_e)
        10^P[8,5]       ,    # 8 energy fraction of Weibel induced magnetic field log(epsilon_B)
        P[9,5]          ,    # 9 PD spectral index 2 < p_FS < 3, p = 2.2
        10^P[10,5]      ,    # 10 number fraction of accelerated electrons log(f_e)
        P[11,5]         ,    # 11 opening angle theta_j = 1.0
        P[12,5]         ,    # 12 viewing angle theta_o = 0.0
        10^P[13, 5]     ,    # 13 energy fraction of accelerated electrons log(epsiron_e,RS)
        10^P[14, 5]     ,    # 14 energy fraction of Weibel induced magnetic field log(epsiron_B,RS)
        P[15, 5]        ,    # 15 PD spectral index 2 < p_RS < 3
        10^P[16, 5]     ,    # 16 number fraction of accelerated electrons log(f_e,RS)
        10^P[17, 5]     ,    # 17 energy fraction of accelerated protons log(epsiron_p)
        10^P[18, 5]     ,    # 18 number fraction of accelerated protons log(f_p)
        10^P[19, 5]     ,    # 19 energy fraction of accelerated protons log(epsiron_p,RS)
        10^P[20, 5]          # 20 number fraction of accelerated protons log(f_p,RS)
    ]

    # Output e/p synchrotron & SSC flux 
    Output=[zeros(Float64,length(t_in),length(nu_in)) for i=1:18]

    # Which is calculated?
    is_calc=[
        true  , # 1 e-syn
        false , # 2 p-syn
        false , # 3 e-SSC
        false , # 4 pp
        false   # 5 pg
    ]

    # Record set parameters and results
    println("recording to txt")
    Record_txt(T, result, P, version*"/record.txt")

    # Model calculation
    println("calculation")
    @time Magnetic_Bullet_Afterglow!(z,DL,t_in,nu_in,Input,Output,is_calc)

    # Draw light curve
    println("drawing LC")
    Lightcurve!(Output[1],Output[2],Output[3],Output[4],Output[5],Output[6], 0)
end

function PlotPosterior!(num_posterior_samples, colo)
    # setting
    fig = plt.figure(figsize=(12,9))
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.set_yscale("log")

    # data
    for m in eachindex(nu_in)
        ax.errorbar(time[m], flux[m], 
            xerr=(time_lower_err[m], time_upper_err[m]), 
            yerr=(flux_lower_err[m], flux_upper_err[m]), 
            fmt=".", color=colo[m])
    end

    # posterior
    posterior_samples = readdlm(version*"/3-post_equal_weights.dat")
    println(size(posterior_samples))  # (# of samples, # of params)
    selected_samples = posterior_samples[rand(1:size(posterior_samples, 1), num_posterior_samples), :]
    for sample in eachrow(selected_samples)
        n = 1
        for i in eachindex(all_param_names)
            if P[i, 1] == 0
                P[i, 5] = sample[n]
                n += 1
            else
                P[i, 5] = P[i, 2]
            end
        end
        Input=Float64[
            10^P[1,5]       ,    # 1 initial energy log(E0 [erg])
            10^P[2,5]       ,    # 2 initial Lorentz factor log(G0 at R0)
            10^P[3,5]       ,    # 3 initial magnetization log(S0)
            3e10*10^P[4,5]  ,    # 4 initial radial width log(D0 [cm])
            10^P[5,5]       ,    # 5 CSM number density log(n0 [cm^(k-3)])
            P[6,5]          ,    # 6 density slope k (ISM:0 --- 2:wind) 
            10^P[7,5]       ,    # 7 energy fraction of accelerated electrons log(epsilon_e)
            10^P[8,5]       ,    # 8 energy fraction of Weibel induced magnetic field log(epsilon_B)
            P[9,5]          ,    # 9 PD spectral index 2 < p_FS < 3, p = 2.2
            10^P[10,5]      ,    # 10 number fraction of accelerated electrons log(f_e)
            P[11,5]         ,    # 11 opening angle theta_j = 1.0
            P[12,5]         ,    # 12 viewing angle theta_o = 0.0
            10^P[13, 5]     ,    # 13 energy fraction of accelerated electrons log(epsiron_e,RS)
            10^P[14, 5]     ,    # 14 energy fraction of Weibel induced magnetic field log(epsiron_B,RS)
            P[15, 5]        ,    # 15 PD spectral index 2 < p_RS < 3
            10^P[16, 5]     ,    # 16 number fraction of accelerated electrons log(f_e,RS)
            10^P[17, 5]     ,    # 17 energy fraction of accelerated protons log(epsiron_p)
            10^P[18, 5]     ,    # 18 number fraction of accelerated protons log(f_p)
            10^P[19, 5]     ,    # 19 energy fraction of accelerated protons log(epsiron_p,RS)
            10^P[20, 5]          # 20 number fraction of accelerated protons log(f_p,RS)
        ]
        
        Output=[zeros(Float64,length(t_in),length(nu_in)) for i=1:18]

        is_calc=[
            true  , # 1 e-syn
            false , # 2 p-syn
            false , # 3 e-SSC
            false , # 4 pp
            false   # 5 pg
        ]

        Magnetic_Bullet_Afterglow!(z, DL, t_in, nu_in, Input, Output, is_calc)
        Syn = Output[1]
        Synp = Output[2]
        SSC = Output[3]
        SynRS = Output[4]
        SynpRS = Output[5]
        SSCRS = Output[6]
        for m in eachindex(nu_in)
            model = (Syn[:,m].+SynRS[:,m].+SSC[:,m].+SSCRS[:,m])
            model = model .*1e26 #mJy=[10^(-26) erg cm^-2 s^-1 Hz^-1]
            ax.plot(t_in, model, color=colo[m], alpha=0.05)  # e-synchrotron FS
        end
    end

    # setting for label
    fig.text(0.04, 0.5, raw"$F_ν$ [erg $cm^{-2}$ $s^{-1}$ $Hz^{-1}$]", va="center", rotation="vertical", fontsize=20) 
    fig.text(0.5, 0.04, raw"Observed Time $t_{\rm obs}$ [s]", ha="center", fontsize=20)
    fig.savefig(version*"/toyevent_sample.pdf",dpi=300)
end

# #=====================================================================================================#
elapsed_time = @elapsed begin
    # Run MultiNest !!
    # ignore signal handling
    pyimport("signal").signal(pyimport("signal").SIGINT, pyimport("signal").SIG_IGN)
    @time result_run = pymultinest.solve(
        LogLikelihood=myloglike!, 
        Prior=myprior!,
        n_dims=length(parameters), 
        outputfiles_basename=version*"/3-", 
        verbose=true,
        n_live_points=450, # default=500
        evidence_tolerance=0.6, # default=0.5
        sampling_efficiency = 0.8 # defalt: none
    )
    pyimport("signal").signal(pyimport("signal").SIGINT, pyimport("signal").SIG_DFL)
    println(" == End the MultiNest ==")

    prefix = version*"/3-"
    filename = "$(prefix)params.json"
    open(filename, "w") do file
        JSON.print(file, parameters, 2)
    end
end
# Output results and Plot lightcurve
Results!(elapsed_time, result_run)
println(" == Finnish MultiNest ==")
PlotPosterior!(100, color)
println(" == CLOSE ==")
exit()