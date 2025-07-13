include("../src/Magglow.jl")
using .Magglow
using PyCall
using PyPlot
using Printf
@pyimport matplotlib.animation as animation
using Revise

# Source Distance
z=1.4607 # redshift
DL=3.3020e28 # luminosity distance [cm]

# Input Observed Time
t_in=10 .^ range(log10(1e1), log10(1e6), length=100)

# Input Observed Frequency
nu_in=10 .^ range(log10(1e9), log10(1e19), length=100)

# Input parameters
Input=Float64[
    6e54   ,    # 1 initial energy E0 [erg]
    12    ,    # 2 initial Lorentz factor G0 at R0
    5e0   ,    # 3 initial magnetization S0
    1400*3e10  ,    # 4 initial radial width D0 [cm], which should be smaller than the initial radius R0
    2e-1    ,    # 5 CSM number density n0 at R0 [/cc]
    2      ,    # 6 density slope k (ISM:0 --- 2:wind)
    0.1    ,    # 7 energy fraction of accelerated electrons epsiron_e
    0.001   ,    # 8 energy fraction of Weibel induced magnetic field epsiron_B
    2.5    ,    # 9 PD spectral index 2 < p < 3
    0.4    ,    # 10 number fraction of accelerated electrons f_e
    0.025    ,    # 11 opening angle theta_j
    0.0    ,    # 12 viewing angle theta_o
    0.1    ,    # 13 energy fraction of accelerated electrons epsiron_e,RS
    0.001   ,    # 14 energy fraction of Weibel induced magnetic field epsiron_B,RS
    2.01    ,    # 15 PD spectral index 2 < p_RS < 3
    1.0    ,    # 16 number fraction of accelerated electrons f_e,RS
    0.1    ,    # 17 energy fraction of accelerated protons epsiron_p,FS
    0.01        # 18 number fraction of accelerated protons f_p,FS
    0.1    ,    # 19 energy fraction of accelerated protons epsiron_p,RS
    0.01        # 20 number fraction of accelerated protons f_p,RS
]

# Output leptonic/hadronic emission
Output=[zeros(Float64,length(t_in),length(nu_in)) for i=1:18]

# Which is calculated?
is_calc=[
    true  , # 1 e-syn
    false , # 2 p-syn
    true  , # 3 e-SSC
    true , # 4 pp
    false   # 5 pg
]

Magnetic_Bullet_Afterglow!(z,DL,t_in,nu_in,Input,Output,is_calc)

function Spectrum!(Syn,Synp,SSC,SynRS,SynpRS,SSCRS)
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

    fig = plt.figure(figsize=(9, 6))
    fig.tight_layout(rect=[0,0,1,0.96])
    plt.subplots_adjust(wspace=0.4, hspace=0.3)
    ax = fig.add_subplot(111)

    function Mov(i)
        ax.cla()

        ax.set_xscale("log")
        ax.set_yscale("log")
        
        a1,=ax.plot(nu_in,1e26.*Syn[i,:],label="FS",color="red",linestyle="dashed")
        # a2,=ax.plot(nu_in,nu_in.*SSC[i,:],label="e-SSC ",color="green",linestyle="dasehd")
        # a3,=ax.plot(nu_in,nu_in.*Synp[i,:],label="p-synchrotron",color="blue",linestyle="dashed")
        # a4,=ax.plot(nu_in,nu_in.*(Syn[i,:].+Synp[i,:].+SSC[i,:]),label="total",color="blue",linestyle="dashed")

        a2,=ax.plot(nu_in,1e26.*SynRS[i,:],label="RS",color="blue",linestyle="dotted")
        # ax.plot(nu_in,nu_in.*SSCRS[i,:],label="e-SSC ",color="green",linestyle="dotted")
        # ax.plot(nu_in,nu_in.*SynpRS[i,:],label="p-synchrotron",color="blue",linestyle="dotted")
        # ax.plot(nu_in,nu_in.*(SynRS[i,:].+SynpRS[i,:].+SSCRS[i,:]),label="total",color="blue",linestyle="dotted")
        
        ax.set_ylim(1e-5,1e1)

        ax.vlines(1e14,1e-24,1e-10)

        # ax.set_ylabel(raw"$νF_ν$ [erg $cm^{-2}$ $s^{-1}$]")
        ax.set_ylabel(raw"$F_ν$ [mJy]")
        ax.set_xlabel(raw"frequency [Hz]")
        
        function eV(x)
            return x * 6.626e-27 / 1.6e-12
        end

        function Hz(x)
            return x / 6.626e-27 * 1.6e-12
        end

        aax = ax.secondary_xaxis("top", functions=(eV, Hz))
        aax.set_xlabel("energy [eV]")

        h1=[a1,a2]
        lab1=[q.get_label() for q in h1]
        fig.legend(h1,lab1,bbox_to_anchor=(0.25, 0.86), loc="upper center", ncol=1, fontsize=13)
        
        time=@sprintf("%e",t_in[i])

        fig.suptitle("Time = "*time*" s")
    end

    Spec = animation.FuncAnimation(fig, Mov, frames = 1 : length(t_in), interval = 100)
    Spec.save("Spec_sample.mp4", writer=animation.FFMpegWriter(fps=10))
end

Spectrum!(Output[1],Output[2],Output[3],Output[4],Output[5],Output[6])