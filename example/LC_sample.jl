include("../src/Magglow.jl")
using .Magglow
using PyCall
using PyPlot

# Source Distance
z=1.0 # redshift
DL=2.1e28 # luminosity distance [cm]

# Input Observed Time
t_in=10 .^ range(log10(1e0), log10(1e7), length=100)

# Input Observed Frequency
nu_in=Float64[
    1e9  , # Radio
    1e14 , # Optical
    1e17   # X-ray
]

# Input parameters
Input=Float64[
    1e54   ,    # 1 initial energy E0 [erg]
    100    ,    # 2 initial Lorentz factor G0 at R0
    5e0   ,    # 3 initial magnetization S0
    300*3e10  ,    # 4 initial radial width D0 [cm], which should be smaller than the initial radius R0
    1e0    ,    # 5 CSM number density n0 at R0 [/cc]
    0      ,    # 6 density slope k (ISM:0 --- 2:wind)
    0.1    ,    # 7 energy fraction of accelerated electrons epsiron_e
    0.001   ,    # 8 energy fraction of Weibel induced magnetic field epsiron_B
    2.2    ,    # 9 PD spectral index 2 < p < 3
    1.0    ,    # 10 number fraction of accelerated electrons f_e
    0.1    ,    # 11 opening angle theta_j
    0.0    ,    # 12 viewing angle theta_o
    0.1    ,    # 13 energy fraction of accelerated electrons epsiron_e,RS
    0.001   ,    # 14 energy fraction of Weibel induced magnetic field epsiron_B,RS
    2.2    ,    # 15 PD spectral index 2 < p_RS < 3
    1.0    ,    # 16 number fraction of accelerated electrons f_e,RS
    0.1    ,    # 17 energy fraction of accelerated protons epsiron_p
    0.01   ,    # 18 number fraction of accelerated protons f_p
    0.1    ,    # 19 energy fraction of accelerated protons epsiron_p,RS
    0.01        # 20 number fraction of accelerated protons f_p,RS
]

# Output leptonic/hadronic emission
Output=[zeros(Float64,length(t_in),length(nu_in)) for i=1:18]

# Which is calculated?
is_calc=[
    true  , # 1 e-syn
    false , # 2 p-syn
    false , # 3 e-SSC
    false , # 4 pp
    false   # 5 pg
]

Magnetic_Bullet_Afterglow!(z,DL,t_in,nu_in,Input,Output,is_calc)

function Lightcurve!(Syn,Synp,SSC,SynRS,SynpRS,SSCRS)
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

    a1,=ax.plot(t_in,nu_in[1].*(Syn[:,1].+Synp[:,1].+SSC[:,1]),label=raw"Radio $10^{9}$ Hz",color="red",linestyle="dashed")
    a2,=ax.plot(t_in,nu_in[2].*(Syn[:,2].+Synp[:,2].+SSC[:,2]),label=raw"Optical $10^{14}$ Hz",color="blue",linestyle="dashed")
    a3,=ax.plot(t_in,nu_in[3].*(Syn[:,3].+Synp[:,3].+SSC[:,3]),label=raw"X-ray $10^{18}$ Hz",color="green",linestyle="dashed")
    ax.plot(t_in,nu_in[1].*(SynRS[:,1].+SynpRS[:,1].+SSCRS[:,1]),label=raw"$3.3\times10^{14}$ Hz",color="red",linestyle="dotted")
    ax.plot(t_in,nu_in[2].*(SynRS[:,2].+SynpRS[:,2].+SSCRS[:,2]),label=raw"$4.8\times10^{14}$ Hz",color="blue",linestyle="dotted")
    ax.plot(t_in,nu_in[3].*(SynRS[:,3].+SynpRS[:,3].+SSCRS[:,3]),label=raw"$7.3\times10^{16}$ Hz",color="green",linestyle="dotted")
    ax.plot(t_in,nu_in[1].*(Syn[:,1].+SynRS[:,1].+SSC[:,1].+SSCRS[:,1]),label=raw"$3.3\times10^{14}$ Hz",color="red")
    ax.plot(t_in,nu_in[2].*(Syn[:,2].+SynRS[:,2].+SSC[:,2].+SSCRS[:,2]),label=raw"$4.8\times10^{14}$ Hz",color="blue")
    ax.plot(t_in,nu_in[3].*(Syn[:,3].+SynRS[:,3].+SSC[:,3].+SSCRS[:,3]),label=raw"$7.3\times10^{16}$ Hz",color="green")

    ax.set_xlim(1e0,1e6)
    ax.set_ylim(1e-27,1e-7)

    fig.supylabel(raw"$νF_ν$ [erg $\rm cm^{-2}$ $\rm s^{-1}$]",fontsize=20)
    fig.supxlabel(raw"Observed Time $t_{\rm obs}$ [s]",fontsize=20)
        
    h1=[a1,a2,a3]
    lab1=[q.get_label() for q in h1]
    fig.legend(h1,lab1,bbox_to_anchor=(0.5, 0.98), loc="upper center", ncol=4, fontsize=13)

    fig.savefig("LC_sample.pdf",dpi=300)
end

Lightcurve!(Output[1],Output[2],Output[3],Output[4],Output[5],Output[6])