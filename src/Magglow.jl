#=====================================================================================================#
# This is the Magnetic Bullet Afterglow (Magglow) code written in Julia
# Detailed information is written in Kusafuka, Obayashi, and Asano (2025)
#=====================================================================================================#
module Magglow
export Magnetic_Bullet_Afterglow!
#=====================================================================================================#
# Physical constants and calculation constants
const e=4.8e-10 # esu
const mpi=140e6 # eV
const mmu=106e6 # eV
const mp=938e6 # eV
const me=511e3 # eV
const c=3e10 # cm/s
const h=6.626e-27 # erg s
const sigma_T=0.665e-24 # cm^2
const k_B=1.38e-16 # erg/K
const alpha=asin(pi/4.0) # pitch angle (isotropic assumption)
const lambda=1.0-(mmu/mpi)^2 # 
const R0=1e13 # initial radius estimated from variability or just assumption (not so important)
const delta_t=0.1
const delta_nu=0.1
const n_theta=150
#=====================================================================================================#
function Magnetic_Bullet_Afterglow!(z,DL,time,freq,In,Out,is_calc) 

    # Observed time
    t=collect(LinRange(0,10,Int(10/delta_t)+1))
    t_obs=zeros(Float64,length(t))
    t_obs.=10 .^t

    # observed frequency
    nu=collect(LinRange(6,28,Int((28-6)/delta_nu)+1))
    nu_obs=zeros(Float64,length(nu))
    nu_obs.=10 .^nu

    # observed flux
    Fnu=[
        zeros(Float64,length(t_obs),length(nu_obs)),    # 1 e-synchrotron    FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 2 p-synchrotron    FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 3 e-SSC            FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 4 e-synchrotron    RS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 5 p-synchrotron    RS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 6 e-SSC            RS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 7 pp e-neutrino    FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 8 pp mu-neutrino   FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 9 pp pi0 gamma     FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 10 pg e-neutrino   FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 11 pg mu-neutrino  FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 12 pg pi0 gamma    FS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 13 pp e-neutrino   RS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 14 pp mu-neutrino  RS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 15 pp pi0 gamma    RS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 16 pg e-neutrino   RS
        zeros(Float64,length(t_obs),length(nu_obs)),    # 17 pg mu-neutrino  RS
        zeros(Float64,length(t_obs),length(nu_obs))     # 18 pg pi0 gamma    RS
    ]

    # MHD time step
    t=collect(LinRange(0,10,500+1))
    t_MHD=zeros(Float64,length(t))
    t_MHD.=10 .^t

    # Forward / Reverse shock properties
    timeFS=[2,1];timeRS=[2,1]
    GFS=zeros(Float64,length(t_MHD));GRS=zeros(length(t_MHD))
    RFS=zeros(Float64,length(t_MHD));RRS=zeros(length(t_MHD))
    NFS=zeros(Float64,length(t_MHD));NRS=zeros(length(t_MHD))
    EFS=zeros(Float64,length(t_MHD));ERS=zeros(length(t_MHD))
    BFS=zeros(Float64,length(t_MHD));BRS=zeros(length(t_MHD))
    DFS=zeros(Float64,length(t_MHD));DRS=zeros(length(t_MHD))

    # Forward & Reverse shock dynamics
    # Fire_Ball!(In[1],In[2],In[5],In[6],In[8],z,time,t_MHD,timeFS,GFS,RFS,NFS,EFS,BFS,DFS)
    Magnetic_Bullet!(In[1],In[2],In[3],In[4],In[5],In[6],In[8],In[14],z,time,t_MHD,timeFS,GFS,RFS,NFS,EFS,BFS,DFS,timeRS,GRS,RRS,NRS,ERS,BRS,DRS)

    # Forward shock radiation
    calculate_radiation!(z,DL,time,In[7],In[8],In[9],In[10],In[17],In[18],In[11],In[12],nu_obs,t_obs,timeFS,t_MHD,NFS,EFS,BFS,RFS,DFS,GFS,is_calc,Fnu[1],Fnu[2],Fnu[3],Fnu[7],Fnu[8],Fnu[9],Fnu[10],Fnu[11],Fnu[12])

    # Reverse shock radiation
    calculate_radiation!(z,DL,time,In[13],In[14],In[15],In[16],In[19],In[20],In[11],In[12],nu_obs,t_obs,timeRS,t_MHD,NRS,ERS,BRS,RRS,DRS,GRS,is_calc,Fnu[4],Fnu[5],Fnu[6],Fnu[13],Fnu[14],Fnu[15],Fnu[16],Fnu[17],Fnu[18])

    # Interpolation
    if is_calc[1] == true
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[1],Out[1])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[4],Out[4])
    end
    if is_calc[2] == true
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[2],Out[2])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[5],Out[5])
    end
    if is_calc[3] == true
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[3],Out[3])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[6],Out[6])
    end
    if is_calc[4] == true
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[7],Out[7])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[8],Out[8])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[9],Out[9])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[13],Out[13])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[14],Out[14])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[15],Out[15])
    end
    if is_calc[5] == true
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[10],Out[10])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[11],Out[11])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[12],Out[12])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[16],Out[16])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[17],Out[17])
        Data_interpolation!(t_obs,time,nu_obs,freq,Fnu[18],Out[18])
    end

    return Out
end

# Forward Shock model - traditional Fire Ball
function Fire_Ball!(E0,G0,n0,k,epsiron_B,z,time,T,timeFS,G,R,N,E,B,DFS)
    A=n0*(3e35)^(k/2.0)
    n1=A*R0^(-k)
    t_dec=((4.0*pi*A*mp*1.6e-12*G0^2)/((3-k)*E0))^(3-k)/c
    @inbounds for i in eachindex(T)
        if i == 1
            R[1]=R0+c*T[1]*sqrt(1.0-1.0/(2.0*G0*G0))
            G[1]=G0
            n1*=(R0/R[1])^k # CSM structure
        else
            R[i] = R[i-1]+c*(T[i]-T[i-1])*sqrt(1.0-1.0/(2.0*G[i-1]^2))
            G[i] = G[i-1] 
            n1*=(R[i-1]/R[i])^k # CSM structure
            if T[i] > t_dec # Blandford-McKee phase
                GB=sqrt(G[i-1]*G[i-1]-1.0)*(T[i-1]/T[i])^(0.5*(3.0-k)) # Gamma Beta
                G[i]=sqrt(1.0+GB*GB)
            end
        end
        N[i] = (4.0*G[i]+3.0)*n1
        E[i] = (G[i]-1)*(N[i]*mp*1.6e-12)
        B[i] = sqrt(8.0*pi*epsiron_B*E[i])
        DFS[i]=(R[i]-R0)/G[i]/(12.0-4.0*k)
        if (1+z)*( T[i] - (R[i]-R0)/c ) < 1e-1*minimum(time) # minimum observed time 
            timeFS[1]=i
        end
        if sqrt(G[i]*G[i]-1.0) > 2
            timeFS[2]=i
        end
    end
    return timeFS,G,R,N,E,B,DFS
end

# Forward & Reverse Shock model - Magnetic Bullet (See Kusafuka & Asano (2024) for details)
function Magnetic_Bullet!(E0,G0,S0,D0,n0,k,epsiron_BFS,epsiron_BRS,z,time,T,timeFS,GFS,RFS,NFS,EFS,BFS,DFS,timeRS,GRS,RRS,NRS,ERS,BRS,DRS)
    A=n0*(3e35)^(k/2.0)
    n1=A*R0^(-k)
    n2=E0/(4*pi*D0*(1+S0)*G0*G0*mp*1.6e-12)/R0^2 # ejecta density
    t_acc=D0/c # acceleration timescale
    t_coast=T[end]
    t_RS=T[end]
    t_BM=T[end]
    delta=0
    delta_RS=0
    D1=D0
    G1=G0
    S1=S0

    @inbounds for i in eachindex(T)
        if 1.0 <= 8.0/3.0*G1*G1*n1/n2/S1 # RS ignition condition
            t_coast=T[max(1,i-1)]
            S1*=-1.0
        elseif delta_RS >= D0 # RS crossing condition
            t_RS=T[max(1,i-1)]
            delta_RS=0.0
        elseif delta >= D1 # Rarefaction catch-up condition
            t_BM=T[max(1,i-1)]
            delta=0.0
        end
        if i == 1
            RFS[1]=R0+c*T[1]*sqrt(1.0-1.0/(2.0*G1*G1))
            RRS[1]=R0+c*T[1]*sqrt(1.0-1.0/(1.0*G1*G1))
            GFS[1]=G1
            n1*=(R0/RFS[1])^k # CSM structure
        else
            RFS[i] = RFS[i-1]+c*(T[i]-T[i-1])*sqrt(1.0-1.0/(2.0*GFS[i-1]^2))
            RRS[i] = RRS[i-1]+c*(T[i]-T[i-1])*sqrt(1.0-1.0/(1.0*GFS[i-1]^2))
            GFS[i] = GFS[i-1] 
            n1*=(RFS[i-1]/RFS[i])^k # CSM structure
            n2=E0/(4*pi*D0*(1+S0)*G1*G0*mp*1.6e-12)/RRS[i]^2 # ejecta density
            if T[i] < t_coast # Acceleration phase
                GFS[i]=min( G0*(1+S0*T[i]/t_acc)^(1/3), G1/sqrt(1.0+2.0*G1*sqrt(n1/n2/(1+S1))) )
                S1=S0*(1+S0*T[i]/t_acc)^(-1/3)
                G1=min( G0*(1.0+S0*T[i]/t_acc)^(1.0/3.0), G0*(1+S0))
                D1+=c*(T[i]-T[i-1])*( sqrt(1.0-1.0/(2.0*GFS[i-1]^2)) - sqrt(1.0-1.0/(GFS[i-1]^2)) )
                timeRS[1]=i+1
            elseif  T[i] >= t_coast && T[i] < t_RS # Coasting phase
                GFS[i]=G1/sqrt(1.0+2.0*G1*sqrt(n1/n2/(1.0+abs(S1))))
                GFS[i]=sqrt(1.0+GFS[i]*GFS[i])
                GRS[i]=GFS[i]
                D1+=c*(T[i]-T[i-1])*( sqrt(1.0-1.0/(2.0*GFS[i-1]^2)) - sqrt(1.0-1.0/(GFS[i-1]^2)) )
                G_rel=0.5*(GRS[i]/G1+G1/GRS[i])
                NRS[i]=n2*( 4.0*G_rel + 3.0 )
                NRS[i],ERS[i]=uRS!(G_rel,abs(S1),NRS[i],ERS[i]) # find RS solution
                delta_RS+=c*(T[i]-T[i-1])*(sqrt(1.0-1.0/G1/G1)-sqrt(1.0-1.0/GRS[i]/GRS[i]))/(1.0-G1*n2/GRS[i]/NRS[i])
                DRS[i]=DRS[i-1]+c*(T[i]-T[i-1])*GRS[i]*( sqrt(1.0-1.0/(G1^2)) - sqrt(1.0-1.0/(GRS[i]^2)) )
                BRS[i] = max(sqrt(8.0*pi*epsiron_BRS*ERS[i]),sqrt(abs(S1)*4.0*pi*NRS[i]*mp*1.6e-12))
                timeRS[2]=i
            elseif T[i] >= t_RS && T[i] < t_BM # Transition phase
                GFS[i]=G1/sqrt(1.0+2.0*G1*sqrt(n1/n2/(1.0+abs(S1))))
                GFS[i]=sqrt(1.0+GFS[i]*GFS[i])
                delta+=c*(T[i]-T[i-1])*( 1.0 - sqrt(1.0-1.0/(2.0*GFS[i-1]^2)) )
            else # Blandford-McKee phase
                GB=sqrt(GFS[i-1]*GFS[i-1]-1.0)*(T[i-1]/T[i])^(0.5*(3.0-k)) # Gamma Beta
                GFS[i]=sqrt(1.0+GB*GB)
            end
        end
        NFS[i] = (4.0*GFS[i]+3.0)*n1
        EFS[i] = (GFS[i]-1.0)*(NFS[i]*mp*1.6e-12)
        BFS[i] = sqrt(8.0*pi*epsiron_BFS*EFS[i])
        DFS[i]=(RFS[i]-R0)/GFS[i]/(12.0-4.0*k)
        if (1+z)*( T[i] - (RFS[i]-R0)/c ) < 1e-1*minimum(time) # minimum observed time 
            timeFS[1]=i
        end
        if sqrt(GFS[i]*GFS[i]-1.0) > 2
            timeFS[2]=i
        end
    end
    return timeFS,GFS,RFS,NFS,EFS,BFS,DFS,timeRS,GRS,RRS,NRS,ERS,BRS,DRS
end

# find RS solution
function uRS!(G,S,n_RS,E_RS)
    A=8.0*G+10.0
    B=-(G+1.0)*(8.0*G*G+4.0*G+6.0)*S-(G-1.0)*(8.0*G*G+18.0*G+11.0)
    C=(G+1.0)*(8.0*G*G+1.0)*S*S+(G*G-1.0)*(10.0*G+6.0)*S+(G+1.0)*(G-1.0)*(G-1.0)
    D=-(G-1.0)*(G+1.0)*(G+1.0)*S*S
    a1=-2.0*B*B*B+9.0*A*B*C-27.0*D*A*A
    a4=-B*B+3.0*A*C
    a3=a1*a1+4.0*a4*a4*a4
    a2=a1+sqrt(a3+0im)
    x1=-B/A/3.0+a2^(1.0/3.0)/3.0/2.0^(1.0/3.0)/A-2^(1.0/3.0)*a4/3.0/a2^(1.0/3.0)/A
    x2=-B/A/3.0-(1-sqrt(3)im)*a2^(1.0/3.0)/6.0/2.0^(1.0/3.0)/A+(1+sqrt(3)im)*a4/3.0/2^(2.0/3.0)/a2^(1.0/3.0)/A
    x3=-B/A/3.0-(1+sqrt(3)im)*a2^(1.0/3.0)/6.0/2.0^(1.0/3.0)/A+(1-sqrt(3)im)*a4/3.0/2^(2.0/3.0)/a2^(1.0/3.0)/A
    sol=[real(sqrt(x1)),real(sqrt(x2)),real(sqrt(x3))]
    minRS=minimum(sol)
    uRS=minimum(setdiff!(sol,[minRS]))
    n_RS*=(G+(sqrt(uRS*uRS+1.0)/uRS)*sqrt(G*G-1.0))/(4.0*G+3.0)
    E_RS=(1.0-0.5*S*(G+1.0)/(uRS*uRS*G+uRS*sqrt(uRS*uRS+1.0)*sqrt(G*G-1.0)))*(G-1.0)*n_RS*mp*1.6e-12
    return n_RS,E_RS
end

# Calculate particle distribution
function particle_distribution!(epsiron_e,epsiron_B,p,f_e,B,N,G,R,gamma_m,gamma_c,PD,gamma_e) 
    Y=zeros(length(gamma_e))

    if gamma_m < gamma_c # slow cooling
        @inbounds for j in eachindex(gamma_e)
            nu_KN=G*(me*1.6e-12)/h/gamma_e[j]
            gamma_KN=sqrt(4.0*pi*(me*1.6e-12/c)*nu_KN/(3.0*e*B))
            if gamma_KN < gamma_m
                eta=0.0
            elseif gamma_KN >= gamma_m && gamma_KN < gamma_c
                eta=(gamma_KN^(3.0-p) - gamma_m^(3.0-p))/(1.0/(p-2.0)*gamma_c^(3.0-p) - gamma_m^(3.0-p)  )
            else
                eta=1.0-((3.0-p)*gamma_c*gamma_KN^(2.0-p))/(gamma_c^(3.0-p) - (p-2.0)*gamma_m^(3.0-p)  )
            end
            Y[j]=(-1.0+sqrt(1.0+4.0*(gamma_m/gamma_c)^(p-2)*eta*epsiron_e/epsiron_B))/2.0
        end
    else # fast cooling
        for j in eachindex(gamma_e)
            nu_KN=G*(me*1.6e-12)/h/gamma_e[j]
            gamma_KN=sqrt(4.0*pi*(me*1.6e-12/c)*nu_KN/(3.0*e*B))
            if gamma_KN < gamma_c
                eta=0.0
            elseif gamma_KN >= gamma_c && gamma_KN <= gamma_m
                eta=(gamma_KN - gamma_c)/((p-1.0)/(p-2.0)*gamma_m - gamma_c )
            else
                eta=1.0-(gamma_m^(p-1.0)*gamma_KN^(2.0-p))/((p-1.0)*gamma_m - (p-2.0)*gamma_c )
            end
            Y[j]=(-1.0+sqrt(1.0+4.0*1.0*eta*epsiron_e/epsiron_B))/2.0
        end
    end

    # root find for gamma_c
    gamma_syn=gamma_c
    del=1e1
    for j in eachindex(gamma_e)
        if abs(gamma_e[j]*(1+Y[j]) - gamma_syn) / gamma_syn <= del
            del = abs(gamma_e[j]*(1+Y[j]) - gamma_syn) / gamma_syn
            gamma_c=max(gamma_syn/(1.0+Y[j]),1.0)
        end
    end

    @inbounds for j in eachindex(gamma_e)
        t_syn=6.0*pi*me*1.6e-12/c/sigma_T/B^2/gamma_e[j]
        t_adi=R/sqrt(G^2-1.0)/c
        t_cool=1.0/((1+Y[j])/t_syn+1.0/t_adi)
        if gamma_e[j] < min(gamma_m,gamma_c)
            PD[j]=0.0
        else
            PD[j]=t_cool/gamma_e[j]*max(gamma_e[j],gamma_m)^(1-p)
        end
    end

    # normalization
    C=0.0
    for j in 1:length(gamma_e)-1
        C+=0.5*(PD[j]+PD[j+1])*(gamma_e[j+1]-gamma_e[j])
    end
    PD .= @. PD*f_e*N/C
    return PD, gamma_e
end

function proton_distribution!(epsiron_p,epsiron_B,p,f_p,B,N,G,R,nu_rad,j_rad,gamma_m,gamma_c,PD,gamma_p) 

    @inbounds for j in eachindex(gamma_p)
        t_pg=cool_pg(gamma_p[j],nu_rad,j_rad,R/(G*c))
        t_pp=1.0/(N*sigma_pp(gamma_p[j]*mp)*c)
        t_syn=6.0*pi*me*1.6e-12/c/sigma_T/B^2/gamma_p[j]*(mp/me)^3
        t_adi=R/sqrt(G^2-1.0)/c
        t_cool=1.0/(1.0/t_syn+1.0/t_pp+1.0/t_pg+1.0/t_adi)
        if gamma_p[j] < min(gamma_m,gamma_c)
            PD[j]=0.0
        else
            PD[j]=t_cool/gamma_p[j]*max(gamma_p[j],gamma_m)^(1-p)
        end
    end

    # normalization
    C=0.0
    for j in 1:length(gamma_p)-1
        C+=0.5*(PD[j]+PD[j+1])*(gamma_p[j+1]-gamma_p[j])
    end
    PD .= @. PD*f_p*N/C
    return PD, gamma_p
end

# Synchrotron function fitting formulae: Fouka & Ouichaoui (2013)
@inline function Bessel_F(y)
    H1=-0.97947838884478688*y-0.83333239129525072*y^(0.5)+0.15541796026816246*y^(1.0/3.0)
    H2=-4.69247165562628882e-2*y-0.70055018056462881*y^(0.5)+1.03876297841949544e-2*y^(1.0/3.0)
    # A1=pi*2^(5.0/3.0)/sqrt(3.0)/gamma(1.0/3.0)*y^(1.0/3.0)
    # A2=sqrt(pi/2.0)*sqrt(y)*exp(-y)
    A1=2.1495282415344787*y^(1.0/3.0)
    A2=1.2533141373155001*sqrt(y)*exp(-y)
    Bessel_F=A1*exp(H1)+A2*(1.0-exp(H2))
end

# Synchrotron polarization function fitting formulae: Fouka & Ouichaoui (2013)
@inline function Bessel_G(y)
    H1=-1.0010216415582440*y+0.88350305221249859*y^(0.5)-3.6240174463901829*y^(1.0/3.0)+0.57393980442916881*y^(0.25)
    H2=-0.2493940736333195*y+0.9122693061687756*y^(0.5)+1.2051408667145216*y^(1.0/3.0)-5.5227048291651126*y^(0.25)
    # A1=0.5*gamma(2.0/3.0)*(y/2.0)^(-2.0/3.0)
    # A2=sqrt(pi/2.0)*sqrt(y)*exp(-y)
    A1=0.67705896971*(y/2.0)^(-2.0/3.0)
    A2=1.25331413732*sqrt(y)*exp(-y)
    Bessel_G=y*(A1*exp(H1)+A2*(1.0-exp(H2)))
end

# Klein-Nishina correction
@inline function KN(g,f)
    q=f/4.0/g/(1-f)
    ans=2.0*q*log(q)+(1.0+2.0*q)*(1.0-q)+0.5*(4.0*g*q)^2/(1.0+4.0*g*q)*(1.0-q)
    if ans>0.0
        KN=ans
    else
        KN=0.0
    end
end

# Cross section of pair annihilation
@inline function sigma_gg(s)
    sigma_gg=3.0*sigma_T/16.0/s/s*( 4.0*(s+1.0-0.5/s)*log( sqrt(s)+sqrt(s-1) ) - 2.0*(s+1)*sqrt(1.0-1.0/s) )
end

# Cross section of pp collision
@inline function sigma_pp(E_p)
    sigma_pp=30*(0.95+0.06*log(E_p/1e9))*1e-27 # cm^2 
end

# Cooling timescale of pg interaction
@inline function cool_pg(gamma_p,nu,jnu,t_dyn)
    # Delta resonance approximation: Shigeo. S. Kimura (2022)
    sigma_pg=5e-28 # cross section cm^2
    E_pk=0.3 # peak energy GeV
    D_pk=0.2 # peak width GeV
    kappa=0.2 # inelasticity
    n=0.0 # cm^-3 GeV^-2
    @inbounds for i in 1:length(nu)-1
        if h*nu[i]/1.6e-3 >= 0.5*E_pk/gamma_p
            n+=(4.0*pi*jnu[i])/(h*nu[i])*t_dyn*(nu[i+1]-nu[i])/(h*nu[i]/1.6e-3)/(h*nu[i]/1.6e-3)
        end
    end
    cool_pg=2.0*gamma_p*gamma_p/(c*sigma_pg*kappa*D_pk*E_pk*n) # s
end

# In the case of Off axis
@inline function phi(theta_o,theta_j,theta)
    if sin(theta_o) <= 0
        phi=pi
    else
        if theta<=theta_j-theta_o
            phi=pi
        elseif theta >= theta_o - theta_j
            phi=acos(min(1.0,max(-1.0,(cos(theta_j)-cos(theta_o)*cos(theta))/(sin(theta)*sin(theta_o)))))
        else
            phi=0.0
        end
    end
end

# Calculate observed flux 
function calculate_radiation!(z,DL,time,epsiron_e,epsiron_B,p,f_e,epsiron_p,f_p,theta_jet,theta_view,nu_obs,t_obs,time_S,t_MHD,N,E,B,Radius,Delta,G,is_calc,Syn,Synp,SSC,PPnue,PPnumu,PPgamma,PGnue,PGnumu,PGgamma)
    gm_norm=(p-2.0)/(p-1.0)*epsiron_e/f_e/(me*1.6e-12)
    gc_norm=6.0*pi*me*1.6e-12/sigma_T
    gM_norm=6.0*pi*e/sigma_T

    @inbounds for i in time_S[1] : time_S[2]
        # Angle range
        mu=collect(LinRange(cos(min(pi, 3.0/G[i]+max(0, theta_view-theta_jet), theta_view+theta_jet)),cos(max(0, theta_view-theta_jet)),n_theta))

        # Electron Distribution
        gamma_M=min(sqrt(gM_norm/B[i]), e*B[i]*Delta[i]/(me*1.6e-12))
        gamma_m=max(gm_norm*(E[i]/N[i]),1.0)
        gamma_c=max(gc_norm/(B[i]*B[i]*Radius[i])*G[i],1.0)
        gamma_len=20*(1+Int(floor(log10(gamma_M))))
        if gamma_len < 1
            continue
        end
        gamma=collect(LinRange(0,log10(gamma_M),gamma_len))
        gamma_e=zeros(Float64,gamma_len)
        @inbounds for j in eachindex(gamma_e)
            gamma_e[j]=10^gamma[j]
        end
        PDe=zeros(Float64,gamma_len)
        particle_distribution!(epsiron_e,epsiron_B,p,f_e,B[i],N[i],G[i],Radius[i],gamma_m,gamma_c,PDe,gamma_e) 

        # Calculate e-synchrotron radiation
        nu=collect(LinRange(6,24,Int((24-6)/delta_nu)+1))
        nu_syn=zeros(Float64,length(nu))
        nu_syn.=10 .^nu
        j_syn=zeros(Float64,length(nu_syn))
        electron_synchrotron!(gamma_e,PDe,Delta[i],B[i],nu_syn,j_syn)
        comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_syn,j_syn,Syn) 

        # Calculate SSC radiation
        if is_calc[3] == true 
            nu=collect(LinRange(12,28,Int((28-12)/delta_nu)+1))
            nu_SSC=zeros(Float64,length(nu))
            nu_SSC.=10 .^nu
            j_SSC=zeros(Float64,length(nu_SSC))
            electron_SSC!(gamma_e,PDe,Delta[i],nu_syn,j_syn,nu_SSC,j_SSC)
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_SSC,j_SSC,SSC) 
        end

        # proton distribution
        if is_calc[2] == true || is_calc[4] == true || is_calc[5] == true 
            gamma_M=sqrt(6.0*pi*e/sigma_T/B[i]*(mp/me)^2)
            gamma_m=max((p-2.0)/(p-1.0)*epsiron_p/f_p/(mp*1.6e-12)*(E[i]/N[i]),1.0)
            gamma_c*=(mp/me)^3
            gamma_len=20*(1+Int(floor(log10(gamma_M))))
            if gamma_len < 1
                continue
            end
            gamma=collect(LinRange(0,log10(gamma_M),gamma_len))
            gamma_p=zeros(Float64,gamma_len)
            gamma_p=10 .^gamma
            PDp=zeros(Float64,gamma_len)
            proton_distribution!(epsiron_p,epsiron_B,2.0,f_p,B[i],N[i],G[i],Radius[i],nu_syn,j_syn,gamma_m,gamma_c,PDp,gamma_p) 
        end

        # Calculate p-synchrotron radiation
        if is_calc[2] == true 
            nu=collect(LinRange(6,28,Int((28-6)/delta_nu)+1))
            nu_synp=zeros(Float64,length(nu))
            nu_synp.=10 .^nu
            j_synp=zeros(Float64,length(nu_synp))
            proton_synchrotron!(gamma_p,PDp,Delta[i],B[i],nu_synp,j_synp)
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_synp,j_synp,Synp) 
        end

        # Calculate pp collision
        if is_calc[4] == true
            nu=collect(LinRange(20,35,Int((35-20)/delta_nu)+1))
            nu_pp=zeros(Float64,length(nu))
            nu_pp.=10 .^nu
            E_pp=nu_pp ./ 1.6e-12 .* h
            j_pp_gamma=zeros(Float64,length(E_pp)) # pi0 gamma emissivity : eV/cc/s/str/eV
            j_pp_nue=zeros(Float64,length(E_pp))   # nu_e emissivity      : eV/cc/s/str/eV
            j_pp_numu=zeros(Float64,length(E_pp))  # nu_mu emissivity     : eV/cc/s/str/eV
            pp_collision!(gamma_p.*mp,PDp./mp,Delta[i],N[i],B[i],E_pp,j_pp_gamma,j_pp_nue,j_pp_numu)
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_pp,j_pp_nue,PPnue) 
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_pp,j_pp_numu,PPnumu) 
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_pp,j_pp_gamma,PPgamma) 
        end

        # Calculate pg interaction 
        if is_calc[5] == true
            nu=collect(LinRange(20,35,Int((35-20)/delta_nu)+1))
            nu_pg=zeros(Float64,length(nu))
            nu_pg.=10 .^nu
            E_pg=nu_pg ./ 1.6e-12 .* h
            j_pg_gamma=zeros(Float64,length(E_pg)) # pi0 gamma emissivity : eV/cc/s/str/eV
            j_pg_se=zeros(Float64,length(E_pg))    # second e emissivity  : eV/cc/s/str/eV
            j_pg_nue=zeros(Float64,length(E_pg))   # nu_e emissivity      : eV/cc/s/str/eV
            j_pg_numu=zeros(Float64,length(E_pg))  # nu_mu emissivity     : eV/cc/s/str/eV
            pg_interaction!(gamma_p.*mp,PDp./mp,nu_syn,j_syn,Delta[i],N[i],B[i],E_pg,j_pg_gamma,j_pg_se,j_pg_nue,j_pg_numu,t_MHD[i+1]-t_MHD[i])
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_pg,j_pg_nue,PGnue) 
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_pg,j_pg_numu,PGnumu) 
            comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,t_MHD[i-1],t_MHD[i],Radius[i-1],Radius[i],G[i],Delta[i],mu,nu_pg,j_pg_gamma,PGgamma) 
        end
    end
    return Syn,Synp,SSC,PPnue,PPnumu,PPgamma,PGnue,PGnumu,PGgamma
end

function electron_synchrotron!(gamma_e,PDe,dr,B,nu_syn,j_syn)
    nu_norm=3.0*e*B*sin(alpha)/(4*pi*me*1.6e-12/c)
    F_syn=sqrt(3.0)*(e*e*e*B*sin(alpha))/(me*1.6e-12)
    @inbounds for j in eachindex(nu_syn) 
        tau=1e-10 
        tau_norm=dr/(8.0*pi*(me*1.6e-12/c/c)*nu_syn[j]*nu_syn[j])
        # comoving synchrotron emissivity
        @inbounds for k in 1:length(gamma_e)-1
            nu_c=(nu_norm*gamma_e[k]*gamma_e[k])
            P=F_syn*Bessel_F(nu_syn[j]/nu_c) 
            j_syn[j]+=0.5*P*(PDe[k]+PDe[k+1])*(gamma_e[k+1]-gamma_e[k])
            tau+=P*(2.0*(PDe[k]/gamma_e[k]+PDe[k+1]/gamma_e[k+1])/2.0*(gamma_e[k+1]-gamma_e[k]) - (PDe[k+1]-PDe[k]))*tau_norm
        end
        # SSA
        j_syn[j]*=(1.0-exp(-1.0*tau))/tau/(4.0*pi)
    end
    return j_syn
end

function proton_synchrotron!(gamma_p,PDp,dr,B,nu_synp,j_synp)
    nu_norm=3.0*e*B*sin(alpha)/(4*pi*mp*1.6e-12/c)
    F_syn=sqrt(3.0)*(e*e*e*B*sin(alpha))/(mp*1.6e-12)
    @inbounds for j in eachindex(nu_synp) 
        tau=1e-10 
        tau_norm=dr/(8.0*pi*(mp*1.6e-12/c/c)*nu_synp[j]*nu_synp[j])
        # comoving synchrotron emissivity
        @inbounds for k in 1:length(gamma_p)-1
            nu_c=(nu_norm*gamma_p[k]*gamma_p[k])
            P=F_syn*Bessel_F(nu_synp[j]/nu_c) 
            j_synp[j]+=0.5*P*(PDp[k]+PDp[k+1])*(gamma_p[k+1]-gamma_p[k])
            tau+=P*(2.0*(PDp[k]/gamma_p[k]+PDp[k+1]/gamma_p[k+1])/2.0*(gamma_p[k+1]-gamma_p[k]) - (PDp[k+1]-PDp[k]))*tau_norm
        end
        # SSA
        j_synp[j]*=(1.0-exp(-1.0*tau))/tau/(4.0*pi)
    end
    return j_synp
end

function electron_SSC!(gamma_e,PDe,dr,nu_syn,j_syn,nu_SSC,j_SSC)
    @inbounds for j in eachindex(nu_SSC)
        taug=1e-10
        # comoving SSC flux
        @inbounds for k in 1:length(gamma_e)-1
            f=h*nu_SSC[j]/gamma_e[k]/(me*1.6e-12)
            n_ssc=0.0
            @inbounds for seed in 1:length(nu_syn)-1
                g=gamma_e[k]*h*nu_syn[seed]/(me*1.6e-12)
                if g/gamma_e[k]^2 < f && f <= 4.0*g/(1.0+4.0*g)
                    n_ph=(4.0*pi*j_syn[seed])/(h*nu_syn[seed])*(dr/c)
                    n_ssc+=(3.0*sigma_T*c*n_ph/4.0/gamma_e[k]^2/nu_syn[seed])*KN(g,f)*(nu_syn[seed+1]-nu_syn[seed])
                    s=(h*nu_SSC[j])*(h*nu_syn[seed])/(2.0*me*1.6e-12*me*1.6e-12)
                    if s > 1.0
                        taug+=sigma_gg(s)*n_ph*delta_nu*nu_syn[seed]*dr
                    end
                end
            end
            j_SSC[j]+=h*nu_SSC[j]*n_ssc*PDe[k]*(gamma_e[k+1]-gamma_e[k])
        end
        # pair annihilation
        j_SSC[j]*=(1.0-exp(-1.0*taug))/taug/(4.0*pi)
    end
    return j_SSC
end

function pp_collision!(E_p,PDp,dr,N,B,E_pp,j_gamma,j_nue,j_numu)
    # calculate pion injection rate
    E_pion=zeros(Float64,length(E_p))
    E_pion=E_p # unit eV
    j_pion=zeros(Float64,length(E_pion)) # injection rate : number/cc/s/eV
    @inbounds for k in eachindex(E_pion)
        @inbounds for j in 1:length(E_p)-1 
            if mpi < E_pion[k] && E_pion[k] < E_p[j] 
                j_pion[k]+=c*N*PDp[j]*sigma_pp(E_p[j])*pp_spec(E_pion[k],E_p[j])*(E_p[j+1]-E_p[j])/E_p[j] 
            end
        end
    end

    # comoving pi0 gamma emissivity
    @inbounds for j in eachindex(E_pp) 
        @inbounds for k in 1:length(E_pion)-1
            if E_pp[j] > E_pion[k]
                continue
            end
            j_gamma[j]+=E_pp[j]*2.0*j_pion[k]*(E_pion[k+1]-E_pion[k])/E_pion[k] # comoving pi0 gamma emissivity
            j_nue[j]+=E_pp[j]*2.0*f_nue(E_pp[j]/E_pion[k])*j_pion[k]*(E_pion[k+1]-E_pion[k])/E_pion[k] # comoving nu_e emissivity
            j_numu[j]+=E_pp[j]*2.0*f_numu(E_pp[j]/E_pion[k])*j_pion[k]*(E_pion[k+1]-E_pion[k])/E_pion[k] # comoving nu_mu emissivity
            if E_pp[j] < lambda*E_pion[k]
                j_numu[j]+=E_pp[j]*2.0/lambda*j_pion[k]*(E_pion[k+1]-E_pion[k])/E_pion[k] # direct decay
            end
        end
        j_gamma[j]*=h/(4.0*pi)/3.0 # eV/cm^2/s/eV -> erg/cm^2/s/str/Hz
        j_nue[j]*=h/(4.0*pi)/3.0   # eV/cm^2/s/eV -> erg/cm^2/s/str/Hz
        j_numu[j]*=h/(4.0*pi)/3.0  # eV/cm^2/s/eV -> erg/cm^2/s/str/Hz
    end
    return j_gamma,j_nue,j_numu
end

function pg_interaction!(E_p,PDp,nu_ph,j_ph,dr,N,B,E_pg,j_gamma,j_nue,j_numu)
    E_ph=nu_ph.*h./1.6e-12

    @inbounds for j in eachindex(E_pg) 
        @inbounds for k in 1:length(E_p)-1
            x=E_pg[j]/E_p[k]
            @inbounds for i in 1:length(E_ph)-1
                eta=4.0*E_ph[i]*E_p[k]/mp/mp
                if eta < 2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
                    continue
                end
                j_gamma[j]+=E_pg[j]*PDp[k]*(4*pi*j_ph[i]/h/nu_ph[i]*dr/c)*phi_gamma(eta,x)*(E_p[k+1]-E_p[k])/E_p[k]*(nu_ph[i+1]-nu_ph[i]) # comoving pi0 gamma emissivity
                j_nue[j]+=E_pg[j]*PDp[k]*(4*pi*j_ph[i]/h/nu_ph[i]*dr/c)*phi_nue(eta,x)*(E_p[k+1]-E_p[k])/E_p[k]*(nu_ph[i+1]-nu_ph[i]) # comoving nu_e emissivity
                j_numu[j]+=E_pg[j]*PDp[k]*(4*pi*j_ph[i]/h/nu_ph[i]*dr/c)*phi_numu(eta,x)*(E_p[k+1]-E_p[k])/E_p[k]*(nu_ph[i+1]-nu_ph[i]) # comoving nu_mu emissivity
            end
        end
        j_gamma[j]*=h/(4.0*pi) # eV/cm^2/s/eV -> erg/cm^2/s/str/Hz
        j_nue[j]*=h/(4.0*pi)   # eV/cm^2/s/eV -> erg/cm^2/s/str/Hz
        j_numu[j]*=h/(4.0*pi)  # eV/cm^2/s/eV -> erg/cm^2/s/str/Hz
    end
    return j_gamma,j_nue,j_numu
end

@inline function pp_spec(E_pion,E_proton) # Kelner & Aharonian (2006) SIBYLL code
    x=E_pion/E_proton
    L=log(E_proton/1e12)
    a=3.67+0.83*L+0.075*L*L
    B_pi=a+0.25
    alpha=0.98/sqrt(a)
    r=2.6/sqrt(a)
    pion_spec=4.0*alpha*B_pi*x^(alpha-1.0)*((1.0-x^alpha)/(1.0+r*x^alpha*(1.0-x^alpha)))^4.0*(1.0/(1.0-x^alpha)+r*(1.0-2.0*x^alpha)/(1.0+r*x^alpha*(1.0-x^alpha)))*(1.0-mpi/E_pion)^0.5
end

@inline function f_numu(x)
    r=1.0-lambda
    if x > r
        f_numu=(3.0-2.0*r)*(9.0*x^2-6.0*log(x)-4.0*x^3-5.0)/9.0/(1.0-r)^2
    else
        h1=(3.0-2.0*r)*(9.0*r^2-6.0*log(r)-4.0*r^3-5.0)/9.0/(1.0-r)^2
        h2=(1.0+2.0*r)*(r-x)*(9.0*(r+x)-4.0*(r^2+r*x+x^2))/9.0/r^2
        f_numu=h1+h2
    end
end

@inline function f_nue(x)
    r=1.0-lambda
    if x > r
        f_nue=2.0*((1.0-x)*(6.0*(1.0-x)^2+r*(5.0+5.0*x-4.0*x^2))+6.0*r*log(x))/3.0/(1.0-r)^2
    else
        h1=2.0*((1.0-r)*(6.0-7.0*r+11.0*r^2-4.0*r^3)+6.0*r*log(r))/3.0/(1.0-r)^2
        h2=2.0*(r-x)*(7.0*r^2+7.0*x*r-2.0*x^2-4.0*(r^3+r^2*x+r*x^2))/3.0/r^2
        f_nue=h1+h2
    end
end

const phi_0_table=[
    [1.1 0.0768 0.544 2.86e-19];
    [1.2 0.106 0.540 2.24e-18];
    [1.3 0.182 0.750 5.61e-18];
    [1.4 0.201 0.791 1.02e-17];
    [1.5 0.219 0.788 1.60e-17];
    [1.6 0.216 0.831 2.23e-17];
    [1.7 0.233 0.839 3.10e-17];
    [1.8 0.233 0.825 4.07e-17];
    [1.9 0.248 0.805 5.30e-17];
    [2.0 0.244 0.779 6.74e-17];
    [3.0 0.188 1.23 1.51e-16];
    [4.0 0.131 1.82 1.24e-16];
    [5.0 0.120 2.05 1.37e-16];
    [6.0 0.107 2.19 1.62e-16];
    [7.0 0.102 2.23 1.71e-16];
    [8.0 0.0932 2.29 1.78e-16];
    [9.0 0.0838 2.37 1.84e-16];
    [10.0 0.0761 2.43 1.93e-16];
    [20.0 0.107 2.27 4.74e-16];
    [30.0 0.0928 2.33 7.70e-16];
    [40.0 0.0772 2.42 1.06e-15];
    [100.0 0.0479 2.59 2.73e-15]
]

const phi_plus_table=[ # positron & anti_numu & numu & nue
    [1.1 0.367 3.12 8.09e-19 0.365 3.09 8.09e-19 1e-10 1e-10 1.08e-18 0.768 2.49 9.43e-19];
    [1.2 0.282 2.96 7.70e-18 0.287 2.96 7.70e-18 0.0778 0.306 9.91e-18 0.569 2.35 9.22e-18];
    [1.3 0.260 2.83 2.05e-17 0.250 2.89 1.99e-17 0.242 0.792 2.47e-17 0.491 2.41 2.35e-17];
    [1.4 0.239 2.76 3.66e-17 0.238 2.76 3.62e-17 0.377 1.09 4.43e-17 0.395 2.45 4.20e-17];
    [1.5 0.224 2.69 5.48e-17 0.220 2.71 5.39e-17 0.440 1.06 6.70e-17 0.31 2.45 6.26e-17];
    [1.6 0.207 2.66 7.39e-17 0.206 2.67 7.39e-17 0.450 0.953 9.04e-17 0.323 2.43 8.57e-17];
    [1.7 0.198 2.62 9.52e-17 0.197 2.62 9.48e-17 0.461 0.956 1.18e-16 0.305 2.40 1.13e-16];
    [1.8 0.193 2.56 1.20e-16 0.193 2.56 1.20e-16 0.451 0.922 1.32e-16 0.285 2.39 1.39e-16];
    [1.9 0.187 2.52 1.47e-16 0.187 2.52 1.47e-16 0.464 0.912 1.77e-16 0.270 2.37 1.70e-16];
    [2.0 0.181 2.49 1.75e-16 0.178 2.51 1.74e-16 0.446 0.940 2.11e-16 0.259 2.35 2.05e-16];
    [3.0 0.122 2.48 3.31e-16 0.123 2.48 3.38e-16 0.366 1.49 3.83e-16 0.158 2.42 3.81e-16];
    [4.0 0.106 2.50 4.16e-16 0.106 2.56 5.17e-16 0.249 2.03 5.09e-16 0.129 2.46 4.74e-16];
    [5.0 0.0983 2.46 5.57e-16 0.0944 2.57 7.61e-16 0.204 2.18 7.26e-16 0.113 2.45 6.30e-16];
    [6.0 0.0875 2.46 6.78e-16 0.0829 2.58 9.57e-16 0.174 2.24 9.26e-16 0.0996 2.46 7.65e-16];
    [7.0 0.0830 2.44 7.65e-16 0.0801 2.54 1.11e-15 0.156 2.28 1.07e-15 0.0921 2.46 8.61e-16];
    [8.0 0.0783 2.44 8.52e-16 0.0752 2.53 1.25e-15 0.140 2.32 1.19e-15 0.0861 2.45 9.61e-16];
    [9.0 0.0735 2.45 9.17e-16 0.0680 2.56 1.36e-15 0.121 2.39 1.29e-15 0.0800 2.47 1.03e-15];
    [10.0 0.0644 2.50 9.57e-16 0.0615 2.60 1.46e-15 0.107 2.46 1.40e-15 0.0723 2.51 1.10e-15];
    [30.0 0.0333 2.77 3.07e-15 0.0361 2.78 5.87e-15 0.0705 2.53 5.65e-15 0.0411 2.70 3.55e-15];
    [100.0 0.0224 2.86 1.58e-14 0.0228 2.88 3.10e-14 0.0463 2.62 3.01e-14 0.0283 2.77 1.86e-14]
]

const phi_minus_table=[ # electron & anti-nue
    [3.0 0.658 3.09 6.43e-19 0.985 2.63 6.61e-19];
    [4.0 0.348 2.81 9.91e-18 0.378 2.98 9.74e-18];
    [5.0 0.286 2.39 1.24e-16 0.31 2.31 1.34e-16];
    [6.0 0.256 2.27 2.67e-16 0.327 2.11 2.91e-16];
    [7.0 0.258 2.13 3.50e-16 0.308 2.03 3.81e-16];
    [8.0 0.220 2.20 4.03e-16 0.292 1.98 4.48e-16];
    [9.0 0.217 2.13 4.48e-16 0.260 2.02 4.83e-16];
    [10.0 0.192 2.19 4.78e-16 0.233 2.07 5.13e-16];
    [30.0 0.125 2.27 1.64e-15 0.135 2.24 1.75e-15];
    [100.0 0.0507 2.63 4.52e-15 0.0770 2.40 5.48e-15]
]

@inline function pg_loglinear_interpolation(table,pos_S,pos_D,pos_B,eta)
    B=0.0;S=0.0;D=0.0
    if eta < table[1,1]
        B=10^(log10(table[2,pos_B])+(log10(table[2,pos_B])-log10(table[1,pos_B]))/(log10(table[2,1])-log10(table[1,1]))*(log10(eta)-log10(table[2,1])))
        S=10^(log10(table[2,pos_S])+(log10(table[2,pos_S])-log10(table[1,pos_S]))/(log10(table[2,1])-log10(table[1,1]))*(log10(eta)-log10(table[2,1])))
        D=10^(log10(table[2,pos_D])+(log10(table[2,pos_D])-log10(table[1,pos_D]))/(log10(table[2,1])-log10(table[1,1]))*(log10(eta)-log10(table[2,1])))
    elseif eta >= table[end,1]
        B=10^(log10(table[end-1,pos_B])+(log10(table[end,pos_B])-log10(table[end-1,pos_B]))/(log10(table[end,1])-log10(table[end-1,1]))*(log10(eta)-log10(table[end-1,1])))
        S=10^(log10(table[end-1,pos_S])+(log10(table[end,pos_S])-log10(table[end-1,pos_S]))/(log10(table[end,1])-log10(table[end-1,1]))*(log10(eta)-log10(table[end-1,1])))
        D=10^(log10(table[end-1,pos_D])+(log10(table[end,pos_D])-log10(table[end-1,pos_D]))/(log10(table[end,1])-log10(table[end-1,1]))*(log10(eta)-log10(table[end-1,1])))
    else
        for i in 1:size(table,1)-1
            if table[i,1] <= eta && eta < table[i+1,1]
                B=10^(log10(table[i,pos_B])+(log10(table[i+1,pos_B])-log10(table[i,pos_B]))/(log10(table[i+1,1])-log10(table[i,1]))*(log10(eta)-log10(table[i,1])))
                S=10^(log10(table[i,pos_S])+(log10(table[i+1,pos_S])-log10(table[i,pos_S]))/(log10(table[i+1,1])-log10(table[i,1]))*(log10(eta)-log10(table[i,1])))
                D=10^(log10(table[i,pos_D])+(log10(table[i+1,pos_D])-log10(table[i,pos_D]))/(log10(table[i+1,1])-log10(table[i,1]))*(log10(eta)-log10(table[i,1])))
            end
        end
    end
    return S,D,B
end

@inline function phi_gamma(eta,x) # Kelner & Aharonian (2009)
    eta_0=2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
    r=mpi/mp
    x_plus=0.5*(eta+r*r+sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)
    x_minus=0.5*(eta+r*r-sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)
    psi=2.5+0.4*log(eta/eta_0)

    s_gamma,delta_gamma,B_gamma=pg_loglinear_interpolation(phi_0_table,2,3,4,eta/eta_0)

    if x <= x_minus 
        phi_gamma=B_gamma*(log(2.0))^psi
    elseif x_minus < x && x < x_plus
        y=(x-x_minus)/(x_plus-x_minus)
        phi_gamma=B_gamma*exp(-s_gamma*(log(x/x_minus))^delta_gamma)*(log(2.0/(1.0+y*y)))^psi
    else
        phi_gamma=0.0
    end
end

@inline function phi_nue(eta,x) # Kelner & Aharonian (2009)
    eta_0=2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
    r=mpi/mp
    x_plus=0.5*(eta+r*r+sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)
    x_minus=0.5*(eta+r*r-sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)/4.0
    psi=2.5+1.4*log(eta/eta_0)

    s_nue,delta_nue,B_nue=pg_loglinear_interpolation(phi_plus_table,11,12,13,eta/eta_0)

    if x <= x_minus 
        phi_nue=B_nue*(log(2.0))^psi
    elseif x_minus < x && x < x_plus
        y=(x-x_minus)/(x_plus-x_minus)
        phi_nue=B_nue*exp(-s_nue*(log(x/x_minus))^delta_nue)*(log(2.0/(1.0+y*y)))^psi
    else
        phi_nue=0.0
    end
end

@inline function phi_numu(eta,x) # Kelner & Aharonian (2009)
    eta_0=2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
    r=mpi/mp
    x_plus=0.5*(eta+r*r+sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)
    if eta/eta_0 <= 2.14
        x_plus=0.427*x_plus
    elseif 2.14 < eta/eta_0 && eta/eta_0 <= 10
        x_plus=(0.427+0.0729*(eta/eta_0-2.14))*x_plus
    end
    x_minus=0.5*(eta+r*r-sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)*0.427
    psi=2.5+1.4*log(eta/eta_0)

    s_numu,delta_numu,B_numu=pg_loglinear_interpolation(phi_plus_table,8,9,10,eta/eta_0)

    if x <= x_minus 
        phi_numu=B_numu*(log(2.0))^psi
    elseif x_minus < x && x <= x_plus
        y=(x-x_minus)/(x_plus-x_minus)
        phi_numu=B_numu*exp(-s_numu*(log(x/x_minus))^delta_numu)*(log(2.0/(1.0+y*y)))^psi
    else
        phi_numu=0.0
    end
end

@inline function phi_anti_numu(eta,x) # Kelner & Aharonian (2009)
    eta_0=2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
    r=mpi/mp
    x_plus=0.5*(eta+r*r+sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)
    x_minus=0.5*(eta+r*r-sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)/4.0
    psi=2.5+1.4*log(eta/eta_0)

    s_numu,delta_numu,B_numu=pg_loglinear_interpolation(phi_plus_table,5,6,7,eta/eta_0)

    if x <= x_minus 
        phi_anti_numu=B_numu*(log(2.0))^psi
    elseif x_minus < x && x < x_plus
        y=(x-x_minus)/(x_plus-x_minus)
        phi_anti_numu=B_numu*exp(-s_numu*(log(x/x_minus))^delta_numu)*(log(2.0/(1.0+y*y)))^psi
    else
        phi_anti_numu=0.0
    end
end

@inline function phi_p(eta,x) # Kelner & Aharonian (2009)
    eta_0=2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
    r=mpi/mp
    x_plus=0.5*(eta+r*r+sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)
    x_minus=0.5*(eta+r*r-sqrt((eta-r*r-2.0*r)*(eta-r*r+2.0*r)))/(1.0+eta)/4.0
    psi=2.5+1.4*log(eta/eta_0)

    s_p,delta_p,B_p=pg_loglinear_interpolation(phi_plus_table,2,3,4,eta/eta_0)

    if x <= x_minus 
        phi_p=B_p*(log(2.0))^psi
    elseif x_minus < x && x < x_plus
        y=(x-x_minus)/(x_plus-x_minus)
        phi_p=B_p*exp(-s_p*(log(x/x_minus))^delta_p)*(log(2.0/(1.0+y*y)))^psi
    else
        phi_p=0.0
    end
end

@inline function phi_e(eta,x) # Kelner & Aharonian (2009)
    eta_0=2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
    r=mpi/mp
    if eta < 4.0*r*(1.0+r)
        return 0.0
    end
    x_plus=0.5*(eta-2.0*r+sqrt(eta*(eta-4.0*r*(1.0+r))))/(1.0+eta)
    x_minus=0.5*(eta-2.0*r-sqrt(eta*(eta-4.0*r*(1.0+r))))/(1.0+eta)/2.0
    if eta/eta_0 > 4.0 
        psi=6.0*(1.0-exp(1.5*(4.0-eta/eta_0)))
    else
        psi=0.0
    end

    s_e,delta_e,B_e=pg_loglinear_interpolation(phi_minus_table,2,3,4,eta/eta_0)

    if x <= x_minus 
        phi_e=B_e*(log(2.0))^psi
    elseif x_minus < x && x <= x_plus
        y=(x-x_minus)/(x_plus-x_minus)
        phi_e=B_e*exp(-s_e*(log(x/x_minus))^delta_e)*(log(2.0/(1.0+y*y)))^psi
    else
        phi_e=0.0
    end
end

@inline function phi_anti_nue(eta,x) # Kelner & Aharonian (2009)
    eta_0=2.0*(mpi/mp)+(mpi/mp)*(mpi/mp)
    r=mpi/mp
    if eta < 4.0*r*(1.0+r)
        return 0.0
    end
    x_plus=0.5*(eta-2.0*r+sqrt(eta*(eta-4.0*r*(1.0+r))))/(1.0+eta)
    x_minus=0.5*(eta-2.0*r-sqrt(eta*(eta-4.0*r*(1.0+r))))/(1.0+eta)/2.0
    if eta/eta_0 > 4.0 
        psi=6.0*(1.0-exp(1.5*(4.0-eta/eta_0)))
    else
        psi=0.0
    end

    s_nue,delta_nue,B_nue=pg_loglinear_interpolation(phi_minus_table,5,6,7,eta/eta_0)

    if x <= x_minus 
        phi_anti_nue=B_nue*(log(2.0))^psi
    elseif x_minus < x && x <= x_plus
        y=(x-x_minus)/(x_plus-x_minus)
        phi_anti_nue=B_nue*exp(-s_nue*(log(x/x_minus))^delta_nue)*(log(2.0/(1.0+y*y)))^psi
    else
        phi_anti_nue=0.0
    end
end

@inline function comoving_to_observed!(z,DL,time,theta_jet,theta_view,nu_obs,t_obs,prev,now,dRadius,Radius,G,dr,mu,nu_rad,j_nu,Flux)
    GB=sqrt(G*G-1.0)
    cont=true
    @inbounds for k in eachindex(mu)
        t_prev=(1+z)*( prev - (dRadius-R0)*mu[k]/c)
        t_now=(1+z)*( now - (Radius-R0)*mu[k]/c)
        @inbounds for t in eachindex(time)
            pos_t=Int(round((log10(time[t]) - log10(t_obs[1]))/delta_t, RoundNearestTiesUp)) + 1
            if t_now < t_obs[pos_t-2] || t_prev > t_obs[pos_t+2]
                cont=false
            else
                cont=true
                break
            end
        end
        if cont==false
            continue
        end
        if t_prev <= 0.0
            continue
        end
        delta_phi=phi(theta_view,theta_jet,acos(mu[k]))
        D=1.0/G/(1.0-GB/G*mu[k])
        F_norm=(1+z)/(DL*DL)*dr*D*D*D*(2.0*delta_phi)*Radius*Radius*abs(mu[2]-mu[1])
        pos_t=Int(round((log10(t_now) - log10(t_obs[1]))/delta_t, RoundNearestTiesUp)) + 1
        pos_dt=Int(round((log10(t_prev) - log10(t_obs[1]))/delta_t, RoundNearestTiesUp)) + 1
        @inbounds for j in eachindex(nu_rad)
            F_obs=F_norm*j_nu[j]
            pos_nu=Int(round((log10(nu_rad[j]*D/(1+z)) - log10(nu_obs[1]))/delta_nu, RoundNearestTiesUp)) + 1
            if pos_nu>=1 && pos_nu <= length(nu_obs) && pos_dt>=1 && pos_t <= length(t_obs) && pos_dt<=pos_t
                if pos_dt == pos_t
                    Flux[pos_t,pos_nu]+=(F_obs*(log10(t_now)-log10(t_prev))/delta_t)
                else
                    for l in pos_dt:pos_t
                        if l == pos_dt
                            Flux[pos_dt,pos_nu]+=(F_obs*(log10(t_obs[pos_dt])+delta_t/2.0-log10(t_prev))/delta_t)
                        elseif l == pos_t
                            Flux[pos_t,pos_nu]+=(F_obs*(log10(t_now)-log10(t_obs[pos_t])+delta_t/2.0)/delta_t)
                        else
                            Flux[l,pos_nu]+=F_obs
                        end
                    end
                end
            end
        end
    end
    return Flux
end

function Data_interpolation!(t_obs,t_data,nu_obs,nu_data,F_model,F_data)
    @inbounds for i in eachindex(t_data)
        pos_t=Int(round((log10(t_data[i]) - log10(t_obs[1]))/delta_t, RoundDown)) + 1
        @inbounds for j in eachindex(nu_data)
            pos_nu=Int(round((log10(nu_data[j]) - log10(nu_obs[1]))/delta_nu, RoundDown)) + 1
            if pos_t>=1 && pos_t < length(t_obs) && pos_nu>=1 && pos_nu < length(nu_obs) 
                F_nu=10^(log10(F_model[pos_t,pos_nu+1])+(log10(F_model[pos_t+1,pos_nu+1])-log10(F_model[pos_t,pos_nu+1]))*(log10(t_data[i])-log10(t_obs[pos_t]))/delta_t)
                F_dnu=10^(log10(F_model[pos_t,pos_nu])+(log10(F_model[pos_t+1,pos_nu])-log10(F_model[pos_t,pos_nu]))*(log10(t_data[i])-log10(t_obs[pos_t]))/delta_t)
                F_data[i,j]=10^(log10(F_dnu)+(log10(F_nu)-log10(F_dnu))*(log10(nu_data[j])-log10(nu_obs[pos_nu]))/delta_nu)
            end
            if F_data[i,j] != F_data[i,j]
                F_data[i,j]=0.0
            end
        end
    end
    return F_data
end
end
#=====================================================================================================#
##### Input Arguments ######

# Input Distance
# z=1.4607 # redshift
# DL=3.3020e28 # luminosity distance [cm]

# Input Observed Time
# t_in=10 .^ range(log10(1e1), log10(1e6), length=100)

# Input Observed Frequency
# nu_in=Float64[
#     1e9,
#     1e14,
#     1e17
# ]

# Input parameters
# Input=Float64[
#     6e54   ,    # 1 initial energy E0 [erg]
#     12    ,    # 2 initial Lorentz factor G0 at R0
#     5e0   ,    # 3 initial magnetization S0
#     1400*c  ,    # 4 initial radial width D0 [cm], which should be smaller than the initial radius R0
#     2e-1    ,    # 5 CSM number density n0 at R0 [/cc]
#     2      ,    # 6 density slope k (ISM:0 --- 2:wind)
#     0.1    ,    # 7 energy fraction of accelerated electrons epsiron_e,FS
#     0.001   ,    # 8 energy fraction of Weibel induced magnetic field epsiron_B,FS
#     2.5    ,    # 9 PD spectral index 2 < p_FS < 3
#     0.4    ,    # 10 number fraction of accelerated electrons f_e,FS
#     0.025    ,    # 11 opening angle theta_j
#     0.0    ,    # 12 viewing angle theta_o
#     0.1    ,    # 13 energy fraction of accelerated electrons epsiron_e,RS
#     0.001   ,    # 14 energy fraction of Weibel induced magnetic field epsiron_B,RS
#     2.01    ,    # 15 PD spectral index 2 < p_RS < 3
#     1.0    ,    # 16 number fraction of accelerated electrons f_e,RS
#     0.1    ,    # 17 energy fraction of accelerated protons epsiron_p,FS
#     0.01        # 18 number fraction of accelerated protons f_p,FS
#     0.1    ,    # 19 energy fraction of accelerated protons epsiron_p,RS
#     0.01        # 20 number fraction of accelerated protons f_p,RS
# ]

# Output leptonic & hadronic radiation flux 
# Output=[
#     zeros(Float64,length(t_in),length(nu_in)),    # 1 e-synchrotron    FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 2 p-synchrotron    FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 3 e-SSC            FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 4 e-synchrotron    RS
#     zeros(Float64,length(t_in),length(nu_in)),    # 5 p-synchrotron    RS
#     zeros(Float64,length(t_in),length(nu_in)),    # 6 e-SSC            RS
#     zeros(Float64,length(t_in),length(nu_in)),    # 7 pp e-neutrino    FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 8 pp mu-neutrino   FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 9 pp pi0 gamma     FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 10 pg e-neutrino   FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 11 pg mu-neutrino  FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 12 pg pi0 gamma    FS
#     zeros(Float64,length(t_in),length(nu_in)),    # 13 pp e-neutrino   RS
#     zeros(Float64,length(t_in),length(nu_in)),    # 14 pp mu-neutrino  RS
#     zeros(Float64,length(t_in),length(nu_in)),    # 15 pp pi0 gamma    RS
#     zeros(Float64,length(t_in),length(nu_in)),    # 16 pg e-neutrino   RS
#     zeros(Float64,length(t_in),length(nu_in)),    # 17 pg mu-neutrino  RS
#     zeros(Float64,length(t_in),length(nu_in))     # 18 pg pi0 gamma    RS
# ]

# Which is calculated?
# is_calc=[
#     true  , # 1 e-syn
#     false , # 2 p-syn
#     false , # 3 e-SSC
#     false , # 4 pp
#     false   # 5 pg
# ]

# Model calculation
# Magnetic_Bullet_Afterglow!(z,DL,t_in,nu_in,Input,Output,is_calc)