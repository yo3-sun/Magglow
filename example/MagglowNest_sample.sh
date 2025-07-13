#!/bin/bash

# grb110213　=====================
# Output file version & memo.
version='sample'
comment="Only 7 Parameter estimation for toy event"

# distance & frequency
redshift=1.0
Luminosity_distance='2.1e28' # cm
freq=(
    "1.0e9"   # Radio
    "4.80e14" # r-band
    "2.42e18" # X-ray
)

# Observed data:
data_file=(
    "../data/toyevent_1e+09_Hz_mJy.txt" 
    "../data/toyevent_4.8e+14_Hz_mJy.txt" 
    "../data/toyevent_2.42e+18_Hz_mJy.txt" 
)
color=("g" "r" "k")

# Define P:["estimate(0) or fix(1)", "fixed(or init) value", "bound_min", "bound_max", "0"]
P_1=("0"    "54"    "51"    "55"    "54")   # 1 initial energy log(E0 [erg])
P_2=("0"    "2.0"   "1.0"   "3.0"   "2.0")  # 2 initial Lorentz factor log(G0 at R0)
P_3=("0"    "0.7"   "-1.0"  "2.0"   "0.7")  # 3 initial magnetization log(S0)
P_4=("0"    "2.5"   "1.0"   "3.0"   "2.5")  # 4 initial radial width log(D0 [cm]), which should be smaller than the initial radius R0
P_5=("1"    "0.0"   "-2.0"  "2.0"   "0.0")  # 5 CSM number density log(n0 at R0 [/cc])
P_6=("1"    "0.0"   "0.0"   "2.0"   "0.0")  # 6 density slope k (ISM:0 --- 2:wind)
P_7=("0"    "-1.0"  "-3.0"  "-0.3"  "-1.0") # 7 energy fraction of accelerated electrons log(epsilon_e)
P_8=("0"    "-3.0"  "-5.0"  "-0.3"  "-3.0") # 8 energy fraction of Weibel induced magnetic field log(epsilon_B)
P_9=("1"    "2.2"   "2.01"  "2.99"  "2.2")  # 9 PD spectral index 2 < p_FS < 3
P_10=("1"   "0.0"   "-2.0"  "0.0"   "0.0")  # 10 number fraction of accelerated electrons log(f_e)
P_11=("1"   "0.1"   "0.001" "1.5"   "0.1")  # 11 opening angle theta_j
P_12=("1"   "0.0"    "0.0"  "1.5"   "0.0")  # 12 Ratio of viewing angle and opening angle beta = theta_o/theta_j
P_13=("0"   "-1.0"  "-3.0"  "-0.3"  "-1.0") # 13 energy fraction of accelerated electrons log(epsiron_e,RS)
P_14=("1"   "-3.0"  "-5.0"  "-0.3"  "-3.0") # 14 energy fraction of Weibel induced magnetic field log(epsiron_B,RS)
P_15=("1"   "2.2"   "2.01"  "2.99"  "2.2")  # 15 PD spectral index 2 < p_RS < 3
P_16=("1"   "0.0"   "-2.0"  "0.0"   "0.0")  # 16 number fraction of accelerated electrons log(f_e,RS)
P_17=("1"   "-1.0"  "-3.0"  "0.0"   "-1.0") # 17 energy fraction of accelerated protons log(epsiron_p,FS)
P_18=("1"   "-2.0"  "-3.0"  "0.0"   "-2.0") # 18 number fraction of accelerated protons log(f_p,FS)
P_19=("1"   "-1.0"  "-3.0"  "0.0"   "-1.0") # 19 energy fraction of accelerated protons log(epsiron_p,RS)
P_20=("1"   "-2.0"  "-3.0"  "0.0"   "-2.0") # 20 number fraction of accelerated protons log(f_p,RS)

# Add freq data_file to args_list
freq_args=("${freq[@]}")
data_file_args=("${data_file[@]}")
col_args=("${color[@]}")

# Julia script
mkdir -p $version
julia ../src/MagglowNest.jl "$version" "$comment" "$redshift" "$Luminosity_distance" "${freq_args[@]}" "${data_file_args[@]}" "${col_args[@]}" "${P_1[@]}" "${P_2[@]}" "${P_3[@]}" "${P_4[@]}" "${P_5[@]}" "${P_6[@]}" "${P_7[@]}" "${P_8[@]}" "${P_9[@]}" "${P_10[@]}" "${P_11[@]}" "${P_12[@]}" "${P_13[@]}" "${P_14[@]}" "${P_15[@]}" "${P_16[@]}" "${P_17[@]}" "${P_18[@]}" "${P_19[@]}" "${P_20[@]}"