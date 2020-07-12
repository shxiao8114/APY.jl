#Calculate Saturated Water Vapor Pressure (Buck Equation, becuase of <0.02% error)
function WVPSat(tempC)
    tempC > 0 ? es = 0.61094 * 1000. .* exp.((18.678 .- (tempC ./ 234.5)) .*
    (tempC ./ (tempC .+ 257.04))) : es = 0.61094 * 1000. .* exp.((23.036 .- (tempC ./ 333.7)) .*
    (tempC ./ (tempC .+ 279.82)))
end

#Calculate Specific Humidity (Mass Water Mixing Ratio in [kg/kg]
function RH2SH(tempC, prs; rh = 100.)
    es = rh * 0.622 * WVPSat(tempC) ./ prs ./100
end

#Virtual temperature, 0.622.*qsh./prs is specific humidity
function TempVir(theta, qsh, prs)
    thetav = theta .* (1.0 .+ 0.622 .* qsh ./ prs)
    return thetav
end

#Calculate Moist Air Constnats
function ActAirConst(qv)
    Cpm = (1 .- qv) .* Cpa .+ qv .* Cpv
    Rm = (1 .- qv) .* Ra .+ qv .* Rv
    return Cpm, Rm
end

#Calculate Potential Temperature
function PotTemp(tempK, prs; qv = 0., inv = 0)
    Cpm, Rm = ActAirConst(qv)
    inv == 0 ? theta = tempK * (PRef / (prs + qv)) ^ (Rm / Cpm) : theta = tempK * (PRef / (prs + qv)) ^ (-Rm / Cpm)
end

#Calculate Latent Energy of H2O(l)->H2O(g)
function Lev(tempC)
    Lv = (2500.8 - 2.36 .* tempC + 1.6e-3 .* tempC .^2 - 6.0e-5 .* tempC .^3) .* 1.0e3
    return Lv
end

#Calculate Equivalent Potential Temperature
function ThetaEqui(tempK, qsh, prs)
    thetaE = (tempK .+ Lev.(tempK .- Ttrip) .* qsh ./ Cpa) .* (PRef ./ prs) .^ (Ra ./ Cpa)
end

function ThetaDry2Wet(thetaK, tempK, qsh)
    thetaE = thetaK .* exp.(Lev.(tempK .- Ttrip) .* qsh ./ (tempK .* Cpa))
end

a = 17.271
b = 237.7 #degC
#Calculate Dew Point Temperature temperature[C] and rh[%]
function gamma(tempC, rh)
    gamma = (a .* tempC ./ (b .+ tempC)) .+ log(rh ./ 100.)
end

using Zygote
function DewPT(tempC, rh; NR = 0)
    if NR == 0
        Td = b .* gamma.(tempC, rh) ./ (a .- gamma.(tempC, rh))
    end
    return Td
end

#Calculate Height of Lifted Condensation Level
using LambertW
function LCLCalc(tempK, prs, rh; h0 = h0)
    es = WVPSat(tempK.-Ttrip)
    qsh = rh .* es / 100.
    qv = 0.622 * qsh ./ prs
    Cpm, Rm = ActAirConst(qv)
    A = Cpm ./ Rm .+ (Cvl .- Cpv) ./ Rv
    B = - (E0v .- (Cvv .- Cvl) .* Ttrip) ./ (Rv .* tempK)
    C = B ./ A
    X = (rh / 100.) .^ (1 ./ A) .* C .* exp.(C)
    Tlcl = C .* tempK ./ (lambertw.(X , -1))
    Hlcl = h0 + Cpm * (tempK - Tlcl) / g
    Plcl = prs * (Tlcl / tempK) ^ (Cpm / Rm)
    return Tlcl, Hlcl, Plcl
end

#Find the vertical index of LCL height
function Findlcl(Hlcl, hmid)
    Indlcl = argmin(abs.(Hlcl .- hmid))
end

#Find the top of Atmospheric Boundary Layer (ABL)
function FindTABL(thetav)
    #find the 
    IndtABL = argmin(abs.(thetav[1:end - 1] .- thetav[end]))
    IndTABL = Int(IndtABL + 0.5 .* (1+(thetav[IndtABL] .- thetav[end]) ./ abs.(thetav[IndtABL] .- thetav[end])))
    return IndTABL
end

#Calculate Convective Available Potential Energy (CAPE)
function CAPE(thetav::Array{Float64,1},Tv::Array{Float64,1},Tlcl::Float64,Hlcl::Float64,hf::Array{Float64,2},prs::Array{Float64,1},qv::Array{Float64,1})
    Indlcl = Findlcl(Hlcl, hf[:,1])
    Cpm, Rm = ActAirConst(qv[Indlcl])
    TvPcl0 = thetav[Indlcl] .* (PRef ./ (prs .+ qv)) .^ (-Rm ./ Cpm)
    Tdiff = TvPcl0 .- Tv
    tfracH = 0.
    h = min(Indlcl,size(thetav)[1]-1)
    while Tdiff[h] > 0.
        tfracH = tfracH + Tdiff[h] / Tv[h] * (hf[h,2] - hf[h,3])
        h -= 1
    end
    CAPEv = tfracH*g
    return CAPEv
end

#Calculate Convective Inhibition(CIN)
function CIN(thetav::Array{Float64,1},Tv::Array{Float64,1},Tlcl::Float64,Hlcl::Float64,hf::Array{Float64,2},prs::Array{Float64,1},qv::Array{Float64,1})
    Indlcl = Findlcl(Hlcl, hf[:,1])
    Cpm, Rm = ActAirConst(qv[Indlcl])
    TvPcl0 = thetav[Indlcl] .* (PRef ./ (prs .+ qv)) .^ (-Rm ./ Cpm)
    Tdiff = TvPcl0 .- Tv
    tfracH = 0.
    h = min(Indlcl,size(thetav)[1]-1)
    while Tdiff[h] < 0.
        tfracH = tfracH + Tdiff[h] / Tv[h] * (hf[h,2] - hf[h,3])
        h -= 1
    end
    CINv = tfracH*g
    return CINv
end