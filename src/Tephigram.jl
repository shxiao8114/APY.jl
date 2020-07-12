using Plots
include("PhyCalc.jl")
include("PhyConst.jl")

function PltTemp()
    label = vcat("Isotherm",fill("",size(TempC)[1]))
    for i in 1:size(TempC)[1]
        plot!(x,x.-TempC[i],label = label[i],color = :gold,
        linewidth = 0.3)
    end
end

function PltDryTheta()
    label = vcat("Dry Adiabat",fill("",size(DryThetaC)[1]))
    for i in 1:size(DryThetaC)[1]
        plot!(x,-x.+DryThetaC[i],label = label[i],color = :red,
        linewidth = 0.3)
        annotate!(x[1]+2,-x[1].+DryThetaC[i].-2,text(string("θ = ",Int(DryThetaC[i])),:red,"left",4))
    end
end

function PresPois(Pres)
    Pois = (PRef / (Pres * 100)) ^ (Ra / Cpa)
end

function PresSlope(Pres)
    Anv = PresPois(Pres) - 1
    Apv = PresPois(Pres) + 1
    slope = Anv / Apv
end


function PltPres()
    label = vcat("Isobar",fill("",size(Pres)[1]))
    for i in 1:size(Pres)[1]
        plot!(x, PresSlope(Pres[i]) .* (x .+ Ttrip),label = label[i] ,color = :purple,
        linewidth = 0.3, linestyle = :dash)
        annotate!(x[end]-2, PresSlope(Pres[i]) .* (x[end] .- 2 .+ Ttrip), text(string(Int(Pres[i])," hPa"),:purple,"left",3))
    end
end

# Find points given a specific humidity and changing temperature
function FindPtsSHumi(SpecHumi)
    x0 = zeros(0)
    y0 = zeros(0)
    for j in 1:size(xcts)[1]
        # to convert into [g/kg]
        if  RH2SH(xcts[j], PRef) .* 1000. < SpecHumi
            continue
        else
        #First Calculate Relative Humidity
            rh = 100. .* (SpecHumi ./ 1000.) ./ (0.622 .* WVPSat(xcts[j]) ./ PRef)
            #println(SpecHumi / 1000., " ", WVPSat(x[j]) / PRef)
            #Calculate Lifted Condensation Level
            Tlcl, Hlcl, Plcl = LCLCalc((xcts[j] .+ Ttrip), PRef, rh; h0 = 0.0)
            #println(rh, " ", Tlcl - Ttrip, " ", Plcl)
            slope = PresSlope(Plcl / 100.)
            append!(x0, (((Tlcl - Ttrip) + slope * Ttrip) / (1. - slope)))
            append!(y0, slope * (x0[end] + Ttrip))
        end
    end
    ind = findall(x->((x > -5)&(x < 90)), y0)
    x0 = x0[ind]
    y0 = y0[ind]
    pts = hcat(x0,y0)
    return pts
end

function PltSHumi()
    for i in 1:size(SpecHumi)[1]
        pts = FindPtsSHumi(SpecHumi[i])
        i == 1 ? plot!(pts[:,1],pts[:,2], label = "Specific Humidity", color = :blue, linestyle = :dashdot, linewidth = 0.3) : plot!(pts[:,1],pts[:,2],label = "",color = :blue, linestyle = :dashdot, linewidth = 0.3)
        if i in vcat([3,4,6,7,8],10:1:14)
            annotate!(pts[1,1],pts[1,2], text(string(Int(SpecHumi[i])),:blue,"bottom",3))
        end
    end
end

using Zygote
function dThetadq(tempK, ThetaVK)
    Lv = Lev(tempK - Ttrip)
    prs = PRef * (tempK ./ ThetaVK) ^ 3.5
    es = RH2SH(tempK - Ttrip, prs)
    dtheta_dq = - Lv .* ThetaVK ./ (Cpa .* tempK)
    return dtheta_dq
end

#Find points given a specific humidity and changing temperature
#(borrowing the idea of ...)
function FindPtsEquivTheta(xcts::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, ThetaV::Float64)
    x1 = zeros(0)
    y1 = zeros(0)
    tempC = xcts
    tempK = tempC .+ Ttrip
    ThetaVK = ThetaV .+ Ttrip
    dthetaV = 0.
    prs = PRef * (tempK[1] ./ ThetaVK) ^ 3.5
    es0 = RH2SH(tempK[1] - Ttrip, prs)
    #start integrating numerically
    for j in 2:1:size(tempC)[1]
        prs = PRef * (tempK[j] ./ ThetaVK) ^ 3.5
        es = RH2SH(tempK[j] - Ttrip, prs)
        dthetaV0 = dThetadq(tempK[j] - 0.05, ThetaVK - 0.05) * (es - es0)
        #if (- tempC[j] .+ ThetaV .+ dthetaV) ./ 2. > - 3.
        ThetaVK = ThetaVK + dthetaV0
        ThetaV = ThetaV + dthetaV0
        #println(ThetaV, " ", dthetaV0)
        append!(x1, (+ tempC[j] .+ ThetaV) ./ 2.)
        append!(y1, (- tempC[j] .+ ThetaV) ./ 2.)
        #end
        es0 = es
    end
    ind = findall(x->((x > -5) & (x < 90)), y1)
    x0 = x1[ind]
    y0 = y1[ind]
    pts = hcat(x0,y0)
    return pts
end

function PltEquivTheta()
    label = vcat("Moist Adiabat",fill("",size(ThetaEquiv)[1]))
    for i in 1:size(ThetaEquiv)[1]
        pts = FindPtsEquivTheta(xcts, ThetaEquiv[i])
        plot!(pts[:,1],pts[:,2],label = label[i],color = :green, linestyle = :dashdot, linewidth = 1)
    end
end

##
#This tephigram is designed
using Plots
function PltStruct(title::String)
    p = plot(xaxis = ("[Celsius]",(-10,40), -10:10:40), yaxis = (" ",(-3,90)),
        size = (300,550), framestyle = :grid, legendfontsize = 6, axis = [])
    plot!([40,40,-10,-10,40],[-3,90,90,-3,-3],color = :black, linewidth = 6, label = "")
    PltTemp()
    PltDryTheta()
    PltPres()
    PltSHumi()
    PltEquivTheta()
    return p
end

#Find x,y from (dew point) tempreature and pressure
function FindPtsTP(tempK, prs)
    x2 = zeros(0)
    y2 = zeros(0)
    for i in 1:size(tempK)[1]
        slp = PresSlope(prs[i]./100)
        x = (tempK[i])/ (1 - slp) .- Ttrip
        append!(x2, x)
        append!(y2, slp * (x + Ttrip))
    end
    pts = hcat(x2,y2)
    return pts
end

#Add Temperature and Dew point tempreature, (Dry Adiabat, Conserved Specific Humidity and Moist Adiabat) onto tephigram
function PltTeph(tempK::Array{Float64,1}, tempdK::Array{Float64,1}, prs::Array{Float64,1}, rh::Array{Float64,1}, hf::Array{Float64,1}, add_adiabat::Bool, ind::Int64, title::String)
    p = PltStruct(title::String)
    #tempdK = DewPT.(tempK .- Ttrip, rh) .+ Ttrip
    pts1 = FindPtsTP(tempK, prs)
    #println(pts1)
    pts2 = FindPtsTP(tempdK, prs)
    #println(pts2)
    plot!(pts1[:,1],pts1[:,2], label = "", marker = :circle, markersize = 2, color = :darkred, linewidth = 0.5)
    plot!(pts2[:,1],pts2[:,2], label = "", marker = :circle, markersize = 2, color = :darkblue, linewidth = 0.5)
    if add_adiabat == true
        Tlcl, Hlcl, Plcl = LCLCalc(tempK[ind], prs[ind], rh[ind]; h0=hf[ind])
        dryT = [Tlcl, tempK[ind]]
        riseP = [Plcl, prs[ind]]
        pts3 = FindPtsTP(dryT, riseP)
        wetT = [Tlcl, tempdK[ind]]
        pts4 = FindPtsTP(wetT, riseP)
        qsh = RH2SH(Tlcl - Ttrip, Plcl; rh = 100.)
        println("qsh = ", qsh)
        ThetaK = PotTemp(Tlcl, Plcl)
        ThetaVK = ThetaDry2Wet(ThetaK, Tlcl, qsh)
        tcts = -80 : 0.1 : round(Tlcl - Ttrip, digits = 1)
        pts5 = FindPtsEquivTheta(tcts, ThetaVK - Ttrip)
        while pts5[end,2] < pts3[1,2]
            ThetaVK += 0.05 + 0.1 * (pts3[1,2] - pts5[end,2])
            pts5 = FindPtsEquivTheta(tcts, ThetaVK - Ttrip)
        end
        #println(ThetaVK - Ttrip)
        #println(pts5[end,1]," ", pts5[end,2])
        plot!(pts3[:,1], pts3[:,2], label = "", color = :red, linewidth = 4)
        plot!(pts4[:,1], pts4[:,2], label = "", color = :blue, linewidth = 2, linestyle = :dash)
        plot!(pts5[:,1], pts5[:,2], label = "", color = :green, linewidth = 4)
    end
    return p
end
