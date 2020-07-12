module APY
    include("DatConst.jl")
    include("DatRead.jl")
    include("PhyCalc.jl")
    include("PhyConst.jl")
    include("Tephigram.jl")
end

inname = "src/"

data = APY.datread(inname,"S202005311110174154203/S202005311110174154203.txt")
data1 = filter((row -> row[:Station_Id_C] == 59316), data)
data1 = filter(row -> row[:Day_Data] == 30, data1)
data1 = filter(row -> row[:Hour_Data] == 0, data1)
data1 = filter(row -> row[:GPH] != 999999.0, data1)
data1 = filter(row -> row[:DPT] != 999999.0, data1)
data1 = sort(data1, (:PRS_HWC))

pres = Array(data1[!,:PRS_HWC])*100
temp = Array(data1[!,:TEM]).+APY.Ttrip.+0.1
Td = Array(data1[!,:DPT]).+APY.Ttrip
rh = APY.WVPSat.(Array(data1[!,:DPT]))./APY.WVPSat.(Array(data1[!,:TEM]))*100
hf = Array(data1[!,:GPH])

names(data1)

APY.PltTeph(temp,Td,pres,rh,hf,true,length(pres),"Title")
savefig("example/Demo.pdf")
