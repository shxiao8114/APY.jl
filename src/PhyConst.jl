#Constants for Physical Calculation
#Specific Capacity of dry air under same volume
Cva = 719.
#Specific Capacity of water vapor under same volume
Cvv = 1418.
#Specific Capacity of liquid water under same volume
Cvl = 4119.
#Specific Capacity of liquid ice same volume
Cvs = 1861.
E0v = 2.374e6
E0s = 0.3337e6
g = 9.81
PRef = 1.e5
PMax = 1.06e5
P0 = 101325.
Ra = 287.04
Rv = 461.
Rearth = 6372339.
Ttrip = 273.16
Cpa = Cva + Ra
Cpv = Cvv + Rv

#Constants or arrays for tephigram plot
x = [-10., 40.]
xcts = -50. :0.1 : 100.
ThetaEquiv = -10. : 10. : 100.
TempC = -90. : 10. : 40.
DryThetaC = -10. : 10. : 130.
Pres = vcat(1050. : -50. :50.)
SpecHumi = vcat(0.5:0.5:1., 2.:1.:5., [7.,10.,12.,15.], 20.:5.:40., 50.:10.:100.)
