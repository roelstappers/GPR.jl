using DataFrames
# import CSV
using Statistics
using Dates
using Plots
using Netatmo
# import PyPlot
using LinearAlgebra
using IterativeSolvers

dtg = Dates.DateTime(2018,05,10,16)
period = Dates.Day(10)
timerange = dtg:Hour(1):dtg+period

latrange  = 59.9:0.01:60  
lonrange  = 10.7:0.01:10.8
latrange  = 59:0.01:60
lonrange  = 7.:0.01:8

df = Netatmo.read(timerange, latrange=latrange, lonrange=lonrange);
println("Finished reading $(size(df,1)) rows")

# Add column to store GP estimated obs
df[!,:pressurehat] .= 0.0

groupedbyid = groupby(df,:id)

sigmat  = 3*60*60   # time correlation length scale 6 hours
sigmao  = 0.2    # hPa

rbf(utc1,utc2) = exp(-(utc1-utc2)^2/(2*sigmat^2))

# estimator(pressure) = 


for station in groupedbyid[1:2]
   K = [rbf(r1[:time_utc],r2[:time_utc]) for r1 in eachrow(station), r2 in eachrow(station)] 
   KpI = copy(K) 
   KpI[diagind(KpI)] .= 1 + sigmao^2
   meanpres = mean(station[!,:pressure])
   pressureanom = station[!,:pressure] .- meanpres
   q, minreslog = minres(KpI,pressureanom,log=true)
   
   pressurehat = K*q .+ meanpres
   for i in 1:nrow(station)
      station[i,:pressurehat] = pressurehat[i]
   end 
   absdiff = abs.(station[!,:pressurehat]-station[!,:pressure])
   println("$(maximum(absdiff)) $(station[1,:id]) $minreslog")
   if any(absdiff .> 3*sigmao)
      println("stationerror")
   end 
end 
scatter(unix2datetime.(groupedbyid[1][!,:time_utc]),groupedbyid[1][!,:pressure],markersize=0.3)
plot!(unix2datetime.(groupedbyid[1][!,:time_utc]),groupedbyid[1][!,:pressurehat])
