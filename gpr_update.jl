using DataFrames
# import CSV
using Statistics
using Dates
using Plots
using Netatmo
# import PyPlot
using LinearAlgebra
using IterativeSolvers

dtg = Dates.DateTime(2018,05,10,10)
period = Dates.Day(2)
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

sigmat  = 0.5*60*60   # time correlation length scale 
sigmao  = 0.1    # hPa

rbf(utc1,utc2) = exp(-(utc1-utc2)^2/(2*sigmat^2))

# estimator(pressure) = 

p=plot(legend=false)
for station in groupedbyid
   K = [rbf(r1[:time_utc],r2[:time_utc]) for r1 in eachrow(station), r2 in eachrow(station)] 
   KpI = copy(K) 
   KpI[diagind(KpI)] .= 1 + sigmao^2
   meanpres = mean(station[!,:pressure])

   meanpres > 1050 && continue   # the values make the plot unreadable 
   meanpres < 1000 && continue 
   pressureanom = station[!,:pressure] .- meanpres
   q, minreslog = minres(KpI,pressureanom,log=true)
   
   pressurehat = K*q .+ meanpres
   for i in 1:nrow(station)
      station[i,:pressurehat] = pressurehat[i]
   end 
   absdiff = abs.(station[!,:pressurehat]-station[!,:pressure])
   println("$(maximum(absdiff)) $(station[1,:id]) $minreslog")
   ind = absdiff .> 3*sigmao


   if any(ind)
      println("stationerror")
      plot!(p,unix2datetime.(station[!,:time_utc]),station[!,:pressurehat])
      scatter!(p,unix2datetime.(station[!,:time_utc]),station[!,:pressure],markersize=0.3,color=:blue)
      scatter!(p,unix2datetime.(station[ind,:time_utc]),station[ind,:pressure],markersize=1,color=:red)
      
   end 
end 
gui()

scatter(unix2datetime.(groupedbyid[1][!,:time_utc]),groupedbyid[1][!,:pressure],markersize=0.3)
plot!(unix2datetime.(groupedbyid[1][!,:time_utc]),groupedbyid[1][!,:pressurehat])
