using DifferentialEquations, DiffEqFlux, Plots
using XLSX, DataFrames
using GalacticOptim, Optim
using BlackBoxOptim
using Evolutionary

using Random
Random.seed!(14);

df = DataFrame(XLSX.readtable("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/code_and_dataset_to_share/dataset/dataset_AAV6_07_04_2021.xlsx", "dataset")...)




function ode_system!(du, u, p, t)
    Xv, GLC, GLN, LAC, AMM, AAV1 = u
    μ_Xv, μ_GLC, μ_GLN, μ_LAC, μ_AMM,kdeg, μ_AAV = p
    du[1] = μ_Xv*Xv  #dXv
    du[2] = -μ_GLC*Xv #dGLC
    du[3] = -μ_GLN*Xv #dGLN
    du[4] = μ_LAC*Xv #+ klac1*GLC + klac2*GLN   #dLAC
    du[5] = μ_AMM*Xv + (kdeg*GLN) #- klac1*GLC  #μ_AMM*Xv - kdeg*du[3] #dAMM #(eq10: dAMM)
    du[6] = Xv*μ_AAV #+ kaav1*GLC*Xv + kaav2*GLN*Xv#dtiter
end


phase="before_transfection"
# phase="after_transfection"

tstart=24.0
tend=84
sampling= 12.0
tgrid=tstart:sampling:tend
data_train=[]
u0=[]

if phase=="before_transfection"

   df_train=offline_bioreactor=[
                      0.3595                    32.1954   5.03      0.111019  0.33       0;     %begin of bioreactor dataset
                       0.709                     30.4191   4.7       3.77463   0.91       0;
                        1.3287                    26.4224   3.97      6.66112   1.32       0;
                         1.2691                    24.091    3.54      7.88232   1.46       0;
                          0.8777                    21.7041   2.96      8.27089   1.63       0.0152;
                           1.49295                   20.4829   2.75      8.49292   1.67       0.3015;
                            1.79725                   18.8177   2.49      9.04802   1.65       2.862;
                             1.86                      20.4829   2.75      8.49292   1.67       4.625;
                              1.75425                   14.3769   2.14      10.9353   1.54       5.77];
    #If the xlsx file is not available use the matrix of data above. it is the same!
    df_train = DataFrame(XLSX.readtable("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/code_and_dataset_to_share/dataset/dataset_test_Data_Alinabioreactor 3L_metabolites_cells count_viral titer.xlsx", "dataset")...)
    data_train=df_train[:,2:7]
    data_train=data_train[1:3,:]
    data_train=Array(data_train)
    data_train=convert(Array{Float64,2}, data_train)
    data_train=transpose(data_train)
    # timespanDataTest=df_test[4:9,1]
    time_bioreactor=[0.0 24.0 43.0 57.0 69.0 77.0 90.0 101.0 114.0]
    tstart = time_bioreactor[1];
    tend= time_bioreactor[3]
    sampling= 12.0
    tgrid=time_bioreactor[1:3]
    u0 = [0.3595                 ,   32.1954 ,  5.03    ,  0.111019 , 0.33       ,0.]
else
    df_train = DataFrame(XLSX.readtable("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/code_and_dataset_to_share/dataset/dataset_test_Data_Alinabioreactor 3L_metabolites_cells count_viral titer.xlsx", "dataset")...)
    data_train=df_train[:,2:7]
    data_train=data_train[4:9,:]
    data_train=Array(data_train)
    data_train=convert(Array{Float64,2}, data_train)
    data_train=transpose(data_train)
    # timespanDataTest=df_test[4:9,1]
    time_bioreactor=[0.0 24.0 43.0 57.0 69.0 77.0 90.0 101.0 114.0]
    tstart = time_bioreactor[4]
    tend=time_bioreactor[9]
    sampling= 12.0
    tgrid=time_bioreactor[4:9]
    # tgrid=tstart:sampling:tend
    u0 = [1.27,	24.09103525,	3.54,	  7.882320289,	1.46,	0.0]
end


p =  [0.0067,0.1070,0.0172,0.02543,0.00001,0.0020,0.0757]

loss_over_step=[]

prob = ODEProblem(ode_system!, u0, (tstart,tend), p)
sol=solve(prob,saveat=tstart:sampling:tend)

println("Status: ",sol.retcode)





if phase=="before_transfection"
    bounds = Tuple{Float64, Float64}[
    (0.0295,0.0310),
    (0.1850,0.1950),
    (0.030,0.0403),
    (0.2500,0.3350),
    (0.00001,0.0002),
    (0.0048,0.005),
    (0.000000000000001,0.00000000000001)]
    # μ_AAV1  = 0.0

else
    bounds = Tuple{Float64, Float64}[
    (0.006,0.007),
    (0.01,0.1070),
    (0.0172,0.035),
    (0.01,0.02543),
    (0.0001,0.0002),
    (0.0015,0.0025),
    (.0600, .07)]
end


counter=0
loss_over_step=[]
function cost(p)
  global counter=counter+1
  probdef = ODEProblem(ode_system!, u0, (tstart,tend), p)
  sol = solve(probdef, AutoTsit5(Rosenbrock23()), saveat = tgrid)
  l=sum(abs2, sol.-data_train)
  push!(loss_over_step,l)
  if counter%100000000==0
    plots=Plots.scatter(tgrid,transpose(data_train),color=[:red  :red :red :red :red :red], label = ["[Xv] Offline Bioreactor1"  "[GLC] Offline Bioreactor1" "[GLN] Offline Bioreactor1" "[LAC] Offline Bioreactor1" "[AMM] Offline Bioreactor1" "[AAV] Offline Bioreactor1"], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size = pdim)
    plot!(sol.t,transpose(sol),color=[:red  :red :red :red :red :red], lw=1.5, label = ["[Xv] Bioreactor1 MM"  "[GLC] Bioreactor1 MM" "[GLN] Bioreactor1 MM" "[LAC] Bioreactor1 MM" "[AMM] Bioreactor1 MM" "[AAV] Bioreactor1 MM" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size =pdim)
    display(plots)
    println("Loss:  ",l, "  epoch: ",counter,"    Parameters:",p)
  end
  return l
end

result = bboptimize(cost;SearchRange = bounds, MaxSteps = 2500, Method = :adaptive_de_rand_1_bin_radiuslimited)

p = result.archive_output.best_candidate
plotly();
# gr();

prob = ODEProblem(ode_system!, u0, (tstart,tend), p)
sol=solve(prob, AutoTsit5(Rosenbrock23()), saveat = tgrid)
pdim=(1200,500)

plots=Plots.scatter(tgrid,transpose(data_train),color=[:red  :red :red :red :red :red], label = ["[Xv] Offline Bioreactor1"  "[GLC] Offline Bioreactor1" "[GLN] Offline Bioreactor1" "[LAC] Offline Bioreactor1" "[AMM] Offline Bioreactor1" "[AAV] Offline Bioreactor1"], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size = pdim)
plot!(sol.t,transpose(sol),color=[:red  :red :red :red :red :red], lw=1.5, label = ["[Xv] Bioreactor1 MM"  "[GLC] Bioreactor1 MM" "[GLN] Bioreactor1 MM" "[LAC] Bioreactor1 MM" "[AMM] Bioreactor1 MM" "[AAV] Bioreactor1 MM" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size =pdim)
display(plots)


pll=plot(loss_over_step,xlabel="Iterations",ylabel="loss",label="VVPP UMKM")#, ylims = (0.0,700)
display(pll)
