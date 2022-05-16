
using CSV
using Turing, DifferentialEquations, StatsPlots, Random,Plots,XLSX,DataFrames
# Set a seed for reproducibility.
using Random
Random.seed!(14);

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
    #if the .xlsx file is not available use the matrix above with the data. it is the same!
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

# time_start=57.0
# time_end=114.0
# tgrid=[57.0 69.0 77.0 90.0 101.0 114.0] #time_start:1:time_end
p =  [0.0067,0.1070,0.0172,0.02543,0.00001,0.0020,0.0757]

prob1 = ODEProblem(ode_system!, u0, (tstart,tend), p)
# predicted1 = solve(prob1,AutoTsit5(Rosenbrock23()),saveat=tgrid)

Turing.setadbackend(:forwarddiff)
@model function fitmb(data, problem)
    μ=0.5
    σ=0.1
    b=0.3
    e=0.8
    σ_p  ~ truncated(Normal(μ, σ),b,e)#InverseGamma(4, 3)

    if phase=="before_transfection"
        μ_Xv    ~ Uniform(0.0295,0.0310)
        μ_GLC   ~ Uniform(0.1850,0.1950)
        μ_GLN   ~ Uniform(0.030,0.0403)
        μ_LAC   ~ Uniform(0.2500,0.3350)
        μ_AMM   ~ Uniform(0.00001,0.0002)
        kdeg    ~ Uniform(0.0048,0.005)
        μ_AAV1  = 0.0

    else
        μ_Xv    ~ Uniform(0.006,0.007)
        μ_GLC   ~ Uniform(0.01,0.1070)
        μ_GLN   ~ Uniform(0.0172,0.035)
        μ_LAC   ~ Uniform(0.01,0.02543)
        μ_AMM   ~ Uniform(0.0001,0.0002)
        kdeg    ~ Uniform(0.0015,0.0025)
        μ_AAV1  ~ Uniform(.0600, .07)
    end

    p = [μ_Xv, μ_GLC, μ_GLN, μ_LAC, μ_AMM,kdeg,μ_AAV1]

    rprob1 = remake(problem, p=p)
    predicted1 = solve(rprob1,AutoTsit5(Rosenbrock23()),saveat=tgrid)
    for  i = 1:length(predicted1)
        data[:,i] ~ MvNormal(predicted1[i], σ_p)#0.5)
    end
end

model = fitmb(data_train, prob1)
chain = sample(model, NUTS(.65),2000)


if phase=="before_transfection"
    μ_Xv = mean(chain["μ_Xv"])
    μ_GLC  = mean(chain["μ_GLC"])
    μ_GLN  = mean(chain["μ_GLN"])
    μ_LAC  = mean(chain["μ_LAC"])
    μ_AMM  = mean(chain["μ_AMM"])
    kdeg  = mean(chain["kdeg"])
    μ_AAV1  = 0.0
    p = [μ_Xv, μ_GLC, μ_GLN, μ_LAC, μ_AMM,kdeg,μ_AAV1]
    # p = [0.0300, μ_GLC, μ_GLN, μ_LAC, μ_AMM,kdeg,μ_AAV1]
else
    μ_Xv = mean(chain["μ_Xv"])
    μ_GLC  = mean(chain["μ_GLC"])
    μ_GLN  = mean(chain["μ_GLN"])
    μ_LAC  = mean(chain["μ_LAC"])
    μ_AMM  = mean(chain["μ_AMM"])
    kdeg  = mean(chain["kdeg"])
    μ_AAV1  = mean(chain["μ_AAV1"])
    p = [μ_Xv, μ_GLC, μ_GLN, μ_LAC, μ_AMM,kdeg,μ_AAV1]
end


prob1 = ODEProblem(ode_system!, u0, (tstart,tend), p)
sol = solve(prob1, AutoTsit5(Rosenbrock23()),saveat=tgrid)

plotly()

display(chain)

#plot chain
display(plot(chain))

#plot the ACF for the chain
display(autocorplot(chain))

if phase=="before_transfection"
    var_μ_Xv = var(chain["μ_Xv"])
    var_μ_GLC  = var(chain["μ_GLC"])
    var_μ_GLN  = var(chain["μ_GLN"])
    var_μ_LAC  = var(chain["μ_LAC"])
    var_μ_AMM  = var(chain["μ_AMM"])
    var_μ_AAV1  = 0.0
    var_kdeg  = var(chain["kdeg"])
    var_p = [var_μ_Xv, var_μ_GLC, var_μ_GLN, var_μ_LAC, var_μ_AMM,var_kdeg,var_μ_AAV1]
else
    var_μ_Xv = var(chain["μ_Xv"])
    var_μ_GLC  = var(chain["μ_GLC"])
    var_μ_GLN  = var(chain["μ_GLN"])
    var_μ_LAC  = var(chain["μ_LAC"])
    var_μ_AMM  = var(chain["μ_AMM"])
    var_μ_AAV1  = var(chain["μ_AAV1"])
    var_kdeg  = var(chain["kdeg"])
    var_p = [var_μ_Xv, var_μ_GLC, var_μ_GLN, var_μ_LAC, var_μ_AMM,var_kdeg,var_μ_AAV1]
end
println("Variance")
display(var_p)


############################################## Analysis ##############################################

# These files come from Matlab codes
# writematrix(offline_bioreactor,'offline_bioreactor_at.csv')
# writematrix(offline_bioreactor,'offline_bioreactor_at.csv')
# writematrix(MC,'EKF_estimation_bioreactor_at.csv')
EKF_estimation_bioreactor_hSPB_at = CSV.read("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/EKF/EKF_Cristovao/EKF_estimation_bioreactor_hSPB_at.csv", DataFrame, header=false)
df_EKF_estimation_bioreactor_at = CSV.read("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/EKF/EKF_Cristovao/EKF_estimation_bioreactor_at.csv", DataFrame, header=false)
df_offline_bioreactor_at = CSV.read("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/EKF/EKF_Cristovao/offline_bioreactor_at.csv", DataFrame, header=false)

EKF_estimation_bioreactor_hSPB_bt = CSV.read("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/EKF/EKF_Cristovao/EKF_estimation_bioreactor_hSPB_bt.csv", DataFrame, header=false)
df_EKF_estimation_bioreactor_bt = CSV.read("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/EKF/EKF_Cristovao/EKF_estimation_bioreactor_bt.csv", DataFrame, header=false)
df_offline_bioreactor_bt = CSV.read("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/EKF/EKF_Cristovao/offline_bioreactor_bt.csv", DataFrame, header=false)

online_data_bioreactor_hspb = CSV.read("/Users/cristovao/Devel/bioprocessesDT/AI4BDT/AAV_production_modeling/EKF/EKF_Cristovao/online_data_bioreactor_hspb.csv", DataFrame, header=false)

# show(df_offline_bioreactor_at)

function RMSE(observed_data,prediction)
    t=observed_data
    y=prediction
    se = (t - y).^2
    mse = mean(se)
    rmse = sqrt(mse)
    return rmse
end

function MSE(observed_data,prediction)
    t=observed_data
    y=prediction
    se = (t - y).^2
    mse = mean(se)

    return mse
end

function SEP(observed_data,prediction)
    sep=RMSE(observed_data,prediction)/maximum(observed_data)
    return sep*100
end


if phase=="before_transfection"
    #estimate the Xv and AAV using the constant parameters to compare with the estimated parameters by EKF
    cu0=[0.2493 , 32.19539379, 5.03000021, 0.111018593,  0.330000013,      0]
    cp=[0.0299,0.1895,0.0350,0.2544,0.0001,0.0049, 0.0]

    ctstart=(34-3)/60.0
    ctend=3231.0/60.0
    csampling= 1.0/60
    ctgrid=ctstart:csampling:ctend

    cprob1 = ODEProblem(ode_system!, cu0, (ctstart,ctend), cp)
    csol = solve(cprob1, AutoTsit5(Rosenbrock23()),saveat=ctgrid)
else
    #estimate the Xv and AAV using the constant parameters to compare with the estimated parameters by EKF
    cu0=[1.0011 , 26.7219  ,  4.0299  ,  7.2925   , 1.5469   ,      0]
    cp=[0.0065,0.0973,0.0213,0.0214,0.0001,0.0020, 0.0644]

    ctstart=3252.0/60
    ctend=6153.0/60
    csampling= 1.0/60
    ctgrid=ctstart:csampling:ctend

    cprob1 = ODEProblem(ode_system!, cu0, (ctstart,ctend), cp)
    csol = solve(cprob1, AutoTsit5(Rosenbrock23()),saveat=ctgrid)
end



gr();
# plotly()

cu0=[1.0011 , 26.7219  ,  4.0299  ,  7.2925   , 1.5469   ,      0]

cp=[0.00998,0.1487,0.0213,0.0718,0.0001,0.0020, 0.1135]
cp=[0.01041,0.1487,0.0213,0.0718,0.0001,0.0020, 0.1135]


ctstart=3252.0/60
ctend=6153.0/60
csampling= 1.0/60
ctgrid=ctstart:csampling:ctend

cprob1 = ODEProblem(ode_system!, cu0, (ctstart,ctend), cp)
ccsol = solve(cprob1, AutoTsit5(Rosenbrock23()),saveat=ctgrid)

# Xv
plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,1],legend=:topleft, color=[:blue], lw=1.5,label = "[Xv] Bioreactor2 EKF", ylabel="[Xv]",xlabel="Time (h)")
plot!([3252:1:6153]/60,Array(online_data_bioreactor_hspb)[3252:6153]/10^6, legend=:topleft,color=[:orange], lw=1.5, label = "[Xv] Online Bioreactor2", ylabel="[Xv]",xlabel="Time (h)")
# plot!(ccsol.t,transpose(ccsol)[:,1],color=[:green ], lw=1.5, label = "[Xv] Bioreactor2 UMKM-P-EKF", ylabel="[Xv]",xlabel="Time (h)")
plot!(csol.t,transpose(csol)[:,1],color=[:red ], lw=1.5, legend=:topleft,label = "[Xv] Bioreactor2 UMKM", ylabel="[Xv]",xlabel="Time (h)")
display(plots)

# AAV
plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,6], legend=:topleft,color=[:blue], lw=1.5,label = "[AAV] Bioreactor2 EKF", ylabel="[AAV]",xlabel="Time (h)")
plot!(csol.t,transpose(csol)[:,6],color=[:red ], legend=:topleft,lw=1.5, label = "[AAV] Bioreactor2 UMKM-P-BI", ylabel="[AAV]",xlabel="Time (h)")
# plot!(ccsol.t,transpose(ccsol)[:,6],color=[:green ], lw=1.5, label = "[AAV] Bioreactor2 UMKM-P-EKF", ylabel="[AAV]",xlabel="Time (h)")
Plots.scatter!([3252-3;4679-3;6144-3]/60,[0;3190000000/10^9;7010000000/10^9], legend=:topleft,color=:orange , label = "[AAV] Offline Bioreactor2" , ylabel="[AAV]" ,xlabel="Time (h)")
display(plots)


# GLC
plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,2],color=[:blue], lw=1.5,label = "[GLC] Bioreactor2 EKF", ylabel="[GLC]",xlabel="Time (h)")
plot!(csol.t,transpose(csol)[:,2],color=[:red ], lw=1.5, label = "[GLC] Bioreactor2 UMKM-P-BI", ylabel="[GLC]",xlabel="Time (h)")
# plot!(ccsol.t,transpose(ccsol)[:,2],color=[:green ], lw=1.5, label = "[GLC] Bioreactor2 UMKM-P-EKF", ylabel="[GLC]",xlabel="Time (h)")
Plots.scatter!([4679-3;6144-3]/60,[22.00943658;17.707466],color=:orange , label = "[GLC] Offline Bioreactor2" , ylabel="[GLC]" ,xlabel="Time (h)")

display(plots)


#LAC
plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,4],color=[:blue],legend=:topleft, lw=1.5,label = "[LAC] Bioreactor2 EKF", ylabel="[LAC]",xlabel="Time (h)")
plot!(csol.t,transpose(csol)[:,4],color=[:red ], legend=:topleft,lw=1.5, label = "[LAC] Bioreactor2 UMKM-P-BI", ylabel="[LAC]",xlabel="Time (h)")
# plot!(ccsol.t,transpose(ccsol)[:,4],color=[:green ], lw=1.5, label = "[LAC] Bioreactor2 UMKM-P-EKF", ylabel="[LAC]",xlabel="Time (h)")
Plots.scatter!([4679-3;6144-3]/60,[9.10352484;11.54593394],color=:orange , legend=:topleft,label = "[LAC] Offline Bioreactor2" , ylabel="[LAC]" ,xlabel="Time (h)")

display(plots)

#





# phase="before_transfection"

pdim=(1200,500)
plotly()
# gr();
if phase=="before_transfection"

    println("EKF estimation for bioreator 1 vs offline bioreator 1 dataset\nRMSEP ---SEP")
    for i=1:6
        # println(df_offline_bioreactor_at[2:end,i],df_EKF_estimation_bioreactor_at[:,i])
        print(  round(RMSE(df_offline_bioreactor_bt[2:end,i],df_EKF_estimation_bioreactor_bt[:,i]); digits = 3))
        println("---",  round(SEP(df_offline_bioreactor_bt[2:end,i],df_EKF_estimation_bioreactor_bt[:,i]); digits = 3))
    end
    println("MM estimation for  bioreator 1 vs offline  bioreator 1 dataset\nRMSEP ---SEP")
    for i=1:6
        # println(df_offline_bioreactor_at[2:end,i],df_EKF_estimation_bioreactor_at[:,i])
        print( round(RMSE(df_offline_bioreactor_bt[2:end,i],transpose(sol)[2:end,i]); digits = 3))
        println("---",  round(SEP(df_offline_bioreactor_bt[2:end,i],transpose(sol)[2:end,i]); digits = 3))
    end


    EKF_estimation_bioreactor_bt=df_EKF_estimation_bioreactor_bt
    EKF_estimation_bioreactor_bt=[0.3595  32.1954  5.03  0.111019  0.33;0.714072608385728 29.8649755992032 4.60059384261552 3.22398363973263 0.900411519555473; 1.31921047098944 26.2383313182102 3.93638542525771 8.00640885173907 1.30227203757861]
    #EKF calibaration
    plots=Plots.scatter(sol.t,transpose(data_train[1:5,:]),color=[:red  :red :red :red :red ], label = ["[Xv] Offline Bioreactor1"  "[GLC] Offline Bioreactor1" "[GLN] Offline Bioreactor1" "[LAC] Offline Bioreactor1" "[AMM] Offline Bioreactor1"], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" ], layout=(2,3),size = pdim)
    plot!(sol.t,transpose(sol),color=[:red  :red :red :red :red ], lw=1.5, label = ["[Xv] Bioreactor1 MM"  "[GLC] Bioreactor1 MM" "[GLN] Bioreactor1 MM" "[LAC] Bioreactor1 MM" "[AMM] Bioreactor1 MM" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" ], layout=(2,3),size =pdim)
    plot!(sol.t[1:end],EKF_estimation_bioreactor_bt,line = 1.5,color=[:blue  :blue :blue :blue :blue], label = ["[Xv] Bioreactor1 EKF"  "[GLC] Bioreactor1 EKF" "[GLN] Bioreactor1 EKF" "[LAC] Bioreactor1 EKF" "[AMM] Bioreactor1 EKF" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" ], layout=(2,3),size = pdim)

    plot!(sol.t[2:end],Array(df_EKF_estimation_bioreactor_bt)[:,1:5],line = 1.5,color=[:blue  :blue :blue :blue :blue], label = ["[Xv] Bioreactor1 EKF"  "[GLC] Bioreactor1 EKF" "[GLN] Bioreactor1 EKF" "[LAC] Bioreactor1 EKF" "[AMM] Bioreactor1 EKF" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" ], layout=(2,3),size = pdim)
    display(plots)

    #EKF test
    # Xv
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,1],color=[:blue], lw=1.5,label = "[Xv] Bioreactor2 EKF", ylabel="[Xv]",xlabel="Time (h)")
        plot!([34-3:1:3231]/60,Array(online_data_bioreactor_hspb)[34-3:3231]/10^6,color=[:orange], lw=1.5, label = "Online Bioreactor2 [Xv]", ylabel="[Xv]",xlabel="Time (h)")
        Plots.scatter!([34-3;1736-3;3182-3]/60,[271326.704120636;501759.15489196;1032084.44843292]/10^6,color=:orange , label = "[XV] Offline Bioreactor2" , ylabel="[Xv]" ,xlabel="Time (h)")
        plot!(csol.t,transpose(csol)[:,1],color=[:red ], lw=1.5, label = "[Xv] Bioreactor2 UMKM", ylabel="[Xv]",xlabel="Time (h)")
        display(plots)

        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,1],color=[:blue], legend=:topleft, lw=1.5,label = "[Xv] Bioreactor2 EKF", ylabel="[Xv]",xlabel="Time (h)")
        plot!([34-3:1:3231]/60,Array(online_data_bioreactor_hspb)[34-3:3231]/10^6,color=[:orange], legend=:topleft, lw=1.5, label = "Online Bioreactor2 [Xv]", ylabel="[Xv]",xlabel="Time (h)")
        display(plots)


# Xv RMSEP and SEP
        EKF_estimation_Xv_bt=[Array(EKF_estimation_bioreactor_hSPB_bt)[1736-3-32,1],Array(EKF_estimation_bioreactor_hSPB_bt)[3182-3-32,1]]
        offline_measurement_Xv_bt=[501759.15489196;1032084.44843292]/10^6

        println("EKF estimation Xv online vs offline dataset\nRMSEP ---SEP")
        print(round(RMSE(offline_measurement_Xv_bt,EKF_estimation_Xv_bt); digits = 3))
        println("---",  round(SEP(offline_measurement_Xv_bt,EKF_estimation_Xv_bt); digits = 3))

        println("EKF estimation Xv UMKM vs offline dataset\nRMSEP ---SEP")
        print(round(RMSE(offline_measurement_Xv_bt,[0.5822039 ,1.196804]); digits = 3))
        println("---",  round(SEP(offline_measurement_Xv_bt,[0.5822039 ,1.196804]); digits = 3))




    # GLC
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,2],color=[:blue], lw=1.5,label = "[GLC] Bioreactor2 EKF", ylabel="[GLC]",xlabel="Time (h)")
        Plots.scatter!(sol.t,transpose(data_train)[:,2],color=:red , label = "[GLC] Offline Bioreactor1" , ylabel="[GLC]" ,xlabel="Time (h)")
        display(plots)
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,2],color=[:blue], lw=1.5,label = "[GLC] Bioreactor2 EKF", ylabel="[GLC]",xlabel="Time (h)")
        display(plots)

    #GLN
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,3],color=[:blue], lw=1.5,label = "[GLN] Bioreactor2 EKF", ylabel="[GLN]",xlabel="Time (h)")
        Plots.scatter!(sol.t,transpose(data_train)[:,3],color=:red , label = "[GLN] Offline Bioreactor1" , ylabel="[GLN]" ,xlabel="Time (h)")
        display(plots)
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,3],color=[:blue], lw=1.5,label = "[GLN] Bioreactor2 EKF", ylabel="[GLN]",xlabel="Time (h)")
        display(plots)
    #LAC
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,4],color=[:blue], lw=1.5,label = "[LAC] Bioreactor2 EKF", ylabel="[LAC]",xlabel="Time (h)")
        Plots.scatter!(sol.t,transpose(data_train)[:,4],color=:red , label = "[LAC] Offline Bioreactor1" , ylabel="[LAC]" ,xlabel="Time (h)")
        display(plots)
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,4],legend=:topleft,color=[:blue], lw=1.5,label = "[LAC] Bioreactor2 EKF", ylabel="[LAC]",xlabel="Time (h)")
        display(plots)
    #AMM
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,5],color=[:blue], lw=1.5,label = "[AMM] Bioreactor2 EKF", ylabel="[AMM]",xlabel="Time (h)")
        Plots.scatter!(sol.t,transpose(data_train)[:,5],color=:red , label = "[AMM] Offline Bioreactor1" , ylabel="[AMM]" ,xlabel="Time (h)")
        display(plots)
        plots=plot([34-3+1:1:3231]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,5],legend=:topleft,color=[:blue], lw=1.5,label = "[AMM] Bioreactor2 EKF", ylabel="[AMM]",xlabel="Time (h)")
        display(plots)
    # #AAV
    #     plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_bt)[:,6],color=[:blue], lw=1.5,label = "EKF [AAV]", ylabel="[AAV]",xlabel="Time (h)")
    #     Plots.scatter!([3252-3;4679-3;6144-3]/60,[0;3190000000/10^9;7010000000/10^9],color=:orange , label = "Offline Bioreactor2 [AAV]" , ylabel="[AAV]" ,xlabel="Time (h)")
    #     display(plots)

    #parameters
        # gr();     cp=[0.0299,0.1895,0.0350,0.2544,0.0001,0.0049, 0.0]



        plots=plot([34-4+1:1:3231]/60,vcat(Array([0.0299,0.1895,0.0350,0.2544,0.0001,0.0049, 0.0])', Array(EKF_estimation_bioreactor_hSPB_bt)[:,7:end]), lw=1.5,label = ["μ_Xv" "μ_GLC" "μ_GLN" "μ_LAC" "μ_AMM" "kdeg" "μ_AAV"], ylabel="Parameters(t)",xlabel="Time (h)")
        display(plots)

        plots=plot(sol.t[1:end],vcat(Array([0.0299,0.1895,0.0350,0.2544,0.0001,0.0049, 0.0])', Array(df_EKF_estimation_bioreactor_bt)[:,7:end]), lw=1.5,label = ["μ_Xv" "μ_GLC" "μ_GLN" "μ_LAC" "μ_AMM" "kdeg" "μ_AAV"], ylabel="Parameters(t)",xlabel="Time (h)")
        display(plots)

else
    println("EKF estimation vs offline dataset\nRMSEP ---SEP")
    for i=1:6
        # println(df_offline_bioreactor_at[2:end,i],df_EKF_estimation_bioreactor_at[:,i])
        print(  round(RMSE(df_offline_bioreactor_at[2:end,i],df_EKF_estimation_bioreactor_at[:,i]); digits = 3))
        println("---",  round(SEP(df_offline_bioreactor_at[2:end,i],df_EKF_estimation_bioreactor_at[:,i]); digits = 3))
    end
    println("MM estimation vs offline dataset\nRMSEP ---SEP")
    for i=1:6
        # println(df_offline_bioreactor_at[2:end,i],df_EKF_estimation_bioreactor_at[:,i])
        print( round(RMSE(df_offline_bioreactor_at[2:end,i],transpose(sol)[2:end,i]); digits = 3))
        println("---",  round(SEP(df_offline_bioreactor_at[2:end,i],transpose(sol)[2:end,i]); digits = 3))
    end

#Bioreactor1 EKF calibaration
    EKF_estimation_bioreactor_at=df_EKF_estimation_bioreactor_at
    EKF_estimation_bioreactor_at=[1.2691000534057617 24.091035248404108 3.54 7.882320288648349 1.46 0.0; 1.29959439242704 22.7160856000123 3.23119708991745 8.17400480941686 1.54244157618275 0.893905061083761; 1.38066150489577 21.7161529959431 3.00413137296532 8.38740330430655 1.59316561992327 1.53778028881235; 1.596863696158 19.7268243789037 2.57618048769347 8.8253602037178 1.66734207756725 2.87235239863282; 1.80074712186719 17.8101360275211 2.17190672187725 9.24951691045143 1.72128323602717 4.1768970541396; 1.86568139434425 15.4891726216375 1.69787989144102 9.77539967921687 1.77450250823754 5.79072831177677]

    plots=Plots.scatter(sol.t,transpose(data_train),color=[:red  :red :red :red :red :red], label = ["[Xv] Offline Bioreactor1"  "[GLC] Offline Bioreactor1" "[GLN] Offline Bioreactor1" "[LAC] Offline Bioreactor1" "[AMM] Offline Bioreactor1" "[AAV] Offline Bioreactor1"], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size = pdim)
    plot!(sol.t,transpose(sol),color=[:red  :red :red :red :red :red], lw=1.5, label = ["[Xv] Bioreactor1 MM"  "[GLC] Bioreactor1 MM" "[GLN] Bioreactor1 MM" "[LAC] Bioreactor1 MM" "[AMM] Bioreactor1 MM" "[AAV] Bioreactor1 MM" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size =pdim)
    plot!(sol.t[1:6],EKF_estimation_bioreactor_at,line = 1.5,color=[:blue  :blue :blue :blue :blue :blue], label = ["[Xv] Bioreactor1 EKF"  "[GLC] Bioreactor1 EKF" "[GLN] Bioreactor1 EKF" "[LAC] Bioreactor1 EKF" "[AMM] Bioreactor1 EKF" "[AAV] Bioreactor1 EKF" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size = pdim)

    # plot!(sol.t[1:6],Array(df_EKF_estimation_bioreactor_at)[:,1:6],line = 1.5,color=[:blue  :blue :blue :blue :blue :blue], label = ["[Xv] Bioreactor1 EKF"  "[GLC] Bioreactor1 EKF" "[GLN] Bioreactor1 EKF" "[LAC] Bioreactor1 EKF" "[AMM] Bioreactor1 EKF" "[AAV] Bioreactor1 EKF" ], ylabel=["[Xv]"  "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]"], layout=(2,3),size = pdim)
    display(plots)

#Bioreactor2 EKF test
# Xv
    plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,1],color=[:blue], lw=1.5,label = "[Xv] Bioreactor2 EKF", ylabel="[Xv]",xlabel="Time (h)")
    plot!([3252:1:6153]/60,Array(online_data_bioreactor_hspb)[3252:6153]/10^6,color=[:orange], lw=1.5, label = "[Xv] Online Bioreactor2", ylabel="[Xv]",xlabel="Time (h)")
    plot!(csol.t,transpose(csol)[:,1],color=[:red ], lw=1.5, label = "[Xv] Bioreactor2 UMKM", ylabel="[Xv]",xlabel="Time (h)")

    # Plots.scatter!([3252-3;4679-3;6144-3]/60,[1032084.44843292;1218415.005;2631881.10122681]/10^6,color=:orange , label = "[XV] Offline Bioreactor2" , ylabel="[Xv]" ,xlabel="Time (h)")
    display(plots)

# GLC
    plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,2],color=[:blue], lw=1.5,label = "[GLC] Bioreactor2 EKF", ylabel="[GLC]",xlabel="Time (h)")
    # Plots.scatter!(sol.t,transpose(data_train)[:,2],color=:red , label = "[GLC] Offline Bioreactor1" , ylabel="[GLC]" ,xlabel="Time (h)")
    Plots.scatter!([4679-3;6144-3]/60,[22.00943658;17.707466],color=:orange , label = "[GLC] Offline Bioreactor2" , ylabel="[GLC]" ,xlabel="Time (h)")

    display(plots)

#GLN
    plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,3],color=[:blue], lw=1.5,label = "[GLN] Bioreactor2 EKF", ylabel="[GLN]",xlabel="Time (h)")
    Plots.scatter!(sol.t,transpose(data_train)[:,3],color=:red , label = "[GLN] Offline Bioreactor1" , ylabel="[GLN]" ,xlabel="Time (h)")
    Plots.scatter!([3252-3;4679-3;6144-3]/60,[0;2.42;4.68],color=:orange , label = "[GLN] Offline Bioreactor2" , ylabel="[GLN]" ,xlabel="Time (h)")

    display(plots)

#LAC
    plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,4],color=[:blue], lw=1.5,label = "[LAC] Bioreactor2 EKF", ylabel="[LAC]",xlabel="Time (h)")
    # Plots.scatter!(sol.t,transpose(data_train)[:,4],color=:red , label = "[LAC] Offline Bioreactor1" , ylabel="[LAC]" ,xlabel="Time (h)")
    Plots.scatter!([4679-3;6144-3]/60,[9.10352484;11.54593394],color=:orange , label = "[LAC] Offline Bioreactor2" , ylabel="[LAC]" ,xlabel="Time (h)")

    display(plots)

#AMM
    plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,5],color=[:blue], lw=1.5,label = "[AMM] Bioreactor2 EKF", ylabel="[AMM]",xlabel="Time (h)")
    Plots.scatter!(sol.t,transpose(data_train)[:,5],color=:red , label = "[AMM] Offline Bioreactor1" , ylabel="[AMM]" ,xlabel="Time (h)")
    Plots.scatter!([3252-3;4679-3;6144-3]/60,[0;0.255;0.255],color=:green , label = "[AMM] Offline Bioreactor2" , ylabel="[AMM]" ,xlabel="Time (h)")

    display(plots)

#AAV
    plots=plot([3252+1:1:6153]/60,Array(EKF_estimation_bioreactor_hSPB_at)[:,6],color=[:blue], lw=1.5,label = "[AAV] Bioreactor2 EKF", ylabel="[AAV]",xlabel="Time (h)")
    plot!(csol.t,transpose(csol)[:,6],color=[:red ], lw=1.5, label = "[AAV] Bioreactor2 UMKM", ylabel="[AAV]",xlabel="Time (h)")
    Plots.scatter!([3252-3;4679-3;6144-3]/60,[0;3190000000/10^9;7010000000/10^9],color=:orange , label = "[AAV] Offline Bioreactor2" , ylabel="[AAV]" ,xlabel="Time (h)")
    display(plots)

#RMSEP and SEP
    UMKM_P_BI_AAV=[1.6544,3.6449]
    UMKM_P_EKF_AAV=[3.0063,6.8479]
    EKF_AAV=[2.6915,7.0700]
    offline_measurement_aav_at=[3.19,7.01]

    UMKM_P_BI_GLC=[24.2208,21.2134]
    UMKM_P_EKF_GLC=[22.7817,17.7487]
    EKF_GLC=[23.0718,17.4199]
    offline_measurement_GLC_at=[22.0094,17.7074]

    UMKM_P_BI_LAC=[7.8844,8.5459]
    UMKM_P_EKF_LAC=[9.2365,11.6667]
    EKF_LAC=[8.8231,11.7060]
    offline_measurement_LAC_at=[9.1035,11.5459]

    println("RMSE FOR AAV - UMKM_P_BI, EKF, UMKM_P_EKF \n")
    print(round(RMSE(offline_measurement_aav_at,UMKM_P_BI_AAV); digits = 3),"  ")
    print(round(RMSE(offline_measurement_aav_at,EKF_AAV); digits = 3),"  ")
    print(round(RMSE(offline_measurement_aav_at,UMKM_P_EKF_AAV); digits = 3),"  ")

    println("RMSE FOR GLC - UMKM_P_BI, EKF, UMKM_P_EKF \n")
    print(round(RMSE(offline_measurement_GLC_at,UMKM_P_BI_GLC); digits = 3),"  ")
    print(round(RMSE(offline_measurement_GLC_at,EKF_GLC); digits = 3),"  ")
    print(round(RMSE(offline_measurement_GLC_at,UMKM_P_EKF_GLC); digits = 3),"  ")

    println("RMSE FOR LAC - UMKM_P_BI, EKF, UMKM_P_EKF \n")
    print(round(RMSE(offline_measurement_LAC_at,UMKM_P_BI_LAC); digits = 3),"  ")
    print(round(RMSE(offline_measurement_LAC_at,EKF_LAC); digits = 3),"  ")
    print(round(RMSE(offline_measurement_LAC_at,UMKM_P_EKF_LAC); digits = 3),"  ")

#parameters
    # gr();
    # vcat(Array([0.0065,0.0973,0.0213,0.0214,0.0001,0.0020, 0.0644])', Array(df_EKF_estimation_bioreactor_at)[:,7:end])

    plots=plot([3252:1:6153]/60,vcat(Array([0.0065,0.0973,0.0213,0.0214,0.0001,0.0020, 0.0644])', Array(EKF_estimation_bioreactor_hSPB_at)[:,7:end]), lw=1.5,label = ["μ_Xv" "μ_GLC" "μ_GLN" "μ_LAC" "μ_AMM" "kdeg" "μ_AAV"], ylabel="Parameters(t)",xlabel="Time (h)")
    display(plots)

    plots=plot(sol.t[1:end],vcat(Array([0.0065,0.0973,0.0213,0.0214,0.0001,0.0020, 0.0644])', Array(df_EKF_estimation_bioreactor_at)[:,7:end]), lw=1.5,label = ["μ_Xv" "μ_GLC" "μ_GLN" "μ_LAC" "μ_AMM" "kdeg" "μ_AAV"], ylabel="Parameters(t)",xlabel="Time (h)")
    display(plots)
end
