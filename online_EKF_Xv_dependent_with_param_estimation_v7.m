

% Equations come from these papers below:
% Paper 1 - The Kalman Filter for the Supervision of Cultivation Processes
% Paper 2 - Development of a Soft-Sensor Based on Multi-Wavelength Fluorescence Spectroscopy and a Dynamic Metabolic Model for Monitoring Mammalian Cell Cultures


%Initialization
clear; close all; clc;
sympref('AbbreviateOutput', false);
%Variable and parameter definition
%Symbols for symbolic math calculations
syms Xv Glc Gln Lac Amm AAV t real
syms uxv uglc ugln ulac uamm kdeg uaav

AAV_data=[
 0.6649                    48.4025   6.18      0.44405   0.12       0.0;
 0.7311                    44.5725   5.63      4.66252   0.59       0.0;
 1.0475                    42.7962   5.12      9.65808   0.95       0.0;
 1.0453                    38.9662   4.39      13.5435   1.23       0.00978;
 1.2085                    35.0807   3.66      15.6528   1.39       0.424;
 1.2261                    29.4189   3.2       16.9849   1.46       1.038;
 1.3959                    33.9151   3.43      22.0915   1.7        1.611;
 1.5437                    29.9185   3.09      23.5346   1.74       1.865;
 2.1302                    25.2004   2.53      20.0933   1.7        2.2415;
 0.6032                    48.236    6.13      0.222025  0.12       0.0;
 0.8325                    45.6271   5.58      4.99556   0.61       0.0;
 0.9847                    42.5187   5.04      9.65808   0.99       0.0;
 0.977                     38.7997   4.33      13.5435   1.29       0.006815;
 1.3077                    37.8006   3.87      16.8739   1.51       0.1355;
 1.3474                    40.2984   4.03      23.2016   1.84       0.543;
 1.5701                    33.8595   3.25      21.2034   1.73       0.952;
 1.6186                    31.5837   3.01      22.9796   1.85       1.0695;
 2.321                     30.7511   2.77      22.9796   1.94       2.0915;
 0.5635                    48.1805   6.28      0.222025  0.12       0.0;
 0.8215                    44.9055   5.67      5.21758   0.6        0.0;
 1.0486                    42.1302   5.21      10.3242   0.97       0.0;
 0.9196                    38.9662   4.4       14.0986   1.25       0.005015;
 1.2261                    25.6444   2.73      11.2123   1.17       0.2475;
 1.3628                    30.5846   3.06      17.0959   1.44       0.7625;
 1.4224                    31.3062   3.14      19.8712   1.62       1.165;
 1.7223                    30.5846   3.18      23.2016   1.84       1.3405;
 2.0376                    28.6418   2.8       22.4245   1.87       1.9955;
 0.7344                    35.5248   2.23      0.222025  0.1        0.0;
 0.6043                    33.0824   1.97      4.66252   0.35       0.0;
 0.9725                    30.085    1.7       9.21403   0.52       0.0;
 0.9792                    28.0313   1.4       12.2114   0.64       0.013;
 1.2504                    19.3166   0.77      10.2131   0.59       0.219;
 1.5348                    24.8118   1.08      17.651    0.73       0.6161;
 1.773                     22.7581   1.03      18.9831   0.77       1.0395;
 2.3618                    20.4267   1.02      20.7593   0.81       1.3365;
 2.6837                    20.1492   0.95      21.6474   0.87       0.411;
 0.6032                    36.0798   9.77      0.44405   0.13       0.0;
 0.6484                    33.1935   9.01      4.4405    0.77       0.0;
 0.9273                    30.8066   8.54      9.10302   1.27       0.0;
 0.8623                    31.6392   8.02      14.6536   1.87       0.01375;
 1.0431                    24.0347   6.45      14.8757   1.96       0.847;
 1.4268                    19.1501   5.21      14.9867   1.86       1.1045;
 1.407                     19.6496   5.52      18.6501   2.15       1.145;
 1.5834                    19.2056   5.74      22.7575   2.53       1.2295;
 1.9748                    18.3175   5.65      23.4236   2.69       1.622;
 0.6418                    61.0582   2.23      0.222025  0.1        0.0;
 0.7763                    57.6167   1.93      4.77353   0.36       0.0;
 1.0155                    56.3956   1.64      9.76909   0.55       0.0;
 0.9174                    47.6809   1.03      12.2114   0.68       0.006755;
 1.2592                    40.798    0.86      12.8774   0.68       0.334;
 1.4621                    53.5092   1.1       21.7584   0.87       0.934;
 1.601                     47.5699   1.03      21.3144   0.84       1.236;
 1.4731                    43.9064   0.89      22.4245   0.89       1.4995;
 2.1225                    44.4615   0.88      23.2016   0.98       1.7655;
 0.5745                    59.5595   9.74      0.222025  0.13       0.0;
 0.7333                    56.6176   9.01      5.21758   0.8        0.0;
 0.9505                    53.3982   8.31      10.4352   1.29       0.0;
 1.0012                    48.68     6.93      14.9867   1.75       0.00915;
 1.3981                    54.5638   7.69      23.0906   2.41       0.582;
 1.6914                    36.0798   4.86      16.5409   1.97       0.491;
 2.0817                    43.9619   5.66      22.3135   2.49       0.4835;
 2.203                     41.2975   5.39      23.8677   2.86       1.0325;
 3.052                     42.1857   5.22      23.9787   3.03       0.305;
 0.3595                    32.1954   5.03      0.111019  0.33       0;     %begin of bioreactor dataset
 0.709                     30.4191   4.7       3.77463   0.91       0;
 1.3287                    26.4224   3.97      6.66112   1.32       0;
 1.2691                    24.091    3.54      7.88232   1.46       0;
 0.8777                    21.7041   2.96      8.27089   1.63       0.0152;
 1.49295                   20.4829   2.75      8.49292   1.67       0.3015;
 1.79725                   18.8177   2.49      9.04802   1.65       2.862;
 1.86                      20.4829   2.75      8.49292   1.67       4.625;
 1.75425                   14.3769   2.14      10.9353   1.54       5.77];  %end of bioreactor dataset

 offline_bioreactor=[
 0.3595                    32.1954   5.03      0.111019  0.33       0;     %begin of bioreactor dataset
 0.709                     30.4191   4.7       3.77463   0.91       0;
 1.3287                    26.4224   3.97      6.66112   1.32       0;
 1.2691                    24.091    3.54      7.88232   1.46       0;
 0.8777                    21.7041   2.96      8.27089   1.63       0.0152;
 1.49295                   20.4829   2.75      8.49292   1.67       0.3015;
 1.79725                   18.8177   2.49      9.04802   1.65       2.862;
 1.86                      20.4829   2.75      8.49292   1.67       4.625;
 1.75425                   14.3769   2.14      10.9353   1.54       5.77];

datasetOnline = csvread("online_data_bioreactor_hspb.csv"); %ONLINE dataset 
datasetOnline = datasetOnline/10^6;

phase ="before_transfection";
phase ="after_transfection";

if phase=="before_transfection";
    offline_bioreactor=offline_bioreactor(1:3,:); %OFFLINE dataset 
else
    offline_bioreactor=offline_bioreactor(4:9,:); %OFFLINE dataset 
end



RUNS=[8 9];%[1 2 3 4 5 6 7 8]
for run = RUNS;%[1 2 3 4 5 6 7 8]
    n = run; % select the index of a run of shakes-flask dataset
    i = 1+((n-1)*9);
    j = 9*n;

    %Variables / Parameters
    if run==8; %
        obs_data=AAV_data(i:end,:);      
        if phase=="before_transfection";
            obs_data = obs_data(1:3,:)'; % Before tran'sfection
            initX = [obs_data(:,1); % initial state (Xv,Glc,Gln,Lac,Amm,AVV)
                0.0299; %initial parameters
                0.1895;
                0.0350;
                0.2544;
                0.0001;
                0.0049;
                0];%0.00157
        else
            obs_data = obs_data(4:9,:)'; % After transfection
            initX = [obs_data(:,1); % initial state (Xv,Glc,Gln,Lac,Amm,AVV)
                0.0065; %initial parameters
                0.0973;
                0.0213;
                0.0214;
                0.0001;%0.0001;
                0.0020;%0.0015;
                0.0644
                 ]  ;          
        end  
        
    else %run 9 ----- 3567 -- 3441 #online dataset
        if phase=="before_transfection";
            initX = [datasetOnline(34); 32.19539379; 5.03000021; 0.111018593;  0.330000013; 0; % initial state (Xv,Glc,Gln,Lac,Amm, AVV)
                0.0299; %initial parameters
                0.1895;
                0.0350;
                0.2544;
                0.0001;
                0.0049;
                0];%0.00157
        else
                %1.0441
        initX = [datasetOnline(3252) ; 26.7219 ;   4.0299  ;  7.2925  ;  1.5469   ;      0; % initial state (Xv,Glc,Gln,Lac,Amm, AVV)
%         initX = [datasetOnline(3252)  ; 26.7205  ;  4.0270  ;  7.3347   ; 1.5466   ;      0; % initial state (Xv,Glc,Gln,Lac,Amm, AVV)
                0.0065; %initial parameters
                0.0973;
                0.0213;
                0.0214;
                0.0001;%0.0001;
                0.0020;%0.0015;
                0.0644
                ];             
        end        
    end
    
    if phase=="before_transfection"
        H2 = [1     0     0     0     0     0     1     1     1     1     1     1     0]; 
        initP = diag([0.0 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00,      1.71e-6,1.53e-6,1.81e-6,2.55e-5,2.97e-9,3.37e-9,0.0]); % initial process Error covarince Matrix
        Q     = diag([0.0006,0.0006,0.0006,0.0006,0.0006,0.0006,     1.71e-7,1.53e-5,1.81e-5,2.55e-4,2.97e-9,3.37e-9,0.0]) ;% process Noise covariance %matrix
    else
        H2 = [1     0     0     0     0     0     0     1     1     1     1     1     1]; % measurement matrix 
        initP = diag([0.0,0.0,0.0,0.0,0.0,0.0,                          7.92e-7, 2.56e-5, 1.05e-5, 9.59e-6, 6.71e-10, 8.71e-8,4.30e-6]); % initial process co-varince 
        Q = diag([0.000006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006,     7.92e-8, 1.56e-5, 1.05e-5, 15.59e-6, 0.11e-8, 0.71e-8,15.30e-6]) ;% process noise covariance %matrix
    
 
    end

    
    init = [initX; initP(:)]; % combined initial value vector for the odesolver
    
    % measurement noise covariance matrix
    R=0.006;

    % Model ODE 
    dS=sym([uxv*Xv;-uglc*Xv;-ugln*Xv;ulac*Xv;uamm*Xv+kdeg*Gln;uaav*Xv;0;0;0;0;0;0;0]);

    %Jacobian of Model with respect to state variables
    F = jacobian(dS,[Xv,Glc,Gln,Lac,Amm,AAV,uxv, uglc, ugln, ulac, uamm, kdeg, uaav]); %Equation 12
    P = sym('P',[13,13]);
    dP = F * P + P*F'+Q; %Equation 11 in paper 1 and Equation 14 in paper 2
    
    %Simulation / State prediction and filtering
    %Assemble all differential equations into a vector of 12 elements %(3x state, 9x P)
    OdeSys = matlabFunction([dS(:);dP(:)],'Vars',{t,[Xv; Glc; Gln; Lac; Amm; AAV; uxv; uglc; ugln; ulac; uamm; kdeg; uaav; P(:)]});

    % Simulate the process from one measurement time to the next:
    MC = zeros(0,13); %store filtered states in these variables
    SimState = zeros(0,13);
    P_state =[];
    SimTime = [];
  
    if run==8 %bioreactor data
        if phase=="before_transfection"
            time_bioreactor=[0 24 43 57 69 77 90 101 114]
            timeE=time_bioreactor(1+1:3);
            timeObsData=time_bioreactor(1:3);
            t0 = time_bioreactor(1);
        else
            time_bioreactor=[0 24 43 57 69 77 90 101 114]
            timeE=time_bioreactor(4+1:9);
            timeObsData=time_bioreactor(4:9);
            t0 = time_bioreactor(4);
        end
    else  %run 9 online data
        if phase=="before_transfection"
            timeE = 34-3:1:3231; % 3231 19-3:1:2919-3;Before tran'sfection
            datasetOnline=datasetOnline(34-3:3231); 
        else
            timeE = 3252:1:6153; % 3043 After transfection
            datasetOnline=datasetOnline(3252:6153);            
        end
        timeE = timeE/60;
        t0 = timeE(1);
        timeObsData = timeE;
        timeE = timeE(2:end);
    end
    

    
    for i = 1:numel(timeE)
        tspan = [t0 timeE(i)];
        [T,state] = ode45(OdeSys, tspan, init); % Solve the system of ODE from t-1 to t (simulate / solve model)
        PS = state(end,1:13)';                  % predicted state from the system of ODE
        P_state = [P_state;PS'];                % Save states for plotting
        if run==9
            MS = datasetOnline(i+1);            % measured state from online dataset   
        else
            MS = obs_data(1,i+1);               % measured state from online dataset      
        end
        P = reshape(state(end,14:end),13,13);   % process error covariance matrix
        K = P*H2'/(H2*P*H2'+ R);                   % kalman gain matrix - Equation 13 in paper1
        FS = PS + K * (MS-PS(1));               % filtered state - Equation 5 and Equation 14 in paper1
        Pfilt = P-K*H2*P;                        % filtered process error covariance matrix - Equation 6 and Equation 15(in paper 1) and Equation 18 (in paper2) 
        init = [FS; Pfilt(:)];                  % new initial condition to be used to solve the system of ODE
        t0 = timeE(i);                          % new starting time for next iteration
        
        % Save intermediate states for plotting
       
        MC = [MC;FS'];
        state(end,1:13) = NaN;
        SimState = [SimState; state(:,1:13)];
        SimTime = [SimTime; T];
    end
    

    
    figure(run);
    counter=1;
    state_names=["[Xv]" "[GLC]" "[GLN]" "[LAC]" "[AMM]" "[AAV]" ];
    for j=[1 2 3 4 5 6]
        subplot(3,3,counter);
       
        plot(timeE,MC(:,j), LineWidth = 2); %plot EKF state variables estimations
        hold on
        
        
        if run==8
            if phase=="before_transfection"
                scatter([0 24 43],offline_bioreactor(:,j)','r','filled'); % plot offline bioreactor bt
                writematrix(MC,'EKF_estimation_bioreactor_bt.csv')
            else
                scatter([57 69 77 90 101 114],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                writematrix(MC,'EKF_estimation_bioreactor_at.csv')
            end
        end        
        
        if run==9
            
            if j==1 %AAV state estimations
                if phase=="before_transfection"
                    plot([1:1:length(datasetOnline)]/60,datasetOnline,LineWidth = 2) %plot xv online bioreactor hspb
                    timeDO=[34-3;1736-3;3182-3]/60;
                    scatter(timeDO,[271326.704120636;501759.15489196;1032084.44843292]/10^6,'g','filled'); %plot xv offline bioreactor hspb bt
                else
                    plot([3252:1:6153]/60,datasetOnline,LineWidth = 2)
                    timeDO=[3252-3;4679-3;6144-3]/60;
                    scatter(timeDO,[1032084.44843292;1218415.005;2631881.10122681]/10^6,'g','filled'); %plot xv offline bioreactor hspb at
                end
            end

            if phase=="before_transfection"
                
                writematrix(offline_bioreactor,'offline_bioreactor_bt.csv')
                writematrix(MC,'EKF_estimation_bioreactor_hSPB_bt.csv')
            
                if j==1
                    scatter([0 24 43],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end
                if j==2
                    scatter([0 24 43],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end
                if j==3
                    scatter([0 24 43],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end
                if j==4
                    scatter([0 24 43],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end
                if j==5
                    scatter([0 24 43],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end
     
            else
                
                writematrix(offline_bioreactor,'offline_bioreactor_at.csv')
                writematrix(MC,'EKF_estimation_bioreactor_hSPB_at.csv')

                if j==1
                    scatter([57 69 77 90 101 114],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end
                if j==2
                    scatter([57 69 77 90 101 114],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                    scatter([4679-3;6144-3]/60,[22.00943658;17.707466],'g','filled'); %plot AVV offline bioreactor hspb           
                end
                if j==3
                    scatter([57 69 77 90 101 114],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end
                if j==4
                    scatter([57 69 77 90 101 114],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                    scatter([4679-3;6144-3]/60,[9.10352484;11.54593394],'g','filled'); %plot AVV offline bioreactor hspb           
                end
                if j==5
                    scatter([57 69 77 90 101 114],offline_bioreactor(:,j)','r','filled'); %plot offline bioreactor at 
                end            
                
               if j==6          
                timeDO=[3252-3;4679-3;6144-3]/60;
                scatter(timeDO,[0;3190000000/10^9;7010000000/10^9],'g','filled'); %plot AVV offline bioreactor hspb           
                end             
            end
           
        end           
      
        xlabel('Time (hours)');
        ylabel(state_names(j));

        counter=counter+1;
          
    end
    
    figure(90);
    plot([1:1:length(P_state(:,7))]/60,P_state(:,7:end),LineWidth = 2);
    xlabel('Time (hours)');
    ylabel('Parameters');
    legend('uxv' ,'uglc', 'ugln', 'ulac', 'uamm' ,'kdeg', 'uaav');
    
   

    
end




