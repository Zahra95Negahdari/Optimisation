clc
clear
close all
% P.iphreeqc.delete();

% Optimization of the Injected Water Composition
% GWA consideration
%% Main Section 
% Using Genethic Algorithm it is aimed to find the best injection salinity
% based on reservoir properties & brine composition

% Calling GA function
global CaseNum;
CaseNum=1;

[x,Fval] = GAOptInjComp();

%% Genetic Algorithm
%  Using a new function to apply integer values
function [x,Fval]=GAOptInjComp()
    nvars = 6;   % Number of decision variables
    populationSize = 150; % Population Size
    generations = 5;   % Number of Generations
    % C_Ca   C_Cl  C_Mg  C_Na  C_SO4   C_HCO3 
    lb = [1 1 1 1 1 1 ]; % Lower Bound
    ub = [200 200 200 200 200 100]; % Upper Bound
%     TimeLimit = 100;
    options = gaoptimset('PopulationSize',populationSize,...
        'Generation',generations);
    IntCon = 1:nvars;
    % Multi-objective genetic algorithm
    [x,Fval,exitflag,output,population,scores] = ga(@OFcn,nvars,[],[],[],[],lb,ub,[],IntCon,options);
end

%% Optimization Procedure
function [OF]= OFcn(a)
    global CaseNum;
    % Changing composition in phreeqc
    xlswrite('Output.xlsx',CaseNum,'GA Result',['A' num2str(CaseNum)]);
    xlswrite('Output.xlsx',a,'GA Result',['B' num2str(CaseNum) ':G' num2str(CaseNum)]);
    SI = Compatibility_FCN(a);
    if SI == 1
        [Time,F1] = phreeqc(a);

        % Running Provided Data File in ECLIPSE
        [OF]=ECL(Time,F1);
    else
        OF = 100;
        xlswrite('Output.xlsx',OF,'GA Result',['H' num2str(CaseNum)]);
    end
% OF = 10;
CaseNum = CaseNum +1;
end
%% Compatibility Function
function SI = Compatibility_FCN(a)
    %%%%%%%% Concentration of Injected Water %%%%%%%%%
    C_Ca = a(1)*50; % in ppm
    C_Cl = a(2)*50; % in ppm
    C_Mg = a(3)*50; % in ppm
    C_Na = a(4)*50; % in ppm
    C_SO4 = a(5)*50; % in ppm
    C_HCO3 = a(6)*10; % in ppm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P.iphreeqc = actxserver('IPhreeqcCOM.Object');
    P.iphreeqc.LoadDatabase('C:\Program Files\USGS\IPhreeqcCOM 3.5.0-14000\database\phreeqc.dat');
    P.iphreeqc.OutputFileOn = true;
    % Defining Injection Solution  
    P.iphreeqc.AccumulateLine('SOLUTION 1 injection brine');
    P.iphreeqc.AccumulateLine(' temp      25');
    P.iphreeqc.AccumulateLine(' pressure  316 # atm');
%     P.iphreeqc.AccumulateLine(' pH        5.7');
    P.iphreeqc.AccumulateLine(' pe        4');
    P.iphreeqc.AccumulateLine(' redox     pe');
    P.iphreeqc.AccumulateLine(' units     ppm');
    P.iphreeqc.AccumulateLine(' density   1');
    P.iphreeqc.AccumulateLine([' Ca  '      sprintf('%f',C_Ca)]);
    P.iphreeqc.AccumulateLine([' Cl  '      sprintf('%f',C_Cl)]);
    P.iphreeqc.AccumulateLine([' Mg  '      sprintf('%f',C_Mg)]);
    P.iphreeqc.AccumulateLine([' Na  '      sprintf('%f',C_Na)]);
    P.iphreeqc.AccumulateLine([' S(6)   '     sprintf('%f',C_SO4)]);
    P.iphreeqc.AccumulateLine([' Alkalinity  '  sprintf('%f',C_HCO3) '  as HCO3']);
    P.iphreeqc.AccumulateLine(' -water    1 # kg');
    P.iphreeqc.AccumulateLine('END');    
    
    % Defining Brine Solution
    P.iphreeqc.AccumulateLine('PHASES');
    P.iphreeqc.AccumulateLine('Calcite');
    P.iphreeqc.AccumulateLine('CaCO3(S) + H+ = Ca+2 + HCO3-');
    P.iphreeqc.AccumulateLine('log_k     1.85');
    
    P.iphreeqc.AccumulateLine('SOLUTION 2  initial water  #Al-Shalabi et al. (2015)');
    P.iphreeqc.AccumulateLine('temp      25');
    P.iphreeqc.AccumulateLine('pressure   316# atm');
    P.iphreeqc.AccumulateLine('pH        7');
    P.iphreeqc.AccumulateLine('pe        4');
    P.iphreeqc.AccumulateLine('redox     pe');
    P.iphreeqc.AccumulateLine('units     ppm');
    P.iphreeqc.AccumulateLine('density   1');
    P.iphreeqc.AccumulateLine('Ca        1904');
    P.iphreeqc.AccumulateLine('Cl        132060');
    P.iphreeqc.AccumulateLine('Mg        2439');
    P.iphreeqc.AccumulateLine('Na        59491');
    P.iphreeqc.AccumulateLine('S(6)      350');
    P.iphreeqc.AccumulateLine('-water    1 # kg');
    P.iphreeqc.AccumulateLine('Alkalinity 354 as HCO3');
    P.iphreeqc.AccumulateLine(' EQUILIBRIUM_PHASES 1');
    P.iphreeqc.AccumulateLine(' Calcite   0     10');
    P.iphreeqc.AccumulateLine(' CO2(g)   -3.5    1');
    P.iphreeqc.AccumulateLine('SAVE solution 2');      
    
    P.iphreeqc.AccumulateLine('END');    
    P.iphreeqc.AccumulateLine(' Mix 1'); 
    P.iphreeqc.AccumulateLine('  1   0.5');  
    P.iphreeqc.AccumulateLine('  2   0.5 ');
    
    P.iphreeqc.AccumulateLine('SELECTED_OUTPUT');
    P.iphreeqc.AccumulateLine(' -saturation_indices    Calcite  Anhydrite Aragonite Dolomite Gypsum Halite');    
    P.iphreeqc.AccumulateLine('END');     
%     P.iphreeqc.GetWarningStringLineCount;
    %Run PHREEQC simulation    %(Type "iphreeqc.Lines" to get simulation input)
    try  P.iphreeqc.RunAccumulated();

        %Fetch Results (Faster than doing GetSelectedOutPutValue each time)
        result = P.iphreeqc.GetSelectedOutputArray();
    %     P.GetSelectedOutputValue(3,4);
        %Clear phreeqc accumulated (input) lines
        P.iphreeqc.ClearAccumulatedLines();

        %Send output back to MATLAB transport simulator
        outputComp = result;
        result(:,1:8) = [];
        result(1,:) = [];
        result = cell2mat(result);
        result = round(result,1);
        SI = find(result>0);
        SI = round(SI,2);
        P.iphreeqc.delete();    
    catch
        SI = 100;
    end
    
    if isempty(SI)
        SI = 1;
    else
        SI = 0;
    end
end

%% Phreeqc  
function [Time, F1] = phreeqc(a)
    global CaseNum;
    
    %%%%%%%% Concentration of Injected Water %%%%%%%%%
    C_Ca = a(1)*50; % in ppm
    C_Cl = a(2)*50; % in ppm
    C_Mg = a(3)*50; % in ppm
    C_Na = a(4)*50; % in ppm
    C_SO4 = a(5)*50; % in ppm
    C_HCO3 = a(6)*10; % in ppm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    P.iphreeqc = actxserver('IPhreeqcCOM.Object');
    P.iphreeqc.LoadDatabase('C:\Program Files\USGS\IPhreeqcCOM 3.5.0-14000\database\phreeqc.dat');
    P.iphreeqc.OutputFileOn = true;
    % defining Surface Master Species & Surface Species
    P.iphreeqc.AccumulateLine('SURFACE_MASTER_SPECIES');
    P.iphreeqc.AccumulateLine('Surfa_p Surfa_pOH');
    P.iphreeqc.AccumulateLine('Surfa_n Surfa_nH');
    P.iphreeqc.AccumulateLine('Surfb_n Surfb_nH #-COO-;-COOH');
    P.iphreeqc.AccumulateLine('Surfb_a Surfb_aH+ #-N site; -NH+  #Korrani et al. (2019)');
    
    P.iphreeqc.AccumulateLine('SURFACE_SPECIES');
    P.iphreeqc.AccumulateLine('Surfa_pOH = Surfa_pOH');
    P.iphreeqc.AccumulateLine('log_k     0');
    P.iphreeqc.AccumulateLine('Surfa_nH = Surfa_nH');
    P.iphreeqc.AccumulateLine('    log_k     0');
    P.iphreeqc.AccumulateLine('Surfa_pOH + H+ = Surfa_pOH2+ # Qiao2014'); 
    P.iphreeqc.AccumulateLine('    log_k     11.8');
    P.iphreeqc.AccumulateLine('    -cd_music 0 1 0 0 0');
    P.iphreeqc.AccumulateLine('SO4-2 + Surfa_pOH2+ = Surfa_pSO4- + H2O # Guo_Haoli 2016'); 
    P.iphreeqc.AccumulateLine('    log_k     2.3');
    P.iphreeqc.AccumulateLine('    -cd_music 0 -1 0 0 0');
    P.iphreeqc.AccumulateLine('CO3-2 + Surfa_pOH2+ = Surfa_pCO3- + H2O # Qiao2015 ');
    P.iphreeqc.AccumulateLine('    log_k     6');
    P.iphreeqc.AccumulateLine('    -cd_music 0 -1 0 0 0');
    P.iphreeqc.AccumulateLine(' Surfa_nH = Surfa_n- + H+ # Qiao2015 ');
    P.iphreeqc.AccumulateLine('    log_k     -5.10');
    P.iphreeqc.AccumulateLine('    -cd_music 0 -1 0 0 0');
    P.iphreeqc.AccumulateLine('Ca+2 + Surfa_n- = Surfa_nCa+ #Qiao2015'); 
    P.iphreeqc.AccumulateLine('    log_k     2.5');
    P.iphreeqc.AccumulateLine('    -cd_music 0 1 0 0 0');
    P.iphreeqc.AccumulateLine('Mg+2 + Surfa_n- = Surfa_nMg+ #Qiao2015 ');
    P.iphreeqc.AccumulateLine('    log_k     2.5');
    P.iphreeqc.AccumulateLine('    -cd_music 0 1 0 0 0');
    % Oil-water interaction
    P.iphreeqc.AccumulateLine('Surfb_aH+ = Surfb_aH+');
    P.iphreeqc.AccumulateLine('log_k 0');
    P.iphreeqc.AccumulateLine('-cd_music 0 0 0');
    P.iphreeqc.AccumulateLine('Surfb_nH = Surfb_nH');
    P.iphreeqc.AccumulateLine('log_k 0');
    P.iphreeqc.AccumulateLine('-cd_music 0 0 0');
    P.iphreeqc.AccumulateLine('Surfb_nH = Surfb_n- + H+  # Qiao2014');
    P.iphreeqc.AccumulateLine('log_k -5');
    P.iphreeqc.AccumulateLine('-cd_music -1 0 0');
    P.iphreeqc.AccumulateLine('Surfb_nH + Ca+2 = Surfb_nCa+ + H+ # Qiao2014');
    P.iphreeqc.AccumulateLine('log_k -3.8');
    P.iphreeqc.AccumulateLine('-cd_music 1 0 0');
    P.iphreeqc.AccumulateLine('Surfb_nH + Mg+2 = Surfb_nMg+ + H+ # Qiao2014');
    P.iphreeqc.AccumulateLine('log_k -4');
    P.iphreeqc.AccumulateLine('-cd_music 1 0 0');
    P.iphreeqc.AccumulateLine('Surfb_aH+ = Surfb_a + H+ #Korrani et al. (2019)');
    P.iphreeqc.AccumulateLine('log_k -6');
    P.iphreeqc.AccumulateLine('-cd_music 0 0 0');
    P.iphreeqc.AccumulateLine('PHASES');
    P.iphreeqc.AccumulateLine('Calcite');
    P.iphreeqc.AccumulateLine('CaCO3(S) + H+ = Ca+2 + HCO3-');
    P.iphreeqc.AccumulateLine('log_k     1.85');
    P.iphreeqc.AccumulateLine('END');
    % Defining Solution #0
    P.iphreeqc.AccumulateLine('SOLUTION 0 injection brine #Al-Shalabi et al. (2015)');
    P.iphreeqc.AccumulateLine(' temp      25');
    P.iphreeqc.AccumulateLine(' pressure  316 # atm');
%     P.iphreeqc.AccumulateLine(' pH        5.7');
    P.iphreeqc.AccumulateLine(' pe        4');
    P.iphreeqc.AccumulateLine(' redox     pe');
    P.iphreeqc.AccumulateLine(' units     ppm');
    P.iphreeqc.AccumulateLine(' density   1');
    P.iphreeqc.AccumulateLine([' Ca  '      sprintf('%f',C_Ca)]);
    P.iphreeqc.AccumulateLine([' Cl  '      sprintf('%f',C_Cl)]);
    P.iphreeqc.AccumulateLine([' Mg  '      sprintf('%f',C_Mg)]);
    P.iphreeqc.AccumulateLine([' Na  '      sprintf('%f',C_Na)]);
    P.iphreeqc.AccumulateLine([' S(6)   '     sprintf('%f',C_SO4)]);
    P.iphreeqc.AccumulateLine([' Alkalinity  '  sprintf('%f',C_HCO3) '  as HCO3']);

%     P.iphreeqc.AccumulateLine(' Ca 325');
%     P.iphreeqc.AccumulateLine(' Cl 9150');
%     P.iphreeqc.AccumulateLine(' Mg 1055');
%     P.iphreeqc.AccumulateLine(' Na 16100');
%     P.iphreeqc.AccumulateLine(' S(6) 4300');
%     P.iphreeqc.AccumulateLine(' Alkalinity  1.2  as HCO3');
  
%     P.iphreeqc.AccumulateLine(' -water    1 # kg');
%     P.iphreeqc.AccumulateLine(' EQUILIBRIUM_PHASES 1');
%     P.iphreeqc.AccumulateLine(' Calcite   0     10');
%     P.iphreeqc.AccumulateLine(' CO2(g)   -3.5    1');
%     P.iphreeqc.AccumulateLine('SAVE solution 0');
    % Defining Solution 1-100
    P.iphreeqc.AccumulateLine('SOLUTION 1-100  initial water  #Al-Shalabi et al. (2015)');
    P.iphreeqc.AccumulateLine('temp      25');
    P.iphreeqc.AccumulateLine('pressure   316# atm');
    P.iphreeqc.AccumulateLine('pH        7');
    P.iphreeqc.AccumulateLine('pe        4');
    P.iphreeqc.AccumulateLine('redox     pe');
    P.iphreeqc.AccumulateLine('units     ppm');
    P.iphreeqc.AccumulateLine('density   1');
    P.iphreeqc.AccumulateLine('Ca        1904');
    P.iphreeqc.AccumulateLine('Cl        132060');
    P.iphreeqc.AccumulateLine('Mg        2439');
    P.iphreeqc.AccumulateLine('Na        59491');
    P.iphreeqc.AccumulateLine('S(6)      350');
    P.iphreeqc.AccumulateLine('-water    1 # kg');
    P.iphreeqc.AccumulateLine('Alkalinity 354 as HCO3');
    P.iphreeqc.AccumulateLine(' -water    1 # kg');
    P.iphreeqc.AccumulateLine(' EQUILIBRIUM_PHASES 1');
    P.iphreeqc.AccumulateLine(' Calcite   0     10');
    P.iphreeqc.AccumulateLine(' CO2(g)   -3.5    1');
    P.iphreeqc.AccumulateLine('SAVE solution 1-100');    
    P.iphreeqc.AccumulateLine('END');
    % Defining Surface 1-100
    P.iphreeqc.AccumulateLine('SURFACE 1-100');
    P.iphreeqc.AccumulateLine('-equilibrate with solution 0');
    P.iphreeqc.AccumulateLine(' -sites_units density');
    P.iphreeqc.AccumulateLine(' -cd_music');
    P.iphreeqc.AccumulateLine(' Surfa_p 3 0.2  10    50');
    P.iphreeqc.AccumulateLine(' Surfa_n 3 0.2  10    50');
    P.iphreeqc.AccumulateLine(' Surfb_n 0.6022     1     10        0.5');
    P.iphreeqc.AccumulateLine(' Surfb_a 0.6022     1     10        0.5');
    P.iphreeqc.AccumulateLine('-capacitance 1e5 0.13');
    P.iphreeqc.AccumulateLine(' -Donnan 3e-10');
    % Reaction Temperature in 1-100
    P.iphreeqc.AccumulateLine('REACTION_TEMPERATURE 1-100');
    P.iphreeqc.AccumulateLine('25');
    % Reaction Pressure 1-100
    P.iphreeqc.AccumulateLine('REACTION_PRESSURE 1-100');
    P.iphreeqc.AccumulateLine('316');
    P.iphreeqc.AccumulateLine('RATES');
    P.iphreeqc.AccumulateLine('Calcite');
    P.iphreeqc.AccumulateLine('-start   #Yutkin et al. (2017)');
    P.iphreeqc.AccumulateLine(' 10 Krxn = 1E-4');
    P.iphreeqc.AccumulateLine(' 20 Ksp  = 1E-12');
    P.iphreeqc.AccumulateLine(' 30 area = 1e-4 #surface area');
    P.iphreeqc.AccumulateLine(' 40 rate = area * Krxn * (((Ksp^(0.5))-ACT("Ca+2") * ACT("CO3-2")) / Ksp^(0.5))');
    P.iphreeqc.AccumulateLine(' 50 moles = rate * Time');
    P.iphreeqc.AccumulateLine(' 60 save moles');
    P.iphreeqc.AccumulateLine(' -end');
    % Defining Kinetics in 1-100
    P.iphreeqc.AccumulateLine('KINETICS 1-100');
    P.iphreeqc.AccumulateLine('Calcite');
    P.iphreeqc.AccumulateLine('-formula   Surfa_pSO4  1  SO4-2   -1   Ca+2  -1');
    P.iphreeqc.AccumulateLine('-m        7e-04');
    P.iphreeqc.AccumulateLine('-m0       7e-04');
    P.iphreeqc.AccumulateLine('-tol      1e-06');
    P.iphreeqc.AccumulateLine('-steps    1');
    P.iphreeqc.AccumulateLine('-step_divide 1');
    P.iphreeqc.AccumulateLine('-runge_kutta   3');
    P.iphreeqc.AccumulateLine('-bad_step_max 100');

    P.iphreeqc.AccumulateLine('INCREMENTAL_REACTIONS True');
    P.iphreeqc.AccumulateLine('SELECTED_OUTPUT');
%     P.iphreeqc.AccumulateLine('-file LOWSAL2.sel');
    P.iphreeqc.AccumulateLine('-step');
    P.iphreeqc.AccumulateLine(['-molalities  SO4-2 CaSO4 CaCO3 Ca+2 Mg+2 CO3-2  Surfa_n-  Surfa_nCa+  Surfa_nH'...
                              '  Surfa_nMg+ Surfa_pCO3-  Surfa_pHCO3  Surfa_pOH  Surfa_pOH2+ Surfa_pSO4-']);
    P.iphreeqc.AccumulateLine('TRANSPORT');
    P.iphreeqc.AccumulateLine('-cells      100');
    P.iphreeqc.AccumulateLine('-shifts     150');
    P.iphreeqc.AccumulateLine('-time_step  688 #second');
    P.iphreeqc.AccumulateLine('-lengths    100*0.001');
    P.iphreeqc.AccumulateLine('-dispersivities   100*0');
    P.iphreeqc.AccumulateLine('-boundary_conditions flux flux');
    P.iphreeqc.AccumulateLine('-thermal_diffusion  1   0');
    P.iphreeqc.AccumulateLine('-punch_cells   1');
    P.iphreeqc.AccumulateLine('-punch_frequency    1');
    P.iphreeqc.AccumulateLine('-print_cells        1');
    P.iphreeqc.AccumulateLine('-print_frequency    1');
    P.iphreeqc.AccumulateLine('#-multi_D  true   1e-11  0.25  .05  1.0');
    % Graph1 - Surfa_pSO4- vs Time
%     P.iphreeqc.AccumulateLine('USER_GRAPH 1 SCM1');
%     P.iphreeqc.AccumulateLine('-headings               Surfa_pSO4- ');
%     P.iphreeqc.AccumulateLine('-axis_titles            "TIME, IN HOURS" "(CONCENTRATION, IN PPB OR MOLAL)" ');
%     P.iphreeqc.AccumulateLine('-chart_title            "SURFACE CONCENTRATIONS"');
%     P.iphreeqc.AccumulateLine('-axis_scale x_axis      0 auto auto auto');
%     P.iphreeqc.AccumulateLine('-axis_scale y_axis      0 auto auto auto');
%     P.iphreeqc.AccumulateLine('-initial_solutions      true');
%     P.iphreeqc.AccumulateLine('-connect_simulations    true');
%     P.iphreeqc.AccumulateLine('-plot_concentration_vs  t');
%     P.iphreeqc.AccumulateLine('  -start');
%     P.iphreeqc.AccumulateLine('10 PUNCH');
%     P.iphreeqc.AccumulateLine('20 x = TOTAL_TIME/3600');
% %     P.iphreeqc.AccumulateLine('#30 PLOT_XY x, MOL("Surfa_p")'); 
%     P.iphreeqc.AccumulateLine('40 PLOT_XY x, MOL("Surfa_pSO4-")');
%     P.iphreeqc.AccumulateLine('#50 PLOT_XY x, MOL("Surfa_pHCO3-")');
%     P.iphreeqc.AccumulateLine('-end');
    % Graph2 F1 vs Time
%     P.iphreeqc.AccumulateLine('USER_GRAPH 2 SCM1');
%     P.iphreeqc.AccumulateLine('-headings               F1');
%     P.iphreeqc.AccumulateLine('-chart_title "F1"');
%     P.iphreeqc.AccumulateLine('-axis_titles "TIME, IN HOURS" "(F1)"');
%     P.iphreeqc.AccumulateLine('-axis_scale x_axis      0 auto auto auto'); 
%     P.iphreeqc.AccumulateLine('-axis_scale y_axis      0 auto auto auto');
%     P.iphreeqc.AccumulateLine('-initial_solutions      true');
%     P.iphreeqc.AccumulateLine('-connect_simulations    true');
%     P.iphreeqc.AccumulateLine('-plot_concentration_vs  t');
%     P.iphreeqc.AccumulateLine('-start');
%     P.iphreeqc.AccumulateLine('10 PUNCH');
%     P.iphreeqc.AccumulateLine('20 x = TOTAL_TIME/3600');
%     P.iphreeqc.AccumulateLine('30 PLOT_XY x,MOL("Surfa_pSO4-") / (TOT("Surfa_p") + TOT("Surfa_n"))');
%     P.iphreeqc.AccumulateLine('-end');
    P.iphreeqc.AccumulateLine('END');

%     P.iphreeqc.GetWarningStringLineCount;
    %Run PHREEQC simulation    %(Type "iphreeqc.Lines" to get simulation input)
    try  P.iphreeqc.RunAccumulated();

        %Fetch Results (Faster than doing GetSelectedOutPutValue each time)
        result = P.iphreeqc.GetSelectedOutputArray();
    %     P.GetSelectedOutputValue(3,4);
        %Clear phreeqc accumulated (input) lines
        P.iphreeqc.ClearAccumulatedLines();

        %Send output back to MATLAB transport simulator
        outputPhreeqc = result;  
        P.iphreeqc.delete();

        for loop_i=4:153
            SurfPN = 0;
            for loop_j=15:23
                x = outputPhreeqc{loop_i,loop_j};
                SurfPN = SurfPN + x;
            end
             F1(loop_i) = outputPhreeqc{loop_i,23}/SurfPN;
             surfpso4(loop_i) = outputPhreeqc{loop_i,23};
             Time(loop_i) = outputPhreeqc{loop_i,5}/3600; % in Day
        end

    %     figure
        subplot(2,2,1);    
        plot(Time,F1);  title(['Case No. = ',num2str(CaseNum)]); xlabel('Time(Hour)'); ylabel('F1(-)'); 
        ylim([0 1]);
        hold on
        subplot(2,2,2);    
        plot(Time,surfpso4); xlabel('Time(Hour)'); ylabel('Surfa_pSO4-(Molal)');
        hold on
    
    catch
        Time = [];
        F1 = [];
    end

   
end

%% Running ECLIPSE
function [OF]=ECL(Time,F1)
    global CaseNum;
    if isempty(Time) == 0 && isempty(F1) == 0
        %% Providing ECL Data-File
        OptEclSim(Time,F1);

        %% Providing Output File
        output = importdata('OptEclSim.RSM');
        output.data(1:2,:)=[];
        output.data(:,1)=[];
        ECL_Output = output.data;
        OF = NPV(ECL_Output);
    else
        OF = 10;
    end
end
 
%% NPV Function
function [OF] = NPV(ECL_Output)
    global CaseNum
    %% ECL_Output
    % Col#1      Col#2      Col#3      Col#4   Col#5
    % Time(Day)  Oil_ProdR  Wat_ProdR  RF      PVinj
    
    %% Ecconomical Parameter Values
    OilPrice = 60; % in $
    LSPrice = 1.5; % in $
    CAPEX = 0.0; % in $
    InterestRate = 0.02;
    %% NPV Calculation
    for loop = 2:length(ECL_Output(:,1))
        delt_Time(loop) = ECL_Output(loop,1)-ECL_Output(loop-1,1);
        OilProd(loop) = ECL_Output(loop,3)-ECL_Output(loop-1,3);
        WatInj(loop) = ECL_Output(loop,4)-ECL_Output(loop-1,4);
    end
        % Dynamic F1
%     output = importdata('OptEclSim.RSM');
%     output = output.data;
%     output(1,:)=[];
%     output(:,1)=[];
%     output(:,9) = (output(1,7)-output(:,7))./output(1,7)*100;
    delt_Time = [0 delt_Time];
    delt_Year = delt_Time/365;
    CashFlow = OilProd.*OilPrice-WatInj.*LSPrice;
    CashFlow = [CAPEX; CashFlow'];
    NPVCalc = CashFlow./(1+InterestRate).^delt_Year';
    NPVtot=NPVCalc(1);
    for loop=2:length(NPVCalc)
        NPVtot(loop) = NPVtot(loop-1) + NPVCalc(loop-1);
    end
    subplot(2,2,3);
    plot(ECL_Output(:,1),ECL_Output(:,5),'-O'); xlabel('Time(Day)'); ylabel('Recovery Factor'); ylim([0 1]);
    hold on    
    subplot(2,2,4);
    Time_new = [0;ECL_Output(:,1)/24];
    OF = 0;
    for loop = 2:length(ECL_Output(:,5))
        OF = OF + (ECL_Output(loop,5)+ECL_Output(loop-1,5))*(ECL_Output(loop,1)-ECL_Output(loop-1,1))/2;
    end
    plot(CaseNum,OF,'-O'); xlabel('Case No.'); ylabel('Average RF'); ylim([0 30]);
%     Time_new = [0;ECL_Output(:,1)/24];
%     plot(Time_new,NPVtot'); xlabel('Time(Day)'); ylabel('NPV($)');
    saveas(gcf,[pwd '\Curves' '\Case_' num2str(CaseNum) '.jpg']);
    close('1');
    %% Objective Function Value
    xlswrite('Output.xlsx',OF,'GA Result',['H' num2str(CaseNum)]);
    OF = 1/OF;
    
end
