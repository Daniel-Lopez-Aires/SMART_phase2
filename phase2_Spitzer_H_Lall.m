%% SEVILLE SPHERICAL TOKAMAK STIMATION
%V2C2 (the winner). PHASE 1, REDUCED SIZE
%
%Daniel López Aires // danlopair@alum.us.es
%Partially merged with Scott's code
%Includes, at the end, eddy plots and stresses stimations

clear 
clc
close all

% %print -depsc2 NOMBREPLOT.eps %this is for saving plots in .eps, without
% saving margins

%#################################################
%######################TFM#########################
%#################################################

%%%PARAMETER TO BE VARIED FOR THE BREAKDOWN!!!!!!!!!!!!!!!!!!!!!
    %1)GAS TYPE

        %i)H    
        %Z_eff=Z_nucleus-Debye shielding. It can be calculated in
        %http://calistry.org/calculate/slaterRuleCalculator , giving:
             Gas_type='H'
             Z_eff=2                                                          %H, accunt Z impurities (Eli)
             C_1=510;                              %[m-1 Torr-1] This is the constant A, I have changed its name
             C_2=1.25e4;                     %[V-1 m-1 Torr-1] This is the constant B, I have changed its name

        %ii)He  
%              Gas_type='He'
%              Z_eff=2 %He
%              C_1=300                                %[m-1 Torr-1] This is the constant A, I have changed its name
%              C_2=3.4e4                        %[V-1 m-1 Torr-1] This is the constant B, I have changed its name
            
        %iii)Ar  
            %Gas_type='Ar'
             %Z_eff=11.85 %Ar
             %C_1=                                %[m-1 Torr-1] This is the constant A, I have changed its name
             %C_2=                                 %[V-1 m-1 Torr-1] This is the constant B, I have changed its name
            
    %2) Greenwald fraction
    
                Gr_fraction=0.15;%0.15;                              %This is to scale the Gr_limit==> <=1
    
    %3) a_eff                                                        %[m] This defines the size of the field null region.
   
             %a_eff=0.05;                                   % little null region WHAT JJ USED, ST25D BASED
             a_eff=0.15;                                        % large null region
            %a_eff=0.3
%%%%%%%%END PARAMETERS TO BE VARIED ON BREAKDWON!!!!!!!!!!!!!

%#################################################
%#####################TFM############
%#################################################


%THINGS TO LOOK AT WHEN MERGING WITH sCOTT'S###############3
% 0) VV,coilset and currents (copy paste it)
% 1) Grid
% 2) Efit geometry

%####################################################






%%%%%%%%%%%%%%%%%  DEFINE DATA OUTPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%
%Taken from the God Scott

%Define figure extension
FigExt = '.png'; 		%'.png','.eps','.pdf'

%Define project and series names
ProjectName = 'S2-000014';			%Define global project name
%SeriesName = 'VaryTauP';		%Define parameter scan series name

%Create global output folders for saved data and figures
ASCIIDir = 'RawData/'; mkdir(ASCIIDir);
FigDir = 'Figures/'; mkdir(FigDir);			%NOT CURRENTLY USED

%Create simulation name based upon relevant run parameters
SimName = 'DefaultSimName';
disp([ 'SimName: ' SimName ]);
disp([ ' ' ]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DEFINE REACTOR GEOMETRY                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Vessel Wall Thickness
VWall_Inboard=0.004;	%Inboard Wall Thickness		[m]
VWall_Outboard=0.008;	%Outboard Wall Thickness	[m]
VWall_Upper=0.015;		%Top Wall Thickness			[m]
VWall_Lower=0.015;		%Bottom Wall Thickness		[m]

%Define Vessel Internal Geometry (Does not include wall thickness)
VesselRMinInner=+0.15+VWall_Inboard;	% R min position [m] 	%Inboard wall 'fixed' by outer edge.
VesselRMaxInner=+0.80;					% R max position [m]
VesselZMinInner=-0.80;					% Z min position [m]
VesselZMaxInner=+0.80;					% Z max position [m]

%Define center points of vessel walls (Inner Geometry + half wall thickness)
ZMinCentre=VesselZMinInner-(VWall_Lower/2);		% Lower Wall 'grows' outwards (-Z direction)
ZMaxCentre=VesselZMaxInner+(VWall_Upper/2);		% Upper Wall 'grows' outwards (+Z direction)
RMinCentre=VesselRMinInner-(VWall_Inboard/2);	% Inboard wall 'grows' inwards (-R direction)
RMaxCentre=VesselRMaxInner+(VWall_Outboard/2);	% Outboard wall 'grows outwards (+R direction)

%%%%%%%%%%%%%%%%%%%%%%%  DEFINE COIL GEOMETRY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Solenoid Geometry and Parameters
nSol = 210;					 		  % Number of Solenoid Windings
RSolInner = 0.115; RSolOuter = 0.145; % Inner and Outer solenoid radii    [m]
RSol = (RSolInner+RSolOuter)/2; % Central radius of solenoid (0.13) [m]
ZMinSol = ZMinCentre;                   % Solenoid Min Z position           [m]
ZMaxSol = ZMaxCentre;                   % Solenoid Max Z position           [m]
ZSol=ZMaxSol+ZMinSol;
SolLength=ZMaxSol-ZMinSol;

%Number of Radial (R) and axial (Z) PF coil windings
nZDiv1=6; nRDiv1=4;
nZDiv2=6; nRDiv2=4;
nZPF1=6; nRPF1=4;
nZPF2=6; nRPF2=4;

%Calculate total number of windings in each coil
nDiv1=nZDiv1*nRDiv1;
nDiv2=nZDiv2*nRDiv2;
nPF1=nZPF1*nRPF1;
nPF2=nZPF2*nRPF2;

turns=[];

iSol=1; iPF1=2; iPF2=3; iDiv1=4; iDiv2=5;
turns(iSol) =nSol; %100 in my tfg
turns(iDiv1) = nDiv1; %8 in my tfg
turns(iDiv2) = nDiv2; %8 in my tfg
turns(iPF1) = nPF1; %24 in my tfg
turns(iPF2) = nPF2; %24 in my tfg

nPF = length(turns);                                %The number of total coils, counting poloidal and inductor


%Define coil total cross-sectional dimensions
width_PF = 0.075;  % Width of the PF coil (m)     (Previously 0.111m) 
height_PF = 0.050; % Height of a the PF coil (m)  (Previously 0.074m)  

%Define central location of coil sets                                          %NOTES  (Outer midplane flange diameter = 176.8mm (180 mm))
R_PF1 = 0.940;  %R position of PF1 (m)	%0.940m     (MINIMUM OF 938mm)
Z_PF1 = 0.200;  %Z Position of PF1 (m)	%0.200m     (MINIMUM OF 233mm)         %Closer together is optimal for +d and -d   (42cm seperation)
R_PF2 = 0.700;  %R Position of PF2 (m)	%0.700m     (MINIMUM OF 938mm)         %Closer to wall is optimal for +d and -d
Z_PF2 = 0.575;  %Z Position of PF2 (m)	%0.575m     (MINIMUM OF 608mm)         %Lower is better for -d, but reduces volume
R_Div1 = 0.250; %R Position of Div1 (m)	%0.250m     (MINIMUM OF 236mm)
Z_Div1 = 0.700; %Z Position of Div1 (m)	%0.700m     (MINIMUM OF 890mm)         %Higher increases aspect ratio (at expense of triang)
R_Div2 = 0.500; %R Position of Div2 (m)	%0.500m     (MINIMUM OF 458mm)
Z_Div2 = 0.700; %Z Position of Div2 (m)	%0.700m     (MINIMUM OF 890mm)         %Lower is better for +d, but reduces volume
  

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    INITIATE VESSEL AND COIL OBJECTS                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global K_B T_prefill
K_B=1.380649*10^(-23);                                          %[JK^-1] Boltzmann constant
T_prefill=373;                                     %[K], temperature of the prefill gas (100º)

 %Creation of the coils%%%%%%%
% the function createVESTPFCircuit (made by Carlos Soria) create two PF
% coils. One in (R, Z) and another in (R, -Z)

%Define Coil density, temperature and resistivity
coil_density = 1;                       % Relative Coil Density      [Arb]
coil_temp = 293.0;                      % Initial Coil Temperature   [K]
resistivity = copper_resistivity_at_temperature(coil_temp);

PF1  = CreateSMARTPFCircuit( 'PF1',R_PF1,Z_PF1,width_PF,height_PF,turns(iPF1),nZPF1,nRPF1,true, coil_temp, resistivity, coil_density);
PF2  = CreateSMARTPFCircuit( 'PF2',R_PF2,Z_PF2,width_PF,height_PF,turns(iPF2),nZPF2,nRPF2,true, coil_temp, resistivity, coil_density);
Div1 = CreateSMARTPFCircuit('Div1', R_Div1, Z_Div1, width_PF,height_PF, turns(iDiv1), nZDiv1,  nRDiv1, true, coil_temp, resistivity, coil_density);
Div2 = CreateSMARTPFCircuit('Div2', R_Div2, Z_Div2, width_PF,height_PF, turns(iDiv2), nZDiv2,  nRDiv2, true, coil_temp, resistivity, coil_density);
    %Create coil set from parameters defined above. (Function made by Carlos Soria)
%Function createVESTPFCircuit creates two PF coils. One in (R, Z) and another in (R, -Z)


%Create array containing number of coil windings - Used to generate coil objects
global coilturns; coilturns=[];
coilturns(iPF1) = nPF1; coilturns(iPF2) = nPF2;
coilturns(iDiv1) = nDiv1; coilturns(iDiv2) = nDiv2;
coilturns(iSol) = nSol; nSolR = 1;

Sol = CreateSMARTSolenoidCircuit('Sol',RSolOuter,RSolInner,ZMaxSol,ZMinSol,coilturns(iSol),nSolR,coil_temp,resistivity,coil_density);

 %%%%%%%%%%%%%%%%%%  INITIATE VACUUM VESSEL FILAMENTS  %%%%%%%%%%%%%%%%%%%

%Define vessel corners, thickness and filament cross-sectional area 
VesselDimensions = [RMinCentre, RMaxCentre, ZMinCentre, ZMaxCentre];       %[m]
WallThickness = [VWall_Upper, VWall_Outboard, VWall_Lower, VWall_Inboard]; %[m]
%Lower filament areas give higher passive current resolution
FilamentArea = 1.50e-4; %(2.5e-4 > A > 1.5e-4 or RZIp M,R matrices fail)   %[m^2]
% ISSUE :: FILAMENT AREA MUST BE CHOSEN VERY CAREFULLY TO ACHIEVE CONVERGENCE - NUMERICAL STABILITY...

%Construct SMART vessel wall filaments ("Static"=fixed fil area, "Diff"=scaled fil area)
[vessel_filament,R_Fil_Array,Z_Fil_Array] = ... 
    CreateRectilinearVessel(VesselDimensions,WallThickness,FilamentArea,"Diff");

%Construct passive vessel components and arrange into a vessel object
%Resistivity and density are set for stainless steel
VesselResistivity = 6.9e-7;  VesselDensity = 7.8e3;	 %[Ohm]; [Kg/m3]; Stainless Steel (GRADE)
global passive; passive = fiesta_passive('STVesselPas',vessel_filament,'g',VesselResistivity,VesselDensity);
global vessel; vessel = fiesta_vessel( 'STVessel',passive);

%%%%%%%%%%%%%%%%%%%%COILSET%%%%%%%%%%%%%%

%Collate global coilset containing Solenoid, PF and Div coil circuits
R_Fil_Array = transpose(R_Fil_Array); Z_Fil_Array = transpose(Z_Fil_Array);     %fiesta_coilset requires a row array
global coilset; coilset = fiesta_coilset('SMARTcoilset',[Sol,PF1,PF2,Div1,Div2],false,R_Fil_Array',Z_Fil_Array');

%Plot of the cross section
    figure;
    set(gca, 'DataAspectRatio', [1,1,1], 'NextPlot', 'add')
    c=plot(coilset);
    set(c, 'EdgeColor', 'k')
    hold on
    c=plot(vessel);
    set(c, 'EdgeColor', 'k')    
    xlabel('R (m)')
    ylabel('Z (m)')
    title('Cross-section')
    %%%OPtionf for tfg
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    %axis([0,1.1,-1.1,1.1]) 
    %print -depsc2 NOMBREPLOT.eps
    Filename = 'CrossSection';
    saveas(gcf, strcat(FigDir,Filename,FigExt));

% @@@@@@@END CREATION OF THE TOKAMAK@@@@@@@@@@@@@@@


%%%%%%%%%%%%%%%%%%%  DEFINE SOL RAMP & COIL CURRENTS  %%%%%%%%%%%%%%%%%%%%%

%Notes:
%Negative coil currents attract the plasma, positive repel the plasma
%Symmetric Solenoid PrePulse and Equil currents aid power supply stability
%Rod current (Irod) sets the toroidal component of the magnetic field.

%Definition of CoilWaveform time intervals:
%time(1)--> All coils and Sol initiate at zero current          Init
%time(2)--> All coils initiate null-field configuration         PrePulse
%time(3)--> All coils maintain null-field configuration         InitRampDown
%time(4)--> Sol ramps down, PF/Div coils init equilibrium       MidRampDown - InitEquil
%time(5)--> Sol completes ramp down, maintain PF/Div coils      EndRampDown - MidEquil
%time(6)--> All coils maintain equilibrium configuration        EndEquil
%time(7)--> All coils and Sol terminate at zero current         Terminate
%%%%%%%
%time(3)-->time(5) lasts timescale TauR (Solenoid Ramp-Down TimeScale)
%time(5)-->time(6) lasts timescale TauP (Pulse/Discharge Timescale)
%%%%%%%

%Solenoid coil currents [kA]		%Phase2     %Phase2NegTri   %Phase2PosTri
I_Sol_Null=+2600;					%+2600;     %+2600;         %+2600;
I_Sol_MidRamp=+0000;				%+0000;     %+0000;         %+0000;
I_Sol_Equil=-1000;                  %-1000;     %-1800;         %-0700;
I_Sol_EndEquil=-1400;           	%-1400;     %-2200;         %-1100;

%PF coil currents (At Equilibrium, time(4,5,6))
I_PF1_Equil=-1100;					%-1100;     %-1100;         %-1100;
I_PF2_Equil=-1100;					%-1100;     %-1100;         %-1100;     (NEG FOR +delta, POS FOR -delta, after efit) 
I_Div1_Equil=+2000;					%+2000;     %-3500;         %+3500;     (HIGH FOR +delta, LOW FOR -delta, before efit)
I_Div2_Equil=+0000;                 %+0000;     %+0000;         %+0000;

%Define number of time-steps (vertices) in the current waveforms
TauB  = 0.016;			% Null Buffer Timescale     [s] Determines null-field buffer
TauR1 = 0.008;			% Breakdown Ramp Timescale  [s] Determines max loop voltage
TauR2 = 0.020;			% PF & Div Ramp Timescale   [s] Determines max PF/Div current ramp
TauR  = TauR1+TauR2;    % Total Ramp Timescale      [s] 
TauP  = 0.100;			% Pulse Timescale      		[s] Determines flat-top timescale
%Time   [Init      PrePulse  InitRampDown  MidRampDown  EndRampDown  MidEquil     Terminate         ];
time =  [-2*TauB   -TauB     0.0           TauR1        TauR         TauR+TauP    TauR+TauP+(2*TauB)];
nTime = length(time);	% Coil Waveform Timesteps	[-]

%Fit any dynamic coil currents, set with 'linear', {pre-ramp, mid-ramp, end-ramp}
I_Sol_MidRamp = FitSolenoidRamp({I_Sol_Null,I_Sol_MidRamp,I_Sol_Equil},time);

%Construct Sol, PF/Div coil current waveforms vertices
%					              %!Null-Field! %!Breakdown!   %!Efit Icoil!
%Time                   [1,       2,                3,            4,             5,             6,             7];
ISol_Waveform =  [0,  I_Sol_Null, I_Sol_Null,   I_Sol_MidRamp, I_Sol_Equil,   I_Sol_EndEquil,0];
IPF1_Waveform =  [0,  NaN,        NaN,          NaN,           I_PF1_Equil,   I_PF1_Equil,   0];
IPF2_Waveform =  [0,  NaN,        NaN,          NaN,           I_PF2_Equil,   I_PF2_Equil,   0];
IDiv1_Waveform = [0,  NaN,        NaN,          NaN,           I_Div1_Equil,  I_Div1_Equil,  0];
    %IDiv1_Waveform = ISol_Waveform;   %IDiv1 in Series with Solenoid
IDiv2_Waveform = [0,  NaN,        NaN,          NaN,           I_Div2_Equil,  I_Div2_Equil,  0];
%%%%%
CoilWaveforms = [ISol_Waveform; IPF1_Waveform; IPF2_Waveform; IDiv1_Waveform; IDiv2_Waveform];

%%%%Loop voltage calc

slope= abs((ISol_Waveform(4)-ISol_Waveform(3))/(time(4)-time(3)));          %[A/t] slope of the ramp down of the Sol

V_loopInner=mu0*turns(iSol)/SolLength*pi*RSolInner^2*slope; %Inner contribution (Rin)
V_loopOuter=mu0*turns(iSol)/SolLength*slope*(-2*pi/(RSolOuter-RSolInner)*(1/3*(RSolOuter^3-RSolInner^3)-RSolInner/2*(RSolOuter^2-RSolInner^2))+pi*(RSolOuter^2-RSolInner^2));                         
            %outer contribution (from Rin to Rout
V_loop=V_loopInner+V_loopOuter


E_Rgeo=V_loop/(2*pi*0.45)            %[V/m] Electric field by Sol only, at RGeo

E= @(R) V_loop./(2*pi*R); %[V/m] Electric field as a function of R, by Sol only!


%%
%%  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@CONFIGURATION OF FIESTA@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%%% Fiesta grid %%%%%%%%%%%%%%%%%%%%
%The grid where FIESTA will run
%Modifying the grid will get different things

R_simulation_limits = [0.01 1.1];                                                    %R limits of the grid(has to contain vessel)
Z_simulation_limits = [-1.3 1.3];         %Z limits of the grid(has to contain vessel)
Grid_size_R = 300;                                                                       %Grid points in R, ST25D based
Grid_size_Z = 251;                                                                       %Grid points in R, ST25D based
Grid = fiesta_grid( R_simulation_limits(1), R_simulation_limits(2),...
    Grid_size_R, Z_simulation_limits(1), Z_simulation_limits(2), Grid_size_Z );

%Extraction of the R and Z value of the points in the grid

rGrid=get(Grid,'r'); %1*200
zGrid=get(Grid,'z'); %1*251
RGrid=get(Grid,'R'); %1*50200, 50200=251*250
ZGrid=get(Grid,'Z'); %1*50200, 50200=251*250

%%% Stimations for the equilibrium%%%%%%%%%%%

%%The excel calculation are introduced here, to avoid using the excel,
%%since I only need few things, that could be easily included in MATLAB.
%%However, there is no change between iterate a few times, or setting
%%general values and carry on

Te=250;                                     %[eV] electron temperature. Calcualted by Eli with my tfg equilibria
Ti=Te*0.1;                                  %[eV] ion temperature
Ip=100e3;                                    %[A] plasma current
RGeo=0.44;                                %[m] This is Rmajor in the excel. From the eq plot!
Rmax=0.7 ;                                 %[m] The max value of R (separatrix), for Z=0. From the eq plot!
a=Rmax-RGeo;                            %[m] minor radius
ZGeo=0;                                      %Major Z, zero due to the symmetry with respect to Z=0
kappa=1.8;                                  %From the eq plot 
A=RGeo/a;%1.909;                       %Aspect ratio
li2=1;                                          %The standard value

Gr_limit=10^20*Ip*10^-6/(pi*a^2*kappa);     % [m^-3] Gr_limit is the plasma density limit (to disrupt)
    
ne=Gr_limit*Gr_fraction;                            %[m^-3] electron density
betaP=3/2*ne*(Te+Ti)/(mu0*Ip/(2*pi*a))^2*2*mu0*1.6*10^-19*kappa;            %pol beta
BT=0.3;                                                      %[T] Toroidal field at Rgeo, plasma geometric centre 
Irod=BT*2*pi*RGeo/mu0;                            %[A] The current needed to achieve BT (toroidal coils)

%Plasma resistance (Spitzer formula)            
log_col=log(12*pi*((8.854*10^-12*1.6*10^-19*Te)^3/(ne*(1.6*10^-19)^6))^(1/2));      %Columb logarithm
plasma_resistance=0.74*Z_eff*1.65*10^-9*log_col/(Te*10^-3)^(3/2); %[Ohm]

 %The value that I was using randomly was 6.9592e-07
                     %Z_eff=4; %VEST VALUE; 4 would be Be==>?¿ (VEST uses
                    %H). May be to take into account impurities form the
                    %wall (A GlobusM article (Plasma formation, 2001) 
                    %about startup says that this affects Zeff)

%%%Plasma model and controls %%%%%%%%%%%%%%%%%%

jprofile = fiesta_jprofile_topeol2( 'Topeol2', betaP, 1, li2, Ip );
control = fiesta_control( 'diagnose',false, 'quiet',false, 'convergence', 1e-5, 'boundary_method',2 );
                            %diagnose was true. If false it do not show
                            %equil plots
%%  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@TARGET EQUILIBRIA@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
config  = fiesta_configuration( 'STV2C2', Grid, coilset);

icoil_equil = fiesta_icoil(coilset);

%Assign equilibrium coil currents to icoil object [kA]
icoil_equil.Sol=CoilWaveforms(iSol,5);         %Solenoid Equilibrium Current at time(5,6)
icoil_equil.PF1=CoilWaveforms(iPF1,5);	%PF1 Equilibrium Current at time(5,6)
icoil_equil.PF2=CoilWaveforms(iPF2,5);	%PF2 Equilibrium Current at time(5,6)
icoil_equil.Div1=CoilWaveforms(iDiv1,5);	%Div1 Equilibrium Current at time(5,6)
    %icoil_equil.Div1=icoil_equil.Sol;         %Div1 in series with Sol
icoil_equil.Div2=CoilWaveforms(iDiv2,5);	%Div2 Equilibrium Current at time(5,6)

%To do the equilibria calc, we can use EFIT algorithm or not:

    %1) NO EFIT..........
    
    %equil = fiesta_equilibrium( 'STV2C2', config, Irod, jprofile, control, [],icoil ); 


    %2)EFIT..........
    %this calculates PF and Divs current given plasma parameters
    %Discovered by Juanjo Toledo Garcia
    
    %[efit_config, signals, weights, index]=fiesta_efit_configuration(config, {'PF1','PF2'}, [0.44, 0, 0.44/1.85 1.8 0.2])
    [efit_config, signals, weights, index]=efit_shape_controller(config, {'PF1','PF2'}, [0.44, 0, 0.44/1.85 1.8 0]);    %DELTA WAS 0.2 BEFORE, NOW TURNED TO 0
        %WHATCH THIS OUT WHEN UPDATING TO SCOTT!!!!!!
    
    % The numbers you give are [Rgeo, Zgeo, a, kappa, delta], Rgeo,Zgeo,a are
    % mandatory.
    %I use the values of the standar shape, to get a similar equil

    equil=fiesta_equilibrium('ST', config, Irod, jprofile, control,efit_config, icoil_equil, signals, weights) %%EFIT!!!
    %It does the case in line 96!! The equil calc is in lin 124

    %close all %to close iterate solver pltos if converge
    
    %Now we have to extract the new currents from the equil, provided that EFIT
    %changed some of them to satisfy the conditions requested:
    icoil_equil=get(equil,'icoil');                         %redefine icoil to save the new current values
    current_post_EFIT=get(icoil_equil,'currents');
    CoilWaveforms(iPF1,5:6) =current_post_EFIT(iPF1);
	CoilWaveforms(iPF2,5:6) =current_post_EFIT(iPF2);
%%%%%%END EFIT................

%%%%%%  Fiesta Plot  %%%%%%%%%%%%%%%%%

 %Plot from the demo in examples of FIESTA folder
%         section_figure=section(equil); %THIS IS A PLOT
%         
         figure;
         plot(equil)        
         hold on
         plot(vessel)
         plot(coilset)
         parametersshow(equil)   %this plots the parameters in the equil
         title('Target equilibria ph2')
         set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
            %Watch out, R0 in the plot is r0_mag, not r0_geom!!
        Filename = 'Target_equilibria';
        saveas(gcf, strcat(FigDir,Filename,FigExt));
        
        clc %to delete all teh warnings that appears

        %Equilibrium parameters

parameters(equil); %its better in this way, because it shows units and coil currents. If you define a variable, it wont do that
param_equil=parameters(equil);                             %Will be used in the null sensors


%% % Make virtual sensors where we want breakdown  %%%%%%%
%Simu2 from 30/6 slides the winner, moved inward same size and shape

                    %InitiateBSensors(EquilParams,length_R,R_centre,Z_centre,length_Z)
sensor_btheta = InitiateBSensors(param_equil,a_eff,0.31); %moved inward

    %r,z of the sensors
r_sensors=get(sensor_btheta,'r'); %size 1*200
z_sensors=get(sensor_btheta,'z'); %size 1*200

    global R_sensor Z_sensor
    [R_sensor,Z_sensor]=meshgrid(r_sensors,z_sensors); %size 200*200
    %(r,z) of the sensors
  
%Plot
        figure;
        plot(equil)        
        hold on
        plot(vessel)
        plot(coilset)
        parametersshow(equil)   %this plots the parameters in the equil
        %title(sprintf('Sensors for simu %d',sen))
        title('Sensors')
        plot(sensor_btheta);
        Filename = 'Sensors';
        %Filename= sprintf('%s_simu_%d',Filename,sen);   
        saveas(gcf, strcat(FigDir,Filename,FigExt));
%%%%%%%%%%%%

%%%%END OF FIESTA EQ@@@@@@@@@@@@@@@@@@@@@@@@@@

%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@RZIp@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%<<<<<<<<<<<<TO DO list
    %*RZIp+EDDYS
    %* Issue with Vp very low and different fromVloop  
    %* Undersntadn Fiesta's RZIp:
        %- Do not match entirely Ati's not/r Lster's JT60S)
        %- Tilde amtrices?
        %- 2 runs??
        %-Inversion to pass from V as input to I?
    %*R,Z eqs have been removed from solver v4, Peter told us
        

%<<<<<<<<<<<<<<<<<<<<<<


% It has the virtual sensors

rzip_config = fiesta_rzip_configuration( 'RZIP', config, vessel, {sensor_btheta} );
[A, B, C, D, curlyM, curlyR, gamma, plasma_parameters, index, label_index, state] = response(rzip_config, equil, 'rp',plasma_resistance);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  Vessel time constant  %%%%%%%%%%%
%This has to deal with the change of eddy currents in the vessel (See the
%article). However, this is not necessary (at least, yet), so I wont use
%it.
 [~,tau_vessel] = eig(curlyR(1:end-3,1:end-3)\curlyM(1:end-3,1:end-3));
 tau_vessel = max(diag(tau_vessel));
 disp([ 'tau_vessel=' num2str(tau_vessel*1e3) 'ms' ]);


%%%%%%  Optimised null  currents%%%%%%%%%%%%
%Extract scaling factors for null-field coil currents - Copied from ST25D Simulation
C_temp = C(end-get(sensor_btheta,'n')+1:end,1:nPF);             %This ic Cn(:,1:nPF) 
                                                  %Cn is the part of the matrix C related to the sensors
                                                  %(see response)
C1 = C_temp(:,iSol);                               %Elements of C_temp(Cn) related to Sol coil

D1_PF1 = C_temp(:,iPF1);                            %Elements of C_temp(Cn) to PF1 coil
D1_PF2 = C_temp(:,iPF2);                            %Elements of C_temp(Cn) to PF2 coil
D1_Div1 = C_temp(:,iDiv1);                          %Elements of C_temp(Cn) to Div1 coil
D1_Div2 = C_temp(:,iDiv2);                          %Elements of C_temp(Cn) to Div2 coil

    %1) Compute all the coil currents (standar)
        D1=C_temp(:,iPF1:end);          %Elements of C_temp(Cn) related to the PF and Div coils
        I_PF_null = -pinv(D1) * (C1*ISol_Waveform(2));	%(PF1,PF2,Div1,Div2)

    %2)Div1 in serie with Sol! :
        %Since we want Div1 in serie with Sol, it should be reomved from the matrix
        %to do the calc:
%         D1=[D1_PF1 D1_PF2 D1_Div2];                     %This is for Div1 in series with Sol, so 
%                                                                     %must not include it in the calculation for the currents                                                                                                                 
%         I_PF_null = -pinv(D1) * (C1*ISol_Waveform(2)+D1_Div1*ISol_Waveform(2));	%Div1 in serie with Sol
%         I_PF_null=[I_PF_null(1); I_PF_null(2); ISol_Waveform(2); I_PF_null(3)]; %To not have problems
         
        %redefining the coil currents

%%%%%%
%Update CoilWaveforms array with null-field values
for i = 1:nPF;
	for j = 1:nTime;
		%Determine if coil 'i' at timestep 'j' requires null-field
		if isnan(CoilWaveforms(i,j));
			CoilWaveforms(i,j) = I_PF_null(i-1);	%i-1 skips Solenoid coil
        end
	end
end
%Note that since Div1 do not have any NaN (is equal to Sol), will not be
%changed!


%%%%%%%%%Construction to do the RZIp calculation%%%%%%%%%%
%The input profile of I_PF currents have to be extended

%Initiate RZip PF Arrays used to calculate plasma and eddy currents
I_PF_input = transpose(CoilWaveforms);  %Coil currents transposed from waveform array
V_PF_input = NaN(nTime,nPF);            %Coil voltages are initiated to zero


%Plot
    figure;
    plot( time*1e3, I_PF_input/(1e3),'*-' );
    xlabel('time (ms)')
    ylabel('I_{{input}} (kA)')
    title('I_{{input}} versus time ph2')
    %title(sprintf('I_{{input}} for simu %d',sen))
    legend('Sol','PF1','PF2','Div1','Div2')
    %%%optinos for tfg
    set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
    set(gca, 'FontSize', 13, 'LineWidth', 0.75);                    %<- Set properties TFG
    Filename = 'Input_currents';
    %Filename= sprintf('%s_simu_%d',Filename,sen);      
    saveas(gcf, strcat(FigDir,Filename,FigExt));
   



%Firstly the input profile is discretized in many intervals
nTime_long = 1000;
time_long = linspace(min(time),max(time),nTime_long);

I_PF_input_long = NaN(nTime_long,nPF);
for iPF=1:nPF
    I_PF_input_long(:,iPF) = interp1(time,I_PF_input(:,iPF),time_long);
end
V_PF_input_long = NaN*I_PF_input_long;                      %it is an I profile, so this is NaN, since it is unkown

%%%%%%%%%
 
%This defines the colours and names of the coils, for the plots
coil_names{iSol} = 'Sol';
coil_names{iDiv1} = 'Div1';
coil_names{iDiv2} = 'Div2';
coil_names{iPF1} = 'PF1';
coil_names{iPF2} = 'PF2';

PF_colors{iSol} = 'Red';
PF_colors{iDiv1} = 'Magenta';
PF_colors{iDiv2} = 'Black';
PF_colors{iPF1} = 'Cyan';
PF_colors{iPF2} = 'Green';


%%RZIP AS FIESTA%%%%%%EXPERIMENTAL%%%%
%{
% %Example as stated on the slides
% t=linspace(0,1,1000)';
% n=get(coilset ,'n');
% v=zeros(length(t),n) ;
% v( t >0.01 & t <0.1 ,1) =100;
% v( t >0.11 & t <0.2 ,2) =100;
% I0=zeros(n,1);
% 
% [T,I_ej] =ode_v(rzip_config , state , A, B, t , v ) ;
% 
% figure;
% subplot(2,4,1)
% plot(T,I_ej(:,index.coilset)*1e-3)
% xlabel('t (s)')
% ylabel('I coilset (kA)')
% title('I coilset (Fiesta example) ')
% subplot(2,4,2)
% plot(T,I_ej(:,index.vessel)*1e-3)
% xlabel('t (s)')
% ylabel('I VV (kA)')
% title('I VV')
% 
% subplot(2,4,3)
% plot(T,I_ej(:,index.coilset)*1e-3)
% xlabel('t (s)')
% ylabel('I coilset (kA)')
% title('I coilset ')
% subplot(2,4,4)
% plot(T,I_ej(:,max(index.vessel)+1)*1e-3)
% xlabel('t (s)')
% ylabel('ZIp (kAm)')
% title('ZIp')
% subplot(2,4,5)
% plot(T,I_ej(:,max(index.vessel)+2)*1e-3)
% xlabel('t (s)')
% ylabel('RIp (kAm)')
% title('RIp')
% subplot(2,4,6)
% plot(T,I_ej(:,max(index.vessel)+2)*1e-3)
% xlabel('t (s)')
% ylabel('Ip (kA)')
% title('Ip')
% 
%     subplot(2,4,7)
%     plot(T,I_ej(:,max(index.vessel)+1)./I_ej(:,max(index.vessel)+3))
%     xlabel('t (s)')
%     ylabel('Z (m)')
%     title('Z (m)')
% 
%     subplot(2,4,8)
%     plot(T,I_ej(:,max(index.vessel)+2)./I_ej(:,max(index.vessel)+3))
%     xlabel('t (s)')
%     ylabel('R (m)')
%     title('R (m)')
% % %%%%%%%%%%%%%%
% % 
% % %Try of doing RZIp as fieta here
% n=get(coilset ,'n');
% V0=zeros(length(I_PF_input_long),1);
% 
%     % I try ode_i, to use I as inputs
%     [T,I] =ode_i(rzip_config , state , curlyM, curlyR, index, time_long' , I_PF_input_long);
% 
%     figure;
%     subplot(4,2,1)
%     plot(T*1e3,I(:,index.coilset)*1e-3)
%     xlabel('t (ms)')
%     ylabel('I coilset (kA)')
%     title('I coilset ')
%     % subplot(3,1,2)
%     % plot(T*1e3,V(:,index.vessel)*1e-3)
%     % xlabel('t (ms)')
%     % ylabel('I VV (kA)')
%     % title('I VV')
%     subplot(4,2,2)
%     plot(T*1e3,sum(I(:,index.vessel),2)*1e-3)
%     xlabel('t (ms)')
%     ylabel('I VV (kA)')
%     title('I VV sum')
%     subplot(4,2,3)
%     plot(T*1e3,I(:,max(index.vessel)+1)*1e-3)
%     xlabel('t (ms)')
%     ylabel('ZIp (kA m)')
%     title('ZIp (kA m)')
% 
%     subplot(4,2,4)
%     plot(T*1e3,I(:,max(index.vessel)+2)*1e-3)
%     xlabel('t (ms)')
%     ylabel('RIp (kAm)')
%     title('RIp (kA m)')
% 
%     subplot(4,2,5)
%     plot(T*1e3,I(:,max(index.vessel)+3)*1e-3)
%     xlabel('t (ms)')
%     ylabel('Ip (kA)')
%     title('Ip')
% 
%     subplot(4,2,6)
%     plot(T*1e3,I(:,max(index.vessel)+1)./I(:,max(index.vessel)+3))
%     xlabel('t (ms)')
%     ylabel('Z (m)')
%     title('Z (m)')
% 
%     subplot(4,2,7)
%     plot(T*1e3,I(:,max(index.vessel)+2)./I(:,max(index.vessel)+3))
%     xlabel('t (ms)')
%     ylabel('R (m)')
%     title('R (m)')
%     
    %%%%Discussion
    %I_VV seems right
    %Ip is not, R is negative, so xD. Maybe there is no contrain of R in
    %the sense that the integration do not stops if the plasma collide
    %with the VV. 
    
    %%%%%END RZIp as Slides, experimental%%%%%%%%%%%%%%%%%
%}    

%%%%%%%RZIp run%%%%%%%%%%%%%
%It is run 2 times. The coils are current driven always, but the plasma changes form
%current driven to voltage driven. The control vector u is defined in
%different ways for plasma current and voltage driven.

%In 1st the run, Ip is zero, and Vp is unknown. u (state vector)
 %contains intensisites The results of the 1st run are:
      %Ip_output zeros, but Vp_output are not null.
      %V_PF_output different to zero, I_PF_output NaN
      %time adaptive (time intervals where the RZIp runs) has negative
      %values. Remember the plasma is created at t=0.
      
%Before the 2nd run, creates Vp_long using Vp_output from the first run,
%but making it zero when time>0 (plasma created at time=0). No Vp have non
%zero values. Ip_log is NaN, and the the second run is made. The output is
%the final output. If show plot is true, I_PF_output=I_PF_input. If not,
%I_PF_output=NaN, since it is not calculated. In the second run, the state
%vector contains I_PF (and gradients), and Vp (non zero).

Ip_long = zeros(size(time_long));
Vp_long = NaN(size(time_long));

[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true );

iTime_plasma = time_adaptive>0;                         %This is a logical, 1 where the condiction is verified, and zero where is not
Vp_output(iTime_plasma) = 0;                                %modified Vp_output where the  previous condition is verified, that is, for time>0. That is, sets Vp=0 when t>0.
Vp_long = interp1( time_adaptive,Vp_output, time_long);         %this finds the Vp for the time intervals time_long, knowing that Vp_ouput corresponds to time_adaptive. Vp_long=0 when time>0.
  
Ip_long = NaN*Vp_long;                                                      %This makes Ip_long unknown (NaN)

[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names', coil_names, 'show_plot',true, 'turns',turns, 'currentScale',1e3, 'PF_colors',PF_colors );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Since the I_Passive contains the eddy at each time for each filament,
%and we want the total eddy, ahve to sum over all the filaments. Since
%each filament is a row, have to sum all the rows

I_Passive_VV=sum(I_Passive,2); 

    %MANUAL PLOTS OF THE RESULTS OF RZIp%%%%%%%

        %%I_PF_output and plasma current
            figure;
            subplot(3,1,1)
            %plot(time_adaptive*1e3,I_PF_output/(1e3))
            %hold on
            plot(time_adaptive*1e3,Ip_output/(1e3))
            ylabel('I (kA)')
            title('Dynamic response SMART phase 1')
            set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
            %set(gca,'YLim',[-5 35]);
            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        %%%vOLTAGE
            subplot(3,1,2)
            plot(time_adaptive*1e3,V_PF_output/(1e3))
            hold on
            plot(time_adaptive*1e3,Vp_output/(1e3))
            ylabel('V (kV)')
            legend('Sol','PF2','PF3','Div1','Div2','Plasma')
            set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
            %set(gca,'YLim',[-1.500 1.500]);
            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        %%%I_PASSIVE VERUS TIME
            subplot(3,1,3)
            plot(time_adaptive*1e3,I_Passive_VV/(1e3))
            xlabel(' time (ms)')
            ylabel('I_{passive} (kA)')
            set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
            %set(gca,'YLim',[-800 800]);
            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
            
            Filename = '_RZIpOutputs';
            %Filename= sprintf('%s_simu_%d',Filename,sen);  
            saveas(gcf, strcat(FigDir,Filename,FigExt));
          %print -depsc2 NOMBREPLOT.eps

          %v{
            figure;           
            plot(time_adaptive*1e3,Ip_output/(1e3))
            ylabel('I_p (kA)')
            xlabel(' time (ms)')
            title('I_p SMART phase 2')
            set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
            %set(gca,'YLim',[-5 35]);
            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
            Filename = 'Ip';
            %Filename= sprintf('%s_Trampdown_%dms',Filename,2*T_ramp_Sol);    
            saveas(gcf, strcat(FigDir,Filename,FigExt));
            
            figure;
            plot(time_adaptive*1e3,I_Passive_VV/(1e3))
            xlabel(' time (ms)')
            ylabel('I_{VV} (kA)')
            title('Net eddy current on VV')
            set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
            Filename = 'I_VV';
            %Filename= sprintf('%s_Trampdown_%dms',Filename,2*T_ramp_Sol);    
            saveas(gcf, strcat(FigDir,Filename,FigExt));
  % }          
          
     %{     
    %HAVE TO CHECK Vp, is not closer to Vloop ==>?¿¿?¿?
    %HAVE TO DIG IN uFinal (x), because it has weird things
        %from 1 to Vessel(n) has eddys
        %last is Ip
        %But the rest are not IPF, are weird things==> ?¿¿?¿?¿?
        %No RIp,ZIp ?¿¿?¿?¿??¿?¿??¿
      
        
        %Plot Vp and Vloop
        
%         figure;
%             plot(time_adaptive*1e3,Vp_output)
%             ylabel('V (V)')
%            xlabel(' time (ms)')
%            hold on
%            plot(time_adaptive*1e3,V_loop*ones(1,length(time_adaptive)),'r*')
%            legend('Vp','Vloop')
%            set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%            title('Loop voltage and Vp (RZIp)')
%            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
 
%save_10ms=[time_adaptive Ip_output I_Passive_VV]
 %save('10ms_time_Ip_IVV','save_10ms')
%}
          
%%% END RZIP@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@BREAKDOWN CALC@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%%%Field%%%%%%%
%To compute the field, will do an equ calc, with no plasma, and the
%breakdown currents. It is vital to include eddy currents, so the
%fiesta_load_assembly will be used to include them.
%What VEST do on its article about startup is plooting the poloidal field
%just at the time when the ramp down is done, and a couple of ms after.
%Could do that. Beware of the plasma, since the issue of including it is
%not yet fulfilled, will try to choose time so that Ip is very little.

%To compute the field inside the vessel, provided that Fiesta returns it in
%the whole grid, have to create a grid for the vessel:

n_pnts_inside=100; %100 is the ideal to have good plots of the fields, but the 
    %L int failures. Need to be 30 for the integration.

r_inside_VV=linspace(VesselRMinInner,VesselRMaxInner,n_pnts_inside); 
z_inside_VV=linspace(VesselZMinInner,VesselZMaxInner,n_pnts_inside);

%Have to do ameshgrid for the interpolation, so each R value has a Z value
global R_in Z_in
[R_in,Z_in]=meshgrid(r_inside_VV,z_inside_VV);

%Some fiesta things are needed for the equil calc:
coilsetVV=fiesta_loadassembly(coilset, vessel);
configVV = fiesta_configuration( 'SMART with VV', Grid, coilsetVV);
  

%%%%Breakdown non optimised(only Sol)%%%%%%%%
 coil_currents_Break_non = zeros(1,nPF);
 coil_currents_Break_non(iSol) = I_Sol_Null; %current at the top previous to ramp down
% 
 icoil = fiesta_icoil( coilset, coil_currents_Break_non );
 equil_non_optimised = fiesta_equilibrium( 'Breakdown non optimised', config, Irod, icoil );
% 
% fileName = 'ST25D_non_optimised';
% legend(gca,'hide');
% % save_to_pdf( gcf, fileName );

%Extraction of the fields:
  [FieldsBreakNon_grid, FieldsBreakNon_VV]=fields(equil_non_optimised);
   
%Now lets plot the fields to deduce the direction of the Sol current   
%{
  %  % Plot Bz Br Bphi   
%     figure; 
%     subplot(1,6,1)
%     plot(equil_non_optimised);
%     
%     subplot(1,6,2)
%     contour(R_in,Z_in,FieldsBreakNon_VV.Br,1000);
%     shading('interp')
%     hold on;
%     plot(vessel)
%     plot(coilset)
%     c=colorbar; %colorbar
%       ylabel(c, 'Bz(T)');
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     %plot(sensor_btheta)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Br Break Non optimized (Sol only)')
%     %   % Plot Br   
%     subplot(1,6,3)
%     contour(R_in,Z_in,FieldsBreakNon_VV.Bz,1000);
%     shading('interp')
%     hold on;
%     plot(vessel)
%     plot(coilset)
%     c=colorbar; %colorbar
%       ylabel(c, 'Bz(T)');
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     %plot(sensor_btheta)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Bz Break Non optimized (Sol only)')
%     
%     subplot(1,6,4)
%     contour(R_in,Z_in,FieldsBreakNon_VV.Bphi,1000);
%     shading('interp')
%     hold on;
%     plot(vessel)
%     plot(coilset)
%     c=colorbar; %colorbar
%       ylabel(c, 'Bz(T)');
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     %plot(sensor_btheta)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Bphi Break Non optimized (Sol only)')
%     
%     subplot(1,6,5)
%     contour(R_in,Z_in,FieldsBreakNon_VV.Bpol,1000);
%     shading('interp')
%     hold on;
%     plot(vessel)
%     plot(coilset)
%     c=colorbar; %colorbar
%       ylabel(c, 'Bz(T)');
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     %plot(sensor_btheta)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Bpol Break Non optimized (Sol only)')
% 
%     subplot(1,6,6)
%     quiver(R_in,Z_in,FieldsBreakNon_VV.Br,FieldsBreakNon_VV.Bz,'AutoScale','on','AutoScaleFactor', 10)
%     hold on;
%     plot(vessel)
%     plot(coilset)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Quiver plot Bpol Break Non optimized (Sol only)')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%So, what I will do is a loop with the time from RZIp, to define the
%current values. The time steps will be
time_loop=[0]*1e-3 %[s] time for calculate things for breakdown
Ip_loop=interpn(time_adaptive,Ip_output,time_loop)  %[A] the plasma current when calculated
                                                                                    %the breakdwon
 IPassive_loop=interpn(time_adaptive,I_Passive_VV*1e-3,time_loop)  %[kA] total eddy
 
%A very rough stimate of the poloidal field of the plasma could be done by
%considering the plasma a circular wire with radius Rgeo. I have used a
%function to compute the field of a circular wire with Biot Savart, to see
%when the plasma field is relevant. 

%Will also compute the 'potential function' of Lazarus to see where the
%breakdown occurs. I divide by the Loop voltage since we do not know it yet
   
for loop=1:length(time_loop)
    
    %Current of the coilset
    coil_currents_break = zeros(1,nPF+get(vessel,'n'));                       %n=308, number of filaments
    coil_currents_break(iSol) = interpn(time_adaptive,I_PF_output(:,iSol),time_loop(loop));
    coil_currents_break(iPF1) =interpn(time_adaptive,I_PF_output(:,iPF1),time_loop(loop));   
    coil_currents_break(iPF2) =interpn(time_adaptive,I_PF_output(:,iPF2),time_loop(loop));   
    coil_currents_break(iDiv1) =interpn(time_adaptive,I_PF_output(:,iDiv1),time_loop(loop));   
    coil_currents_break(iDiv2) =interpn(time_adaptive,I_PF_output(:,iDiv2),time_loop(loop)); 

    %Passive current
    %Fuck, have to do a loop to define the current of each filament :(
    
    for i=1:get(vessel,'n')
        coil_currents_break(nPF+i)=interpn(time_adaptive,I_Passive(:,i),time_loop(loop));   %eddy currents, have to define each filament
        %coil_currents_break(nPF+i)=0; %NULLING EDDY CURRENTS
    end
    
    %With the currents, the icoil object can be created, and with it the
    %equil
    icoilVV_break = fiesta_icoil( coilsetVV, coil_currents_break );
    equil_break = fiesta_equilibrium( 'optimised null (eddys in)', configVV, Irod, icoilVV_break );
 
    %Manual plot
    psi_null=get(equil_break,'Psi_vac'); %fiesta field


    %
    psi_null=get(psi_null,'data','2D'); %this retunrs 2D data, so no need to reshape!
                    %251*200
    %{  
% PLOT WHOLE VV
%  RGrid2D=get(Grid,'R','2D'); %251*200
%  ZGrid2D=get(Grid,'Z','2D'); %251*200
% %        
%  figure;
%     contour(RGrid2D,ZGrid2D,log10(psi_null),5000)   
%     %shading('interp') %this is to make the transition between values continuous,
%          %instedad of discontinuously between pixels
%     hold on
%     plot(vessel)
%     plot(coilset) 
%     %colormap(Gamma_II)
%     plot(vessel)
%     colorbar %colorbar
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%       [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
%      xlabel('R (m)')
%      ylabel('Z (m)')
%      title(sprintf('log10(Psi)  at t=%d ms (iter %d/%d) 5000c',time_loop(loop)*1e3,loop,length(time_loop)))
%}      
  
psi_null_interpn=  @(r,z) interpn(zGrid,rGrid,psi_null,z,r,'mikama');        
psi_null_ins_VV=psi_null_interpn(R_in,Z_in);    %contour plot!!!
    %Plot
        figure; 
        contour(R_in,Z_in,psi_null_ins_VV,500)
        hold on;
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'Psi(W)');  
                %Try to show max and min of colorbar !!!!!!!!!!!!!
        %t=get(c,'Limits');
        %set(c,'Ticks',linspace(t(1),t(2),10));
        view(2) %2D view
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
        [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        axis equal
        xlabel('R (m)')
        ylabel('Z (m)')
        title(sprintf('Psi  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))
        %title(sprintf('B_{pol} t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'Psi';
        %Filename= sprintf('%s_simu_%d',Filename,sen);      
        saveas(gcf, strcat(FigDir,Filename,FigExt));
     
 
%Extraction of the fields (Earth inside)!
    [FieldsBreak, FieldsBreakNoEarth]=fields(equil_break); %To include Earths!!!
    %[thisDoIncludeEarths,FieldsBreak]=fields(equil_break); %To NOT include Earths!!!   
    
%Plots!

  %Bpol 
        figure; 
        %contourf(R_in,Z_in,log10(abs(Bpol_ins_vessel)),'EdgeColor','none');
        %surf(R_in,Z_in,log10(abs(Bpol_ins_vessel)),'EdgeColor','none'); shading('interp') %this is to make the transition between values continuous,
        %instedad of discontinuously between pixels
        contourf(R_in,Z_in,log10(FieldsBreak.VV.Bpol),'ShowText','on')
        hold on;
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'log10(Bpol (T))');  
                %Try to show max and min of colorbar !!!!!!!!!!!!!
        %t=get(c,'Limits');
        %set(c,'Ticks',linspace(t(1),t(2),10));
        view(2) %2D view
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
        [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        axis equal
        xlabel('R (m)')
        ylabel('Z (m)')
        title(sprintf('B_{pol}  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))
        %title(sprintf('B_{pol} t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'Bpol';
        %Filename= sprintf('%s_simu_%d',Filename,sen);      
        saveas(gcf, strcat(FigDir,Filename,FigExt));

%Bphi
        figure; 
        contourf(R_in,Z_in,log10(FieldsBreak.VV.Bphi),'ShowText','on')
        hold on;
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'log10(Bphi (T))');  
                %Try to show max and min of colorbar !!!!!!!!!!!!!
        %t=get(c,'Limits');
        %set(c,'Ticks',linspace(t(1),t(2),10));
        view(2) %2D view
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
        [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        axis equal
        xlabel('R (m)')
        ylabel('Z (m)')
        title(sprintf('B_{phi}  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))
        %title(sprintf('B_{phi} t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'Bphi';
        %Filename= sprintf('%s_simu_%d',Filename,sen);      
        saveas(gcf, strcat(FigDir,Filename,FigExt));        
        
   %Lloyd
           
     Lloyd=E(R_in).*FieldsBreak.VV.Bphi./FieldsBreak.VV.Bpol;
         figure;
        contourf(R_in,Z_in,log10(Lloyd),'ShowText','on')
        hold on
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'log10(E*Bphi/Bpol (V/m))');
        view(2) %2D view
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        xlabel('R (m)')
        ylabel('Z (m)')
        %axis([0,1,0,1])        %For reduced size and compare to VEST plot
        title(sprintf('Lloyd criteria  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        %title(sprintf('Lloyd criteria t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'LLoyd';
        %Filename= sprintf('%s_simu_%d',Filename,sen);    
        saveas(gcf, strcat(FigDir,Filename,FigExt));
               
   %%L CALC FORMULAE%%%%%%%%%%%%

    %The field in the whole vessel will be used, since it is the field that is
    %used in the field line integrator.

    %Minimun with veesel interp 
    Bpolmin=min(min(FieldsBreak.sensor.Bpol));                                         %minimun of poloidal field
    [Bpolmin_index_row Bpolmin_index_column]=find(FieldsBreak.sensor.Bpol==Bpolmin);     %indexes  

        %%%%%THIS IS TO COMPARE

        
    %And the toroidal field at the center of the null region is:
    Bphi_centerNull=FieldsBreak.interpn.Bphi((min(r_sensors)+max(r_sensors))/2,...
    (min(z_sensors)+max(z_sensors))/2);

%     %The connective length without averaging is
%     L_no_aver=0.25*a_eff*Bphi_centerNull/Bpolmin                            %[m]

    %A more precise way to calculate L is by averaging in the poloidal field.
    %As stated by TCV thesis and the newer ITER article, Bpol is obtained by
    %averaging on the surface of the null region. The toroidal field is not
    %averaged. To do this, if I use the sensor field, it will be much much more
    %easier, so I will do it.

    
    Bpolmin_av=mean(mean(FieldsBreak.sensor.Bpol));              %to compute the mean inside the sensor region
    
    a_sensor=(max(get(sensor_btheta,'r'))-min(get(sensor_btheta,'r')))/2; %[m]minor radius of the null region
    R_major_sensor=min(get(sensor_btheta,'r'))+a_sensor;            %[m] Major raidus  of the sensorregion
    L_aver(loop)=0.25*a_sensor*Bphi_centerNull/Bpolmin_av                %[m] L with the average pol field
    Campos_L(loop)=Bphi_centerNull/Bpolmin_av; %field to compute L emp
    
    L_emp_Scott=0.25*R_major_sensor*Bphi_centerNull/Bpolmin_av   
    
    %REMEMBER THAT THIS IS COMPUTED INSIDE NULL REGION, WHILE THE LOWER
    %BPOL IS NOT IN THAT REGION DUE TO EDDYS (AND EARTH)!!!!!!!!!!!!!!!!!!
    %%NEED TO THINK ON THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    %% L calc by field line integration###########################
% { 
         %Lazarus paper 1998, they compute connective length by avergaing on
            %9 lines, the line with Bpol min, and the 8 surroundings. 
         %However can compute the lines in all the VV
      
   IntMethod='Lp';         % ''Phi'or 'Lp' to switch the integration mode 

        %Grid
        %I redefine the grid to compute the connection length, for less computer
        %demands (time). Will label it with an L at the end!

        n_pnts_insideL=30; %30 previously          %100 is the ideal to have good plots of the fields, but the L int failures. 

        r_inside_VVL=linspace(VesselRMinInner,VesselRMaxInner,n_pnts_insideL); 
        z_inside_VVL=linspace(VesselZMinInner,VesselZMaxInner,n_pnts_insideL);

        %Lets do a meshgrid, will be needed
        [r_ins_VVL,z_ins_VVL]=meshgrid(r_inside_VVL,z_inside_VVL);
        
        %Plot
            figure;
            plot(r_ins_VVL,z_ins_VVL,'r.')
            hold on
            plot(vessel)
            axis equal
            xlabel('R (m)')
            ylabel('Z (m)')
            title(sprintf('meshgrid for the integration with %d^2 points',n_pnts_insideL))

        %To remove the points at the VV:
         %r_inside_VVL=r_inside_VVL(2:end-1);
         %z_inside_VVL=z_inside_VVL(2:end-1); 
        
       %Now, if the coils are inside, the grid points there have to be removed because ode 'se raya', and spent too long time
       
       Rin_coils=[R_PF2 R_Div1 R_Div2]-width_PF/2*ones(1,nPF-2); %inner R of coilset (no SOl and PF1 (outside))
       Rout_coils=[R_PF2 R_Div1 R_Div2]+width_PF/2*ones(1,nPF-2); %outer R of coilset (no SOl and PF1 (outside))
       Zdown_coils=[Z_PF2 Z_Div1 Z_Div2]-height_PF/2*ones(1,nPF-2); %lowerZ of coilset (no SOl and PF1 (outside))
       Zup_coils=[Z_PF2 Z_Div1 Z_Div2]+height_PF/2*ones(1,nPF-2); %higher Z of coilset (no SOl and PF1 (outside))    
            %They are good
        
       %Lets rehape the meshgrid to do the loop to remove points
       r_ins_VVL=reshape(r_ins_VVL,length(r_ins_VVL)^2,1);
       z_ins_VVL=reshape(z_ins_VVL,length(z_ins_VVL)^2,1);
       
    
        for co=1:length(Rin_coils) %at each iter, removes the grid points inside the coils
              StoreRZ=[0 0]; %initialization of stored grid points
            for i=1:length(r_ins_VVL) %Have to check each point
                Point=[r_ins_VVL(i) z_ins_VVL(i)]; %grid point to test            
                
                switch sign(Point(2)) %First lets check if Z><0
                    
                    case 1 %z>0
                
                    if Point(1)<Rin_coils(co) | Point(1)>Rout_coils(co) %R out of the coil
                                                                                                   %All Z are good
                               StoreRZ=[StoreRZ; Point]; %store of good points                                                                                           
                    
                    elseif  Point(1)>Rin_coils(co) | Point(1)<Rout_coils(co) %R inside of the coil
                            if Point(2)<Zdown_coils(co) | Point(2)>Zup_coils(co)  %Z out coil
                                StoreRZ=[StoreRZ; Point]; %store of good points
                            end
                    end
                        
                    case -1 %z<0
                
                    if Point(1)<Rin_coils(co) | Point(1)>Rout_coils(co) %R out of the coil
                                                                                                   %All Z are good
                               StoreRZ=[StoreRZ; Point]; %store of good points                                                                                        
                    
                    elseif  Point(1)>Rin_coils(co) | Point(1)<Rout_coils(co) %R inside of the coil
                            if Point(2)>-Zdown_coils(co) | Point(2)<-Zup_coils(co)  %Z out coil
                                StoreRZ=[StoreRZ; Point]; %store of good points
                            end
                    end  
                        
                end             
            end  
            StoreRZ=StoreRZ(2:end,:); %remove first row, the initialization one
            r_ins_VVL=StoreRZ(:,1); %to use the grid created on the following coil loop
            z_ins_VVL=StoreRZ(:,2); %to use the grid created on the following coil loop          
        end             
     
     %Lets reshape it again to do the contour plots later (both are
     %vectors)
     r_ins_VVL_contour=reshape(r_ins_VVL,floor(length(r_ins_VVL)/2),[]); %arbitrary reshape
     z_ins_VVL_contour=reshape(z_ins_VVL,floor(length(z_ins_VVL)/2),[]); %arbitrary reshape
                
        %Plot
                figure;
                subplot(1,2,1)
                plot(StoreRZ(:,1),StoreRZ(:,2),'r.')
                hold on
                plot(vessel)
                axis equal
                title(sprintf('iter %d',co))   
                subplot(1,2,2)
                plot(StoreRZ(:,1),StoreRZ(:,2),'r.')
                hold on
                plot(vessel)
                plot(coilset)
                axis equal
                title(sprintf('iter %d',co))    
  
    %%%End grid
    
    %%%Begin integration
    
    L_max=1000;                             %[m] max L value for the integration; when L achieves
                                                            %this value, the integration stops. 
                                                            %Iter 55,85,86,81 achieves about 10,000, spending 4h
                                     
    event_colission_wall=@(t,y) Colission_wall(t,y,VesselRMaxInner,...
            VesselRMinInner,VesselZMaxInner,VesselZMinInner,FieldsBreak.interpn.Br,...
            FieldsBreak.interpn.Bz,FieldsBreak.interpn.Bphi,L_max); 
    options = odeset('OutputFcn',@ode_progress_bar,'Events',event_colission_wall,'AbsTol',1e-10,'RelTol',1e-6); 
                                    %I include a Fiesta funciton to show the progress of the ode
  
   %Both integrators are programmed, so to swich between them will use a
   %switch (xD)
   tic          %to measure time the intergration takes
    switch IntMethod %switch for the starting points (change the number of inputs)
            
         case  'Lp' %Phi as integration method
              %1) Starting points y0
                 y0=[0 0 0 0 0];                    %(r0,z0,L0,phi0,U0) starting points
                    %note it has to be r z L for using the same event function

                    for i=1:length(r_ins_VVL)        
                         points=[r_ins_VVL(i) z_ins_VVL(i) 0 0 0];      %r z L phi U        
                            %U(0)=0 (arbitrary)  
                         y0=[ y0; points];                          
                    end 
                %I have the additional point 0 0 0 form the begining, that i can remove
                %easily with
                y0=y0(2:end,:);
                
               %2)Independt variable values, t0
                 t_values=1000; %1000            %Max value of the independent variable
                 n_iter_t=30000000; %3000000 on s1-14                 %Integer, number of values for tSpan
                 tSpan=linspace(0,t_values,n_iter_t);            %the range of values of independant variable
                                     
                %3)Odefun
                 odefun= @(Lp, rzLphiU) Field_LineIntegrator_Lp(Lp,rzLphiU,FieldsBreak.interpn.Br,...
                    FieldsBreak.interpn.Bz,FieldsBreak.interpn.Bphi);                
                
                %4) Integration
                        %%%%%%%%SINGLE FIELD LINE TRACER and plotter

                        %need to find i for the chosen R,Z value in r0_z0_L0_U0.
                        %I= 85 for a line inside, 49 for a max L outside, 152 for the
                        %outward arm (Z>0). 135 for the outward Z<0 line. 64 for the upper
                        %arm
        
                        i=49 %looked in the y0
                        [t_fieldline, y_fieldline]=ode45(odefun,tSpan,y0(i,:),options);    
    
                        %To save the last values of R,Z,L
                        L_single=y_fieldline(end,3);         %L      
        
                        %Plot of one of the line
                            figure;
                            plot3(y0(i,1)*cos(y0(i,4)),y0(i,1)*sin(y0(i,4)),y0(i,2),'k*','LineWidth',3)
                            hold on;
                            plot3(vessel);
                            plot3(coilset);
                            plot3(y_fieldline(:,1).*cos(y_fieldline(:,4)),y_fieldline(:,1).*sin(y_fieldline(:,4)),...
                                y_fieldline(:,2),'r.','LineWidth',3)
                                xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');  
                            plot3(y_fieldline(end,1).*cos(y_fieldline(end,4)),y_fieldline(end,1).*sin(y_fieldline(end,4)),...
                                y_fieldline(end,2),'g*','LineWidth',3)
                            title(sprintf('Single field line integration Lp, L=%3.2f m ',L_single))
                            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
                            %legend('Starting point (Point with less Bpol)','Field line',...
                        %%%END ONE LINE TRACER
   
                        for i=1:length(y0)
                            fprintf('Iter %d de %d',i,length(y0))
                            [t_fieldline, y_fieldline]=ode45(odefun,tSpan,y0(i,:),options);        %ode15s Carlos
    
                            %To save the last values of the dependants variables
                            y_end(i,1)=y_fieldline(end,1);           %R
                            y_end(i,2)=y_fieldline(end,2);           %Z
                            y_end(i,3)=y_fieldline(end,3);           %L           
                            y_end(i,4)=y_fieldline(end,4);           %Phi
                            y_end(i,5)=y_fieldline(end,5);           %U
                        end                 
                        U_int=reshape(y_end(:,5),size(r_ins_VVL_contour,1),size(r_ins_VVL_contour,2));            
        
        case 'Phi' %Lp as integration method
              %1) Starting points y0
                y0=[0 0 0 0]; %r0,z0,L0,U0
                %note it has to be r z fro using the same event function

                for i=1:length(r_ins_VVL)        
                        points=[r_ins_VVL(i) z_ins_VVL(i) 0 0];  %r z L phi U        
                            %U(0)=0 (arbitrary)  
                        y0=[ y0; points];                          
                end
                %I have the additional point 0 0 0 form the begining, that i can remove
                %easily with
                y0=y0(2:end,:);
            
               %2)Independt variable values, t0
                 t_values=1000; %100 to little                            %Cycles(toroidal turns)
                 n_iter_t=300000;         %3000                         %Integer, number of values for tSpan
                 tSpan=linspace(0,2*pi*t_values,n_iter_t);    %the range of values of independant variable
               
               %3)Odefun 
                 odefun= @(phi, rzLU) Field_LineIntegrator(phi,rzLU,FieldsBreak.interpn.Br,...
                    FieldsBreak.interpn.Bz,FieldsBreak.interpn.Bphi);               
                %4) Integration
                        %%%%%%%%SINGLE FIELD LINE TRACER and plotter

                        %need to find i for the chosen R,Z value in r0_z0_L0_U0.
                        %I= 85 for a line inside, 49 for a max L outside, 152 for the
                        %outward arm (Z>0). 135 for the outward Z<0 line. 64 for the upper
                        %arm
        
                        i=49 %looked in the y0
                        [t_fieldline, y_fieldline]=ode45(odefun,tSpan,y0(i,:),options);    
    
                        %To save the last values of R,Z,L
                        L_single=y_fieldline(end,3);         %L      
        
                        %Plot of one of the line
                        figure;
                        plot3(y0(i,1)*cos(t_fieldline(1)),y0(i,1)*sin(t_fieldline(1)),y0(i,2),'k*','LineWidth',3)
                        hold on;
                        plot3(vessel);
                        plot3(coilset);
                        plot3(y_fieldline(:,1).*cos(t_fieldline(:)),y_fieldline(:,1).*sin(t_fieldline(:)),...
                            y_fieldline(:,2),'r.','LineWidth',3)
                        xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');  
                        plot3(y_fieldline(end,1).*cos(t_fieldline(end)),y_fieldline(end,1).*sin(t_fieldline(end)),...
                            y_fieldline(end,2),'g*','LineWidth',3)
                        title(sprintf('Single field line integration phi, L=%3.2f m ',L_single))
                        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
                        %%%END ONE LINE TRACER
                
                        for i=1:length(y0)
                            fprintf('Iter %d de %d',i,length(y0))
                            [t_fieldline, y_fieldline]=ode45(odefun,tSpan,y0(i,:),options);        %ode15s Carlos
    
                            %To save the last values of the dependants variables
                            y_end(i,1)=y_fieldline(end,1);           %R
                            y_end(i,2)=y_fieldline(end,2);           %Z
                            y_end(i,3)=y_fieldline(end,3);           %L           
                            y_end(i,4)=y_fieldline(end,4);           %U
                        end 
                        U_int=reshape(y_end(:,4),size(r_ins_VVL_contour,1),size(r_ins_VVL_contour,2));
        
    end         
   Time_Run_Integrator=toc     %Time spent by the integrator
   
   L_int=reshape(y_end(:,3),size(r_ins_VVL_contour,1),size(r_ins_VVL_contour,2));
    
   %Calc of average L on the null region
        index= y0(:,1)<=max(get(sensor_btheta,'r')) & y0(:,1)>=min(get(sensor_btheta,'r')) &...
            y0(:,2)<=max(get(sensor_btheta,'z')) & y0(:,2)>=min(get(sensor_btheta,'z')); 
                    %Index of R,Z inside null
        L_int_row=y_end(:,3); %[m] L in row form
        L_int_null=mean(L_int_row(index))
   
   %%%Storing of the non colliding starting points
    %To store start points that do not collide: first I get the index of both R
    %and Z, but together, since they do not collide if oth R and Z are greater
    %than the min value, and lower than the greatest value
    
    RZ_store_index=y_end(:,1)<VesselRMaxInner & ...
         y_end(:,1)>VesselRMinInner & y_end(:,2)<VesselZMaxInner...
         & y_end(:,2)>VesselZMinInner; %100*5==> error, has to ve vector,, not matrix!!!
                        
    RZ_no_collide=[y0(RZ_store_index,1) y0(RZ_store_index,2)];    
           
    %However, this is not perfect, when including in the grid the points in the wall,
    %something extrange happens, some of them are store in the non colliding points
    %thought they do not collide since the starting point is also the ending points
    %(you get like some stars just in the VV, but not all, only some of them)
    %To remove it:
     Index=RZ_no_collide(:,1)<VesselRMaxInner & ...
          RZ_no_collide(:,1)>VesselRMinInner & RZ_no_collide(:,2)<VesselZMaxInner...
          & RZ_no_collide(:,2)>VesselZMinInner;
          RZ_no_collide=[RZ_no_collide(Index,1) RZ_no_collide(Index,2)];
                   
       
     %%%Contour plots
      %1) L
        figure;
        contourf(r_ins_VVL_contour,z_ins_VVL_contour,L_int,'ShowText','On')
        %surf(r_ins_VVL_contour,z_ins_VVL_contour,L_int,'EdgeColor','none'), shading('interp')
        view(2)
        hold on
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'g.--')
        plot(RZ_no_collide(:,1),RZ_no_collide(:,2),'m*')
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        colormap(Gamma_II)
        c=colorbar; %colorbar
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        ylabel(c, 'L(m)');
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('L','Field null region','Non colliding points')
        title(sprintf('L  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'L_cont';
        %Filename= sprintf('%s_simu_%d',Filename,sen);     
        saveas(gcf, strcat(FigDir,Filename,FigExt));       

        %1D to 2D plot
        figure;
        tri = delaunay(y0(:,1),y0(:,2));
        plot(y0(:,1),y0(:,2),'.')
        [r,c] = size(tri); %number of triangles there
        trisurf(tri, y0(:,1), y0(:,2),y_end(:,3),'FaceAlpha',1,'EdgeColor','none'), shading('interp');
        view(2)
        colorbar
        hold on                     %this is to make the transition between values continuous,                                                        %instedad of discontinuously between pixels
        colormap(Gamma_II)
        plot3([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],...
            ones(1,5)*max(y_end(:,3)),'g.--')
        plot3(RZ_no_collide(:,1),RZ_no_collide(:,2),max(y_end(:,3))*ones(length(RZ_no_collide(:,2)),1),'m*')
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        colormap(Gamma_II)
        c=colorbar; %colorbar
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        ylabel(c, 'L(m)');
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('L','Field null region','Non colliding points')
        title(sprintf('L  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'L';        
        saveas(gcf, strcat(FigDir,Filename,FigExt));  
        
        figure;
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'g.--')
        hold on
        plot(RZ_no_collide(:,1),RZ_no_collide(:,2),'m*') 
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);        
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');      
        colormap(Gamma_II)
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('Field null region','Non colliding points')
        title(sprintf('L  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'L_noL';                
        saveas(gcf, strcat(FigDir,Filename,FigExt));  
        
      %2) Pseudo potential U/Vloop
        figure;
        contourf(r_ins_VVL_contour,z_ins_VVL_contour,U_int,10)
        %surf(r_insVV_noLimit,z_insVV_noLimit,U_int), shading('interp')
        hold on
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'U/Vloop');
        view(2) %2D view
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        xlabel('R (m)')
        ylabel('Z (m)')
        title(sprintf('Pseudo potential  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        %title(sprintf('Pseudo potential at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'Pseudo_contour';
        %Filename= sprintf('%s_simu_%d',Filename,sen);     
        saveas(gcf, strcat(FigDir,Filename,FigExt));
        
        %1D to 2D plot
                %1D to 2D plot
        figure;
        tri = delaunay(y0(:,1),y0(:,2));
        plot(y0(:,1),y0(:,2),'.')
        [r,c] = size(tri); %number of triangles there
        trisurf(tri, y0(:,1), y0(:,2),y_end(:,5),'FaceAlpha',1,'EdgeColor','none'), shading('interp');
        view(2)
        colorbar
        hold on                     %this is to make the transition between values continuous,                                                        %instedad of discontinuously between pixels
        colormap(Gamma_II)
        plot3([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],...
            ones(1,5)*max(y_end(:,3)),'g.--')
        plot3(RZ_no_collide(:,1),RZ_no_collide(:,2),max(y_end(:,3))*ones(length(RZ_no_collide(:,2)),1),'m*')
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        colormap(Gamma_II)
        c=colorbar; %colorbar
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        ylabel(c, 'U/Vloop');
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('U/Vloop','Field null region','Non colliding points')
        title(sprintf('U/Vloop  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'Pseudo';        
        saveas(gcf, strcat(FigDir,Filename,FigExt));  
        
        
      %%%3)[Experimental] E_rel plot, to predict where the gas breaks down
        
        p_test=5*10^-5; %Tor
        E_RZmin=C_2*p_test./(log(C_1*p_test*L_int)); %E min, Paschen, but 2D        
        E_RZmin(E_RZmin<0)=NaN; %when Emin<0, there is no breakdwon, so NaN not
                    %to plot it
        
        figure;
        contourf(r_ins_VVL_contour,z_ins_VVL_contour,U_int./L_int*V_loop./E_RZmin)
        %surf(r_insVV_noLimit,z_insVV_noLimit,U_int), shading('interp')
        hold on
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'E_{rel}');
        view(2) %2D view
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        xlabel('R (m)')
        ylabel('Z (m)')
        %title(sprintf('Pseudo potential  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        %title(sprintf('E_{rel} at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        Filename = 'U_L';
        %Filename= sprintf('%s_simu_%d',Filename,sen);   
        saveas(gcf, strcat(FigDir,Filename,FigExt));
% }

end    %end loop time (not in use right now)!!!!

        %{ 
Plot Vloop
  
%   V_loop=@(t_rd) pi*RSol^2*mu0*turns(iSol)/(2*ZMaxSol)*2*I_Sol_max./t_rd;
%         
%   figure;
%   plot(5*1e-3,V_loop(5*1e-3),'*','MarkerSize',10)
%   hold on
%   plot(10*1e-3,V_loop(10*1e-3),'*','MarkerSize',10)
%   plot(20*1e-3,V_loop(20*1e-3),'*','MarkerSize',10)
%   plot(50*1e-3,V_loop(50*1e-3),'*','MarkerSize',10)
%   plot(linspace(5,50)*1e-3,V_loop(linspace(5,50)*1e-3),'k')
%   legend('5ms','10ms','20ms','50ms')
%   xlabel('t_{rd} (ms)')
%   ylabel('V_{loop} (V)')      
%     title('V_{loop} at t=0ms vs ramp-down time')
%     set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%         Filename = 'V_loop_vs_ramp_down'; 
%         saveas(gcf, strcat(FigDir,Filename,FigExt));                
%}        


%% %Earths field:
        [Field_EarthGrid Field_Earth]=EarthField(R_in,Z_in)
        
     
%%%%End test plots field line%%%%%%%
    
        
    figure;
     subplot(1,3,1)
    contour(R_in,Z_in,FieldsBreak.VV.Bz,1000);
    shading('interp')
    hold on;
    plot(vessel)
    plot(coilset)
    c=colorbar; %colorbar
      ylabel(c, 'Bz(T)');
    view(2) %2D view
    plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
    [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
    %plot(sensor_btheta)
    xlabel('R (m)')
    ylabel('Z (m)')
title(sprintf('B_z at t=%d ms (iter %d/%d) 1000c',time_loop(loop)*1e3,loop,length(time_loop)))
    %   % Plot Br   
	subplot(1,3,2)
    contour(R_in,Z_in,FieldsBreak.VV.Br,2000);
    shading('interp')
    hold on;
    plot(vessel)
    plot(coilset)
    c=colorbar; %colorbar
      ylabel(c, 'Br(T)');
    view(2) %2D view
    plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
    [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
    %plot(sensor_btheta)
    xlabel('R (m)')
    ylabel('Z (m)')
title(sprintf('B_r  at t=%d ms (iter %d/%d) 1000c',time_loop(loop)*1e3,loop,length(time_loop)))
    
subplot(1,3,3)
        contour(R_in,Z_in,FieldsBreak.VV.Bphi,1000);
    shading('interp') %this is to make the transition between values continuous,
    %instedad of discontinuously between pixels
    colormap(Gamma_II)
    hold on;
    plot(vessel)
    plot(coilset)
    c=colorbar; %colorbar
      ylabel(c, 'Bphi(T)');
    view(2) %2D view
    plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
     xlabel('R (m)')
    ylabel('Z (m)')
     title(sprintf('B_phi at t=%d ms (iter %d/%d) 1000c',time_loop(loop)*1e3,loop,length(time_loop)))
      
 %%%%%%End plots B and L%%%%%%%%%    
        

%% ADDITIONAL FEATURES (EDDYS, FORCES, INTRODUCING EDDYS ON EQUIL (TRY), ETC


% %% RE-DOING EQUILIBRIA CALC WITH EDDYS (FAIL)
%{
    %Now that we havr compute the eddys, we could re do all the calc, to
    %get an equil with those eddys, with the new equil do the RZip again,
    %and, it the new eddys do not change much, could accept the result.
    %Yes, this is a iterative process, but with a huge consumption of
    %computer resources==> :)
     
     %coil currents for the new equilibria calc
     
    coil_currents_eddys = zeros(1,nPF+get(vessel,'n'));                       %n=308, number of filaments
    coil_currents_eddys(iSol) =	CoilWaveforms(iSol,5);             %Do not change, in principle
    coil_currents_eddys(iPF1) =CoilWaveforms(iPF1,5);  
    coil_currents_eddys(iPF2) =CoilWaveforms(iPF2,5);  
    coil_currents_eddys(iDiv1) =CoilWaveforms(iDiv1,5);  
    coil_currents_eddys(iDiv2) =CoilWaveforms(iDiv2,5);            %May need to be changed

    %To select the eddy currents, provided that at the equil time the eddys are amx, can
    %just look at the max:
    
    [max_eddy,index_max_eddy]=max(I_Passive_VV);                %the total eddy, the max value
    coil_currents_eddys(nPF+1:end)=I_Passive(index_max_eddy,:);       %eddy currents, have to define each filament
    
    icoil_eddy = fiesta_icoil( coilsetVV, coil_currents_eddys );
    
    %%%1) EFIT
    [efit_configVV, signalsVV, weightsVV, indexVV]=efit_shape_controller(configVV, {'PF1','PF2'}, [0.44, 0, 0.44/1.85 1.8 0.2])
    % The numbers you give are [Rgeo, Zgeo, a, kappa, delta], Rgeo,Zgeo,a are
    % mandatory.
    %I use the values of the standar shape, to get a similar equil

    equil_eddy=fiesta_equilibrium('Target+eddys', configVV, Irod, jprofile, control,efit_configVV, icoil_eddy, signalsVV, weightsVV) %%EFIT!!!
    %It does the case in line 96!! The equil calc is in lin 124

    %Now we have to extract the new currents from the equil, provided that EFIT
    %changed some of them to satisfy the conditions requested:
    icoil_eddy=get(equil_eddy,'icoil');
    current_post_EFIT=get(icoil_eddy,'currents');
    coil_currents_eddys(iPF1) =current_post_EFIT(iPF1);
    coil_currents_eddys(iPF2) =current_post_EFIT(iPF2);
    coil_currents_eddys(iDiv1) =current_post_EFIT(iDiv1);
    coil_currents_eddys(iDiv2) =current_post_EFIT(iDiv2);
%No need of redefine the Sol current of course. Actually, Div 1 and 2 are
%not neccessary, since I am not changing them


%     %%2) NO EFIT
%             equil = fiesta_equilibrium( 'STV2C2', config, Irod, jprofile, control, [],icoil );

                %This is not an option, since it does not
                %convnerge with the values without eddys!!!!!!!!!!
                             
  %Plot of the equil with the eddys!
        %section_figure=section(equil); %THIS IS A PLOT
        figure;
        plot(equil_eddy)
        parametersshow(equil_eddy) %this plots the parameters in the equil
        hold on
        plot(vessel)
        plot(coilset)
        title('Target equilibria eddys included')


%%% Make virtual sensors where we want breakdown  %%%%%%%

% %Initiate virtual B-field sensors centered on Rgeo
% sensor_btheta = InitiateBSensors(param_equil,a_eff);
% 
% %(r,z) of the sensors
% 
% r_sensors=get(sensor_btheta,'r'); %size 1*200
% z_sensors=get(sensor_btheta,'z'); %size 1*200
% 
% global R_sensor Z_sensor
% [R_sensor,Z_sensor]=meshgrid(r_sensors,z_sensors); %size 200*200
  
%Plot
%         figure;
%         plot(equil)        
%         hold on
%         plot(vessel)
%         plot(coilset)
%         parametersshow(equil)   %this plots the parameters in the equil
%         title('Target equilibria')
%         plot(sensor_btheta);

%%%%END OF FIESTA EQ WITH EDDYS@@@@@@@@@@@@@@@@@@@@@@@@@@


%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@RZIp WITH EDDYS@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% It has the virtual sensors

%Redefinition of the icoil in the equil with eddys to not have the vessel

% equil_eddy=set(equil_eddy,configVV,'icoil',fiesta_icoil( coilset, coil_currents_eddys(iSol:iDiv2)));
% 
% %%%%DO NOT WORK :)
% 
% %%%%
% rzip_configVV = fiesta_rzip_configuration( 'RZIP with eddys', configVV, vessel, {sensor_btheta} );
%     %have to use the confing without the vessel, which makes sense
% [A, B, C, D, curlyM, curlyR, gamma, plasma_parameters, index, label_index, state] = response(rzip_configVV, equil_eddy, 'rp',plasma_resistance);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%  Optimised null  %%%%%%%%%%%%
% %Copy from ST25D stimation 
% C_temp = C(end-get(sensor_btheta,'n')+1:end,1:nPF);
% C1 = C_temp(:,1);
% D1 = C_temp(:,2:end);

%%%%%THIS DO NOT WORK, SINCE THE COILSET ALSO CONTAINS THE VESSEL, SO
%%%%%RESPONDE DO NOT WORK PROPERLY (eig(A) contains NaN or Inf!!!!!!!!!!!!!!!!!!!!!!ç
%}

%%
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@FORCES AND EDDY CURRENT PLOTS@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%% EDDY CURRENTS ON THE VESSEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the filament variables r and z
ptmp = get(vessel,'passives');
ftmp = get(ptmp,'filaments');
Rfil = get(ftmp(:),'r');        %dim 1*number of filaments; number of filaments=nf
Zfil = get(ftmp(:),'z');        %dim 1*number of filaments
%Note Rfil and Zfil are xaccum and yaccum, but transposed!!

%dim I_Passive= 3811*nf, and time adaptative is
%3811*1. So, I_Passive contains the eddy current of the nf filaments at
%each instant of time.

%Okay. We do have the eddy current on each filament an on every time on
%I_Passive. I_passive on the figures is sum over each filament of the eddy
%current, to get the total eddy current induced upon each time. We can not
%sum for every time the eddy current, since it varys its sign, it also
%ceases during certain amounts of tim, so it can not be done.  But there is
%no necesssity, since I_passive contains the eddy current at any time of
%the filament, and this will also provide the force upon each instant; i
%only need to pick up the greatest

sizeIpas=size(I_Passive); %number of time intervals*number of filaments

for i=1:sizeIpas(2) % for each filament
    %To decide if the largest value is positive or engative, i could chose
    %the max and the min, and compare its absolute value, and store the
    %greates
 
    positive=max(I_Passive(:,i));
    negative=min(I_Passive(:,i));
    
    if abs(positive)> abs(negative)
        
  I_Passive_fil(i)=positive;
    else
        I_Passive_fil(i)=negative;
    end
end
%This have just stored the largest values of the eddy currents on all the
%filaments.

%Plot
    figure;
    scatter3(Rfil,Zfil,I_Passive_fil/(1e3),100,I_Passive_fil/(1e3),'filled')
    hold on
    plot(coilset)
    view(2) %2D view
    c=colorbar; %colorbar
     ylabel(c, 'I_{VV} (kA)');
    xlabel(' R (m)')
    ylabel('Z (m)')
    %zlabel('I (A)')
    axis([0,1.03,-1.1,1.1]) %for the tfg EDDY Y FORCES
    title('Max Eddy current in each VVs filament')
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    grid off
    %print -depsc2 NOMBREPLOT.eps

%%%Plot for each instant
% for i=1:5%length(time_adaptive)
% figure;
% acumm=I_Passive(i,:)+acumm
%  scatter3(RR,ZZ,acumm,100,acumm,'filled')
%  view(2) %para verlo en 2D
% xlabel(' R (m)')
% ylabel('Z (m)')
% zlabel('I (A)')
% axis([0,1,-1.1,1.1]) %for the tfg
% title('sprintf(iter %d,i)')
% set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
% end

%% %%%%FORCES UPON A CROSS-SECTION OF THE VESSEL%%%%%%%%%%%
%For doing that, we firsly obtain the field, from the equilibrium:

%Extraction of the fields:
  [Fields]=fields(equil);

%The field at the vessel is:

Br_VV=Fields.interpn.Br(Rfil,Zfil);
Bz_VV=Fields.interpn.Bz(Rfil,Zfil);
Bphi_VV=Fields.interpn.Bphi(Rfil,Zfil);
        
B_vessel=[Br_VV' Bphi_VV' Bz_VV']; 
        %size(number of filaments*3), each row is the vector field on one filament

%Current
%The maximun currrent on each filament is I_Passive_fil 
%(size 1*number of filaments), whose unit vector is fi, so the vector
%would be in cylindrical

I_passive_fil_vec=I_Passive_fil'*[0 1 0]; %size number of filaments*3

%The force upon the cross-scetion of filament would then be

Force_fil=cross(I_passive_fil_vec,B_vessel); %size number of filaments*3, 
%so this is the force acting upon every filament

%The force upon all the filament would be 2piR*Force_fil_cross. R is stores
%in RR, which contains all the R values in a vector form with number of fil components. 
%Force_fil_cross is a vector of 3 components. It would be difficult to
%multiply them, but we do not need to, right now, because to obtain the
%pressure R cancles out, since the areas are 2piR*anchura (or altura). We
%assimilate the 3D filament as a 2D filament, so that it has no width in
%the R axis, s that its surface is 2piR*altura

%%%%Stresses
%This is th eplot of the pressures, much more better
%THIS ARE STRESSES, NOT PRESSURES!!!!

stress_R=(Force_fil(:,1))/(height_PF); 
%have the sign, to indicate wether goes to the inside or to the outside
stress_Z=(Force_fil(:,3))/(height_PF);

stress_Z_max=max(abs(stress_Z))
stress_R_max=max(abs(stress_R))

%An increase of 0.3Pa on R is found when adding B Earth!!


%plot
figure; 
 scale_factor=1%2*10^5; %graphic needs to be scaled
quiver(R_Fil_Array,Z_Fil_Array,stress_R/scale_factor,stress_Z/scale_factor,'color',[1 0 0],'AutoScale','on')
hold on;
plot(coilset)
plot(vessel);
xlabel('R (m)')
ylabel('Z (m)')
axis([-0.1,1.05,-1.3,1.3]) %for the tfg  EDDY Y FORCES
axis equal
set(gca,'XLim',[-0.5 1.5]);
title('Stresses on the vessel for phase 2')
set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%print -depsc2 NOMBREPLOT.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PASCHEN CURVE PLOT
 

C_1=[510 300]; %C_1 constant [ m^-1 Tor^-1] for H and He
C_2=[1.25e4 3.4e4];  %C_2 constant [V m^-1 Tor^-1]
Gas_type=["H_2","He"];
 p=linspace(1e-5,1e-3,100000);                           %[Tor]pressure of the prefill gas. size>1000 because if
                                                                            %not, It do not work properly
 
Emin= @(L,p,C1,C2) C2*p./log(C1*p*L); %Min E as Paschen state

        figure;
        %subplot(1,2,1)
        loglog(p,Emin(20,p,C_1(1),C_2(1)),'-','LineWidth', 0.95)
        hold on
        loglog(p,Emin(50,p,C_1(1),C_2(1)),'-','LineWidth', 0.95)
        loglog(p,Emin(70,p,C_1(1),C_2(1)),'-','LineWidth', 0.95)
        loglog(p,Emin(100,p,C_1(1),C_2(1)),'-','LineWidth', 0.95)
        %VEST
        loglog(p,3/(2*pi*0.36)*ones(1,length(p)),'b--','LineWidth', 0.95) 
        loglog([2e-5 2e-5],[10^-1 10^5],'b--','LineWidth', 0.95)
        loglog([3e-5 3e-5],[10^-1 10^5],'b-.','LineWidth', 0.95)
        %GlobusM
        loglog(p,4.5/(2*pi*0.36)*ones(1,length(p)),'k--','LineWidth', 0.95)
        loglog(p,8/(2*pi*0.36)*ones(1,length(p)),'k--','LineWidth', 0.95)
        loglog([3e-5 3e-5],[10^-1 10^5],'k--','LineWidth', 0.95)
        loglog([6e-5 6e-5],[10^-1 10^5],'k--','LineWidth', 0.95)
        %
        loglog(p,E(param_equil.rin)*ones(1,length(p)),'r-','LineWidth', 1.15)   
        loglog(p,E_Rgeo*ones(1,length(p)),'r--','LineWidth', 1.15) 
        xlabel('Prefill pressure (Torr)')
        ylabel('E_{min} (V/m)')
        legend('L=20m','L=50m','L=70m','L=100m','','VEST','','','GlobusM','','','SMART E(Rin)','SMART E (Rgeo)')
        title(sprintf('Paschen curve, Gas=%s',Gas_type(1))); %d for numbers
        set(gca, 'FontSize', 13); %<- Set properties TFG
        Filename='Paschen_complete';
        saveas(gcf, strcat(FigDir,Filename,FigExt));
    
 
 %%%Estimation of the time when the avalancha has ended (begin
 %%%burn-through)

%Imagine the values choosen are
%E_choosen=5 %V/m
%p_choosen=10^-5 %Torr (multiplo de Pa)

    %Plot3 varying E,p
    Ep=linspace(0.1,5e1,1000);  %[Tor]
    p_single=linspace(1e-5,1e-3,1000);  %[Tor]
    
    %Have to a meshgrid with p and E so that all the values of E are combined
    %with the p, and viceversa, which is was I want to see.
    [p,Ep]=meshgrid(p_single,Ep);
    
    %To plot this, will do as a I did with L inside the vessel. Since there
    %are problems with negative values, with substitute them with NaN. Have
    %to evaluate the function, which want only a single value of E and p
    %for not crashing or giving wrong things
    
    %Will do a loop to plot it for different L values
    clear tau
    tic
    L_plot=[L_aver/2 L_aver 2*L_aver 4*L_aver]             %[m]
    
    for i_L=1:length(L_plot)
    
        for i_p=1:length(p)
        
            for i_E=1:length(Ep)
            
                tau(i_p,i_E)=tau_ava(Ep(i_p,i_E),p(i_p,i_E),L_plot(i_L),C_1(1),C_2(1)); %[s] time 
            end
        end
      time_Ep_loop=toc  
        %Plot
        figure;
        contourf(p,Ep,log10(tau*1e3),10,'ShowText','on');
        %contourf(p,E,tau);
        shading('interp') %this is to make the transition between values continuous,
                                    %instedad of discontinuously between pixels
       %colormap(Gamma_II)
       %colormap('hot') 
       hold on;
        plot(p_single,E(param_equil.rin)*ones(length(p_single),1),'r-','LineWidth',1.15)
        plot(p_single,E_Rgeo*ones(length(p_single),1),'r--','LineWidth',1.15)
        c=colorbar; %colorbar
        ylabel(c, 'log10(tau_{ava} (ms))');
        view(2) %2D view
        xlabel('p (Tor)')
        ylabel('E (V/m)')
        set(gca,'Xscale','log') %para poner el eje x en log   
        set(gca,'Yscale','log') %para poner el eje x en log 
        grid on
        title(sprintf('Avalanche time for L=%3.1f m, ph2',L_plot(i_L)))
        %title(sprintf('time radiation wal[ms]) for L=%d m',L_plot(i_L)))
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG      
        Filename = 'tau_ava';
        Filename= sprintf('%s_L_%d',Filename,i_L);    
        saveas(gcf, strcat(FigDir,Filename,FigExt));

    end

  %Also, assuming the ionization fraction of 5%=0.05, the plasma current when
  %the line radiation peaks could be founded:
  
  E_ej=3%[V/m] example of field
  j_rad_wall=2*1.6*10^(-19)*1/(K_B*T_prefill)*E_ej*0.05 % [A/m^2] plasma current
                               %density when line rad peakd. Note that p cancels out
%asuming a surface of 2pi R, can find the current:
I_rad_wall= j_rad_wall*2*pi*0.2 %Pasma current [A], Assuming circular plasma of R=0.2m
    
 p_single=linspace(1e-5,1e-3,100000);  %[Tor]
 figure;
    loglog(p_single,Emin(20,p_single,C_1(1),C_2(1)),'LineWidth',0.95) 
    hold on
    %
    loglog(p_single,Emin(50,p_single,C_1(1),C_2(1)),'LineWidth',0.95) 
    loglog(p_single,Emin(70,p_single,C_1(1),C_2(1)),'LineWidth',0.95)
    loglog(p_single,Emin(100,p_single,C_1(1),C_2(1)),'LineWidth',0.95)
    loglog(p_single,E(param_equil.rin)*ones(length(p_single),1),'r-','LineWidth',1.15)
    loglog(p_single,E_Rgeo*ones(length(p_single),1),'r--','LineWidth',1.15)
    xlabel('Prefill pressure (Torr)')
    ylabel('E_{min} (V/m)')
    legend('L=20m','L=50m','L=70m','L=100m','E(Rin)','E(Rgeo)')
    title(sprintf('Paschen curve ph2, Gas=%s',Gas_type(1))); %d for numbers
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    axis([10^-5 10^-3 10^-1 50])
    Filename = 'Paschen_ava';
    %Filename= sprintf('%s_L_%3.1f',Filename,L_plot(i_L));    
    saveas(gcf, strcat(FigDir,Filename,FigExt));    

%%%%%%End loop paschen plots%%%%%%%%

    
        
  %Plot of the mean free path%%%%%
  %The mean free path is lambda=1/C_1*p. Can easily plot it for the
  %interesting gases, H, He and Ar.
  
  lambda= @(C1,p) 1./(C1.*p); 
  figure;
  plot(p_single,lambda(C_1(1),p_single))
  %hold on;
  %plot(p,lambda(C_1(2),p))
  xlabel('p (Tor)')
  ylabel('\lambda (m)')
  set(gca,'Xscale','log') %para poner el eje x en log   
  title('Mean free path \lambda versus p')
  legend(Gas_type(1))%,Gas_type(2))
  set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
  Filename = 'Lambda';
  %Filename= sprintf('%s_L_%3.1f',Filename,L_plot(i_L));    
 saveas(gcf, strcat(FigDir,Filename,FigExt));    
  %% SAVING THINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1)Current data
    %Now we save the waveform, the plasma current and the passive current

    %Create subdirectory for coil current related data
    icoilDir = strcat(ASCIIDir,'icoil_Data/'); mkdir(icoilDir);

    %I_PF_input=I_PF_output (readapted to time_adaptive), provided that
    %on RZIp, the show plot option must is true. If not, I_PF_output=NaN,
    %since it is not calculated, V_PF is calculated.
    %We save it in A*turn

    Filename = strcat(icoilDir,'Coil-waveform.txt');

    %1.1) I_PF
    a=I_PF_output(:,iSol).*turns(iSol);
    b=I_PF_output(:,iPF1).*turns(iPF1);
    c=I_PF_output(:,iPF2).*turns(iPF2);
    d=I_PF_output(:,iDiv1).*turns(iDiv1);
    e=I_PF_output(:,iDiv2).*turns(iDiv2);

    fileID=fopen(Filename,'w');
    fprintf(fileID,'%6s %6s %6s %6s %6s %6s\r\n', 'Time (s)','I_Sol (A turns)','I_PF1 (A turns)','I_PF2 (A turns)','I_Div1 (A turns)','I_Div2 (A turns)');
    fprintf(fileID,'%1.12f %0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[time_adaptive';a'; b'; c'; d'; e']);
    fclose(fileID) %Its mandatory to close the file, to avoid problems

    %1.2)Ip, plasma current, in A
    Filename = strcat(icoilDir,'Iplasma.txt');
    fileID=fopen(Filename,'w');
    fprintf(fileID,'%6s %6s\r\n', 'Time (s)','I_plasma (A)');
    fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive';Ip_output']);
    fclose(fileID)

    %1.3) I_passive, the eddy curretns, in A

    Filename = strcat(icoilDir,'Ipassive.txt');
    fileID=fopen(Filename,'w');
    fprintf(fileID,'%6s %6s\r\n', 'Time (s)','I_passive (A)');
    fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive';I_Passive_VV']);
    fclose(fileID)
    
%2) Eq related
    %Create subdirectory
%     EqDir = strcat(ASCIIDir,'Equil_Data/'); mkdir(EqDir);
%     
%     Filename = strcat(EqDir,'Geqsdk.txt');
%     geqdsk_write_BUXTON(config, equil, Filename)
%      clc %to remove all the warnings that this file makes  
%%THIS DO NOT WORK!!!!!





%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@FUNCTIONS!!!!!!!!!!!!@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@







%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@VV creation@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Constructs a set of rectilinear vessel wall filaments of variable thickness
%Set FilamentArea to 0.0 to autocompute filament area
%Set WallNorm to 'Static' to enable variable filament area - not recommended
function [vessel_filaments,R_Fil_Array,Z_Fil_Array]=...
    CreateRectilinearVessel(VesselDimensions,WallThickness,FilamentArea,WallNorm)
    
    %Initiate any required data or control arrays
    VesselFaces = ["Horizontal", "Vertical", "Horizontal", "Vertical"];
    R_Fil_Array = [];  Z_Fil_Array = [];
    dR_Fil_Array = []; dZ_Fil_Array = [];
	
    %Unpack wall corners into single parameters
    RMinCentre = VesselDimensions(1); RMaxCentre = VesselDimensions(2); %[m]
    ZMinCentre = VesselDimensions(3); ZMaxCentre = VesselDimensions(4); %[m]
 	
    %Determine constant filament area if not otherwise specified
    if FilamentArea == 0.0  
        FilamentArea = max(WallThickness)^2;    %[m^2]
    end
	
    %Define four vessel vertices going clockwise from inboard upper (V1)
    Vertice1 = [RMinCentre ZMaxCentre];	%Top Left
    Vertice2 = [RMaxCentre ZMaxCentre];	%Top Right
    Vertice3 = [RMaxCentre ZMinCentre];	%Bottom Right
    Vertice4 = [RMinCentre ZMinCentre];	%Bottom Left
    WallCorners = [[Vertice1]; [Vertice2]; [Vertice3]; [Vertice4]];
 	
    %Define four wall edges going clockwise from inboard upper vertex (V1)
    %Top and bottom walls 'own' their vertices and are sized to match radial wall extents
    WallEdge1 = [RMinCentre+WallThickness(4)/2, RMaxCentre-WallThickness(3)/2];  %Upper Wall
    WallEdge2 = [ZMaxCentre, ZMinCentre];  %Outboard Wall
    WallEdge3 = [RMaxCentre-WallThickness(3)/2, RMinCentre+WallThickness(4)/2];  %Lower Wall
    WallEdge4 = [ZMinCentre, ZMaxCentre];  %Inboard Wall
    WallEdges = [[WallEdge1]; [WallEdge2]; [WallEdge3]; [WallEdge4]];
    
    %Define normalisation factors for each wall to achieve filament area
    WallNormFactor1 = FilamentArea/(WallThickness(1)^2);
    WallNormFactor2 = FilamentArea/(WallThickness(2)^2);
    WallNormFactor3 = FilamentArea/(WallThickness(3)^2);
    WallNormFactor4 = FilamentArea/(WallThickness(4)^2);
    WallNormFactors = [WallNormFactor1, WallNormFactor2, WallNormFactor3, WallNormFactor4];

    %If requested, set variable filament area (not recommended)
    if WallNorm == "Static"
        WallNormFactors = [1, 1, 1, 1];
    end
    
    %For each wall, create a linspace of filament coordinates in (R,Z)
    %and a corresponding linspace of scaled width and height values.
    for i=1:length(VesselFaces)
        
        %If vessel face is horizontal, scale width of filament to maintain area
        if VesselFaces(i) == "Horizontal"
            Height = WallThickness(i);                               %Define thickness parallel to wall direction
            Width = WallThickness(i)*WallNormFactors(i);             %Define thickness perpendicular to wall direction
            NumFil = (RMaxCentre-RMinCentre+2*WallThickness(i))/WallThickness(i);
            NumFil = floor( NumFil/WallNormFactors(i) );             %Scale number of filaments to maintain total width

            R_fil = linspace(WallEdges(i,1),WallEdges(i,2),NumFil);  %Create evenly spaced array of wall R coordinates
            Z_fil = WallCorners(i,2)*ones(1,NumFil);                 %Create evenly spaced array of wall Z coordinates
            dR_Wall = linspace(Height,Height,NumFil);                %Create Axial wall thickness array of size NumFil
            dZ_Wall = linspace(Width,Width,NumFil);                  %Create Radial wall thickness array of size NumFil

            %Top and bottom walls 'own' their vertices - i.e. full width of filaments
            R_Fil_Array = [R_Fil_Array, R_fil];                      %Append filament R coordinates to array
            Z_Fil_Array = [Z_Fil_Array, Z_fil];                      %Append filament Z coordinates to array
            dR_Fil_Array = [dR_Fil_Array, dR_Wall];                  %Append filament radial widths to array
            dZ_Fil_Array = [dZ_Fil_Array, dZ_Wall];                  %Append filament axial heights to array

        %If vessel face is vertical, scale height of filament to maintain area
        elseif VesselFaces(i) == "Vertical"
            Height = WallThickness(i)*WallNormFactors(i);            %Define thickness parallel to wall direction
            Width = WallThickness(i);                                %Define thickness perpendicular to wall direction
            NumFil = (ZMaxCentre-ZMinCentre+2*WallThickness(i))/WallThickness(i);
            NumFil = floor( NumFil/WallNormFactors(i) );             %Scale number of filaments to maintain total height

            Z_fil = linspace(WallEdges(i,1),WallEdges(i,2),NumFil);  %Create evenly spaced array of wall Z coordinates
            R_fil = WallCorners(i,1)*ones(1,NumFil);                 %Create evenly spaced array of wall R coordinates
            dR_Wall = linspace(Width,Width,NumFil);                  %Create Axial wall thickness array of size NumFil
            dZ_Wall = linspace(Height,Height,NumFil);                %Create Radial wall thickness array of size NumFil

            %Radial walls don't own the vertices - remove first and last filaments
            R_Fil_Array = [R_Fil_Array, R_fil(2:end-1)];             %Append filament R coordinates to array (removing corners)
            Z_Fil_Array = [Z_Fil_Array, Z_fil(2:end-1)];             %Append filament Z coordinates to array (removing corners)
            dR_Fil_Array = [dR_Fil_Array, dR_Wall(2:end-1)];         %Append filament radial widths to array (removing corners)
            dZ_Fil_Array = [dZ_Fil_Array, dZ_Wall(2:end-1)];         %Append filament axial heights to array (removing corners)
        end 
    end

    %Construct vessel wall employing passive FIESTA filaments using position arrays
    for i=1:length(R_Fil_Array)
        vessel_filaments(i) = fiesta_filament(R_Fil_Array(i),Z_Fil_Array(i),dR_Fil_Array(i),dZ_Fil_Array(i),1,0,0);
    end
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Sol creation@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Function is a little dodgy at the moment - SolTurnsR > 1 requires CoilWaveforms Modifications

function [Sol_Circuit]=...
    CreateSMARTSolenoidCircuit(SolName,RMaxSol,RMinSol,ZMaxSol,ZMinSol,SolTurnsZ,SolTurnsR,SolTemp,SolResistivity,SolDensity)

    %SolTurnsZ :: Solenoid axial windings (true solenoid turns)
    %SolTurnsR :: Solenoid radial windings (proxy for radial resolution)

    %Initiate Z-coordinates for solenoid filaments and calculate width/height
    Z_filaments = linspace(ZMinSol,ZMaxSol,SolTurnsZ); clear('coil_filaments');
    
    %Calculate solenoid filament width, height and radial spacing
    SolFilWidth = (RMaxSol-RMinSol)/SolTurnsR;    %Filament width  [m]
    SolFilHeight = (2*ZMaxSol)/SolTurnsZ;             %Filament height [m]
    SolFilRadii = linspace(RMinSol+SolFilWidth/2, RMaxSol-SolFilWidth/2, SolTurnsR );   %[m]
    
    %Construct central solenoid filaments - solenoid is treated as 'vessel wall' with nonzero current
    for i=1:length(SolFilRadii)
        %Save each vertical filament stack as a seperate coil - moving from inner radius outwards
        for iFilament=1:SolTurnsZ
            Sol_filaments(iFilament) = fiesta_filament( SolFilRadii(i), Z_filaments(iFilament), SolFilWidth, SolFilHeight ); 
        end
        Sol_Coil(i) = fiesta_coil( 'Sol_Zcoil', Sol_filaments, 'Blue', SolResistivity, SolDensity );
    end
    
    %Compile vertical coil stacks into a single circuit
    Sol_Circuit = fiesta_circuit( SolName, ones(1,length(SolFilRadii)), Sol_Coil(:) );
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@PF creation@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [PF_circuit] = CreateSMARTPFCircuit(label,Rc,Zc,DR,DZ,nt, nZ, nR, simmetry, coil_temperature, resistivity, density)
% createVestPFCircuit(Rc,Zc,W,H,nt)
% Rc, Zc          coordinates of coil center
% DR, DZ          coil width and height
% nt              number of turns
% nv, nh          coils are distributed in an array of nv x nh coils

if rem(nt, nZ*nR) > 0
    warning('number of turns is not a multiple of nZ*nR');
    PF_circuit = false;
    return
end

c1 = create_coil( Rc,Zc,DR,DZ,nt, nZ, nR, coil_temperature, resistivity, density);
c1 = set(c1,'label','unique');
if simmetry
    c2 = create_coil( Rc,-Zc,DR,DZ,nt, nZ, nR, coil_temperature, resistivity, density);
    c1 = set(c1,'label','up');
    c2 = set(c2,'label','down');
    PF_circuit = fiesta_circuit(label,[1 1],[c1 c2]);
else
    PF_circuit = fiesta_circuit(label,[1],[c1]);    
end

end

%%%Asociateed function for creating PF circuit!!

function c=create_coil(Rc,Zc,DR,DZ,nt, nZ, nR, coil_temperature, resistivity, density)

turnsPerCoil = nt/(nZ*nR); 
[Zcoils,dz] = divideIntoIntervals( Zc-0.5*DZ, Zc+0.5*DZ, nZ);
[Rcoils,dr] = divideIntoIntervals( Rc-0.5*DR, Rc+0.5*DR, nR);
nFilament = nZ*nR;
for i=nZ:-1:1
    for j=nR:-1:1
        filament( nFilament) = fiesta_filament(Rcoils(j),Zcoils(i), dr, dz, turnsPerCoil, 0, 0);
        nFilament = nFilament - 1 ;
    end
end
%c = fiesta_coil('',filament);
c = fiesta_coil('',filament, 'Blue', resistivity, density ); %this may rsolve whats wrong?
end

%Associated function to create PF circuit

function [x,dx] = divideIntoIntervals(xi,xf,n)
% divide the interval [xi,xf] in n parts
% output: centers and width of the intervals
x = linspace(xi,xf,2*n+1);
dx = x(3)-x(1);
x = x(2:2:end);
end


%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Odefun for field line tracing LP@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Field line integrator function with Lp
     %this solves the field line eq, using Lp (poloidal length) 
         %as an independent variable
     %way more faster than with phi (5 min when Lmax=3000m, 15 inside points)
        
    function [results]=Field_LineIntegrator_Lp(Lp,rzLphiU,Br_interpn,Bz_interpn,Bphi_interpn)
    %rzphiLU=[r z phi L U]
    %Lp= poloidal length (have to write capital L so it not apperas as
    %internsity I). just tchange the variables in the phi equations
    
    %First, the field needs to be evaluated at the point (r,phi,z):
    
    Br_eval=Br_interpn(rzLphiU(1),rzLphiU(2));
    Bphi_eval=Bphi_interpn(rzLphiU(1),rzLphiU(2));
    Bz_eval=Bz_interpn(rzLphiU(1),rzLphiU(2)); 
    Bpol_eval=sqrt(Br_eval^2+Bz_eval^2);
    
    %With the field, the eq to solve is:
    
    dr_dLp=Br_eval/Bpol_eval;
    dphi_dLp=rzLphiU(1)*Bphi_eval/Bpol_eval;
    dz_dLp=Bz_eval/Bpol_eval;
    dlength_dLp=sqrt(1+(Bphi_eval/Bpol_eval)^2);
    dU_Vloop_dLp=1/(2*pi*rzLphiU(1))*dlength_dLp; %pseudo potential U/V_loop
    
    results=zeros(5,1); %column vector to group the results
    results(1)=dr_dLp;
    results(4)=dphi_dLp;
    results(2)=dz_dLp;
    results(3)=dlength_dLp;
    results(5)=dU_Vloop_dLp;
    
    end
    
    %% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Odefun for field line tracing PHI@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
     function [results]=Field_LineIntegrator(phi,rzLU,Br_interpn,Bz_interpn,Bphi_interpn)
    %rzL=[r z L U]
   
    %First, the field needs to be evaluated at the point (r,phi,z):
    
    Br_eval=Br_interpn(rzLU(1),rzLU(2));
    Bphi_eval=Bphi_interpn(rzLU(1),rzLU(2));
    Bz_eval=Bz_interpn(rzLU(1),rzLU(2));    
    Bpol_eval=sqrt(Br_eval^2+Bz_eval^2);
    
    %With the field, the eq to solve is:
    
    dr_dphi=rzLU(1)*Br_eval/Bphi_eval;
    dz_dphi=rzLU(1)*Bz_eval/Bphi_eval;
    length=rzLU(1)*sqrt(Bphi_eval^2+Bpol_eval^2)/Bphi_eval;
    U_Vloop=1/(2*pi*rzLU(1)); %pseudo potential U/V_loop
    
    results=zeros(4,1); %column vector to group the results
    results(1)=dr_dphi;
    results(2)=dz_dphi;
    results(3)=length;
    results(4)=U_Vloop;
    
  end	
    
%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Event function for the ode@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   %Only one function can be introduced into ode, so this event has to have
   %all the conditions
   
   function [rz_value isterminal direction]=Colission_wall(phi,rzL,...
       VesselRMaxPoint,VesselRMinPoint,VesselZMaxPoint,VesselZMinPoint,...
       Br_interp,Bz_interp,Bphi_interp,L_max) 
   
   %when rz_value is zero, stop the integration. To implement the 4
   %possible colission, we could do like the product of each, since when one of
   %them is zero, the product will be zero, and also to define row vectors
   %for isterminal, direction, and rz_value. Will do this second option
   
   isterminal=[1 1 1 1 1];                                                  %to stop the integration
   direction= [0 0 0 0 0];                                                  %to not worry about the sloop 
   
   up_colission=rzL(2)-VesselZMaxPoint; 
   down_colission=rzL(2)-VesselZMinPoint;
   out_colission=rzL(1)-VesselRMaxPoint; 
   in_colission=rzL(1)-VesselRMinPoint;
   
   %For the max condition of L, the field needs to be introduced:
   
    Br_eval=Br_interp(rzL(1),rzL(2));
    Bphi_eval=Bphi_interp(rzL(1),rzL(2));
    Bz_eval=Bz_interp(rzL(1),rzL(2));
   
   L_lim=rzL(3)-L_max;                                                  %[m] Maximum L value=20
   
   rz_value=[up_colission down_colission out_colission in_colission L_lim];
   
    %Have checked that if I dont use the minR condition, it will impige 
    %in the upper wall, which was was happens at the beggining, when
    %I dont have the inner wall condition
   end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Sensors creation@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Initiate virtual sensors within null-field region
function SensorBTheta=InitiateBSensors(EquilParams,length_R,R_centre,Z_centre,length_Z)

    %Determine number of sensors (constant for now) in each direction,  
    NumSensors = 10;         %10*10=100
    
    %Geo variables
    RGeo = EquilParams.r0_geom; ZGeo=EquilParams.z0_geom;
    
    switch nargin
        
        case 2 %Square region centered at (RGeo,ZGeo=0)
            length_Z=length_R;

        BP_virt_R = linspace(RGeo-length_R,RGeo+length_R,NumSensors); 
        BP_virt_Z = linspace(ZGeo-length_Z,ZGeo+length_Z,NumSensors);

        
        case 3 %square region centered elsewhere (Z=0)
            length_Z=length_R;
    
        BP_virt_R = linspace(R_centre-length_R,R_centre+length_R,NumSensors); 
        BP_virt_Z = linspace(ZGeo-length_Z,ZGeo+length_Z,NumSensors);
        
        case 4 %square region centered elsewhere (both R and Z)
            length_Z=length_R;
            
        BP_virt_R = linspace(R_centre-length_R,R_centre+length_R,NumSensors);
        BP_virt_Z = linspace(Z_centre-length_Z,Z_centre+length_Z,NumSensors);        

        
        case 5 %reactangular region centered elsewhere (both R and Z)
        
        BP_virt_R = linspace(R_centre-length_R,R_centre+length_R,NumSensors);
        BP_virt_Z = linspace(Z_centre-length_Z,Z_centre+length_Z,NumSensors); 
        
    end
    %Create null field region grid
    [BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
    BP_virt_R = BP_virt_R(:)';
    BP_virt_Z = BP_virt_Z(:)';

    %Create sensors over the null field region
    BP_virt_theta = zeros(1,length(BP_virt_R));
    nSensors = length(BP_virt_theta);

    %Create the array of sensors, sepparate radial (Br) and vertical (Bz)
    BP_virt_names = {};
    for iSensor=1:nSensors
        BP_virt_names{iSensor} = ['Radial Bp Virtual Sensor #' num2str(iSensor) ];
    end

      %Now Bz sensors
    for iSensor=nSensors+1:2*nSensors
        BP_virt_names{iSensor} = ['Vertical Bp Virtual Sensor #' num2str(iSensor) ];
    end
    
    %R and Z of Dim[1*200] and cyclical (i.e. BP_virt_R[201] == BP_virt_R [1]) 
    BP_virt_R = [BP_virt_R  BP_virt_R];
    BP_virt_Z = [BP_virt_Z  BP_virt_Z];
                        %Both size 1*200. It is replicated, so element 101=element 1 
    BP_virt_theta = [BP_virt_theta  BP_virt_theta+pi/2];    %size 1*200. The first 100 have 0, and the second has pi/2 

%Taken from a sensor Btheta function:
    %     theta=0.00 --> sensor is pointing in the R direction
    %     theta=pi/2 --> sensor is pointing in the Z direction
%This confirm you are creating both sensors to measure Br and Bz, by
%setting the angle os the sensors. Could also use sensorBZ and sensorBr
%functions, of course.
    
    %Initiate sensors with sensor names over defined region
    SensorBTheta = fiesta_sensor_btheta('BSensors',BP_virt_R,BP_virt_Z,BP_virt_theta,BP_virt_names);
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Extraction of the fields@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%1)Field of things
function [Field_Earth Field_NoEarth]=fields(equil)
%function [Field_grid, Field_VV, Field_interpol, Field_sensor]=fields(equil)
    %beware of the sensors, where they are not defined (yet), its field must not
    %be called!
    %Will add the Earths field! Will be useful when the fields are low, for
    %example at breakdown.
    
    global R_in Z_in R_sensor Z_sensor %to obtain global variables

    %Fiesta fields:
    Br=get(equil,'Br'); 
    Bz=get(equil,'Bz'); 
    Bphi=get(equil,'Bphi_vac'); %not alwais
    
    %Fields data (data in all the grid)
    Br_data = get(get(equil,'Br'),'data');     %Note this is 200*251, GridR*GridZ dimension
    Bz_data = get(get(equil,'Bz'),'data');
    Bphi_data = get(get(equil,'Bphi_vac'),'data');
            %Was Bphi, not always is vac, although the plasma field is
            %negligible in comparison to TF coils field.

    %%2D reshape of the data, to interpolate things
    
    rGrid = get(get(get(equil,'Br'),'grid'),'r');
    zGrid = get(get(get(equil,'Br'),'grid'),'z');
    
    Br_data = reshape( Br_data, length(zGrid), length(rGrid)); %251*200
    Bz_data = reshape( Bz_data, length(zGrid), length(rGrid)); %251*200
    Bphi_data = reshape( Bphi_data, length(zGrid), length(rGrid)); %251*200
    Bpol_data=sqrt(Br_data.^2+Bz_data.^2);
   
    
    %Interpolation vectors
    Br_interp = @(r,z) interpn(zGrid,rGrid,Br_data,z,r,'mikama');
    Bz_interp = @(r,z) interpn(zGrid,rGrid,Bz_data,z,r,'mikama');
    Bphi_interp = @(r,z) interpn(zGrid,rGrid,Bphi_data,z,r,'mikama');
    
    %Finally, the fields inside the vessel are
    Br_VV=Br_interp(R_in,Z_in);
    Bz_VV=Bz_interp(R_in,Z_in);
    Bphi_VV=Bphi_interp(R_in,Z_in);
    Bpol_VV=sqrt(Br_VV.^2+Bz_VV.^2);
    
    %And inside the sensor region:
    Br_sens=Br_interp(R_sensor,Z_sensor);
    Bz_sens=Bz_interp(R_sensor,Z_sensor);
    Bphi_sens=Bphi_interp(R_sensor,Z_sensor);
    Bpol_sens=sqrt(Br_sens.^2+Bz_sens.^2);
    
    %To store them, will create several structures, one for the data, other
    %for the VV data and other for the sensor data
    
    Field_grid.Br=Br_data;
    Field_grid.Bz=Bz_data;
    Field_grid.Bphi=Bphi_data;
    
    Field_VV.Br=Br_VV;
    Field_VV.Bz=Bz_VV;
    Field_VV.Bphi=Bphi_VV;
    Field_VV.Bpol=Bpol_VV;    

    Field_sensor.Br=Br_sens;
    Field_sensor.Bz=Bz_sens;
    Field_sensor.Bphi=Bphi_sens;
    Field_sensor.Bpol=Bpol_sens;
    
    Field_interpol.Br=Br_interp;
    Field_interpol.Bz=Bz_interp;
    Field_interpol.Bphi=Bphi_interp;
    
    %Final structure englobating all
    Field_NoEarth.grid=Field_grid;
    Field_NoEarth.VV= Field_VV;
    Field_NoEarth.sensor= Field_sensor;    
    Field_NoEarth.interpn=Field_interpol;
    
        %%%%%Earths field%%%%%%%%%%%%%
    
    %In seville, with coordinates 37º23'N 5°59′00″W (Wikipedia, spanish)=
    %37+23/60ºN 5+59/60º W=37.3833ºN 5.9833ºW, the components are
    
    X_Earth=27316.6e-9;                 %[T] N-S component, >0 for N, <0 for S
    Y_Earth=-423.4e-9;                      %[T] E-W component, <0 for W, >0 for E
    Z_Earth=33833.4e-9;                 %[T] vertical component, <0 for U, >0 for D
    
    %To create the vectors, I have problems for the r and phi directions, since its
    %magnitude vary when moving the toroidal angle, because Fiesta considers
    %axysymmetry, but Earths field is not axysymmetric. What I will do as an
    %approximation is take the average value of the field over all the phi angles. 
    %This value, the mean, will be used for the r and phi components. 
    %THe z components is norproblematic
    
        BrEarth=0; %initialization
        BphiEarth=0; %initialization
        ang_Earth=linspace(0,2*pi,100);
        
        for i=1:length(ang_Earth)
            %Addition at each step (debug)
            add_Br(i)=X_Earth*cos(ang_Earth(i))+(-Y_Earth)*sin(ang_Earth(i));       
            add_Bphi(i)=-X_Earth*sin(ang_Earth(i))+(-Y_Earth)*cos(ang_Earth(i));
                                                %-Y because Seville is in the West (W)
            BrEarth=BrEarth+add_Br(i);
            BphiEarth=BphiEarth+add_Bphi(i);
            
        end
         
        %The r and phi components are the mean:       
        BrEarth=BrEarth/length(ang_Earth); %mean
        BphiEarth=BphiEarth/length(ang_Earth); %mean      
        
        %The vertical component is
        BzEarth=-Z_Earth;
        BpolEarth=sqrt(BrEarth^2+ BzEarth^2);
    
        %Now lets create the a grid with this constant field values;
       
        BzEarth=BzEarth*ones(length(zGrid), length(rGrid));        
        BrEarth= BrEarth*ones(length(zGrid), length(rGrid));
        BphiEarth=BphiEarth*ones(length(zGrid), length(rGrid));    
        BpolEarth=sqrt(BrEarth.^2+ BzEarth.^2);
    
        
    %%%%%%%%END EARTHS FIELD%%%%%%%%
    %Now will create the same as above but with Earths field
    
    %Interpolation vectors 
    Br_interp = @(r,z) interpn(zGrid,rGrid,Br_data+BrEarth,z,r,'mikama');
    Bz_interp = @(r,z) interpn(zGrid,rGrid,Bz_data+BzEarth,z,r,'mikama');
    Bphi_interp = @(r,z) interpn(zGrid,rGrid,Bphi_data+BphiEarth,z,r,'mikama');
    
    %Finally, the fields inside the vessel are
    Br_VV=Br_interp(R_in,Z_in);
    Bz_VV=Bz_interp(R_in,Z_in);
    Bphi_VV=Bphi_interp(R_in,Z_in);
    Bpol_VV=sqrt(Br_VV.^2+Bz_VV.^2);
    
    %And inside the sensor region:
    Br_sens=Br_interp(R_sensor,Z_sensor);
    Bz_sens=Bz_interp(R_sensor,Z_sensor);
    Bphi_sens=Bphi_interp(R_sensor,Z_sensor);
    Bpol_sens=sqrt(Br_sens.^2+Bz_sens.^2);
    
    %To store them, will create several structures, one for the data, other
    %for the VV data and other for the sensor data
    
    Field_grid.Br=Br_data;
    Field_grid.Bz=Bz_data;
    Field_grid.Bphi=Bphi_data;
    
    Field_VV.Br=Br_VV;
    Field_VV.Bz=Bz_VV;
    Field_VV.Bphi=Bphi_VV;
    Field_VV.Bpol=Bpol_VV;    

    Field_sensor.Br=Br_sens;
    Field_sensor.Bz=Bz_sens;
    Field_sensor.Bphi=Bphi_sens;
    Field_sensor.Bpol=Bpol_sens;
    
    Field_interpol.Br=Br_interp;
    Field_interpol.Bz=Bz_interp;
    Field_interpol.Bphi=Bphi_interp;
    
        %Final structure englobating all
    Field_Earth.grid=Field_grid;
    Field_Earth.VV= Field_VV;
    Field_Earth.sensor= Field_sensor;    
    Field_Earth.interpn=Field_interpol;
    
    %Yes, everything is all right, you are not rewritting things
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Earth's field@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    %Not in use, but its helpful to see the Earths field alone!

function [Field_Grid Field_Value] =EarthField(RGrid,ZGrid)

    %In seville, with coordinates 37º23'N 5°59′00″W (Wikipedia, spanish)=
    %37+23/60ºN 5+59/60º W=37.3833ºN 5.9833ºW, the components are
    
    X_Earth=27316.6e-9;                 %[T] N-S component, >0 for N, <0 for S
    Y_Earth=-423.4e-9;                      %[T] E-W component, <0 for W, >0 for E
    Z_Earth=33833.4e-9;                 %[T] vertical component, <0 for U, >0 for D
    
    %To create the vectors, I have problems for the r and phi directions, since its
    %magnitude vary when moving the toroidal angle, because Fiesta considers
    %axysymmetry, but Earths field is not axysymmetric. What I will do as an
    %approximation is take the average value of the field over all the phi angles. 
    %This value, the mean, will be used for the r and phi components. 
    %THe z components is norproblematic
    
        BrEarth=0; %initialization
        BphiEarth=0; %initialization
        ang_Earth=linspace(0,2*pi,100);
        
        for i=1:length(ang_Earth)
            %Addition at each step (debug)
            add_Br(i)=X_Earth*cos(ang_Earth(i))+(-Y_Earth)*sin(ang_Earth(i));       
            add_Bphi(i)=-X_Earth*sin(ang_Earth(i))+(-Y_Earth)*cos(ang_Earth(i));
                                                %-Y because Seville is in the West (W)
            BrEarth=BrEarth+add_Br(i);
            BphiEarth=BphiEarth+add_Bphi(i);
            
        end
         
        %The r and phi components are the mean:       
        BrEarth=BrEarth/length(ang_Earth); %mean
        BphiEarth=BphiEarth/length(ang_Earth); %mean      
        
        %The vertical component is
        BzEarth=-Z_Earth;
        BpolEarth=sqrt(BrEarth^2+ BzEarth^2);
    
        %Now lets create the a grid with this constant field values;
       
        BzEarthG=BzEarth*ones(length(ZGrid), length(RGrid));        
        BrEarthG= BrEarth*ones(length(ZGrid), length(RGrid));
        BphiEarthG=BphiEarth*ones(length(ZGrid), length(RGrid));    
        BpolEarthG=sqrt(BrEarthG.^2+ BzEarthG.^2);

    %Store of the values
    
    Field_Grid.Br=BrEarthG;
    Field_Grid.Bz=BzEarthG;
    Field_Grid.Bphi=BphiEarthG;
    Field_Grid.Bpol=BpolEarthG;

    Field_Value.Br=BrEarth;                     %single values of the fields
    Field_Value.Bz=BzEarth;
    Field_Value.Bphi=BphiEarth;
    Field_Value.Bpol=BpolEarth;
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Fit Sol ramp@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Specify or extrapolate coil current for a given current ramp timescale
function CoilRampCurrent=FitSolenoidRamp(CoilRampCurrents,TimeVertices)

    %Extract relevent coil currents
    CoilStartCurrent = CoilRampCurrents{1};
    CoilRampCurrent = CoilRampCurrents{2};
    CoilEndCurrent = CoilRampCurrents{3};

    %Extrapolate a linear ramp-down:
    if strcmp(CoilRampCurrent, 'Linear');
        %Maintain a linear solenoid ramp-down from time(4), through time(5) to time (6)
        %Apply a linear fit to the solenoid ramp-down profile between PrePulse to Equil
        [coef] = polyfit([TimeVertices(3), TimeVertices(5)], [CoilStartCurrent, CoilEndCurrent], 1);
        %Extrapolate solenoid current when PF and Div coils reach equilibrium values, at time(4)
        CoilRampCurrent = (coef(1)*TimeVertices(4)) + coef(2);
   
    %Employ user specified value if requested
    elseif isfloat(CoilRampCurrent);
        CoilRampCurrent = double(CoilRampCurrent);
    end
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@avalanche time@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [tau_radwall]=tau_ava(E,p,L,C_1,C_2)
    %This is a function to calculate the time the avalanche ends
    
    global T_prefill K_B
    
    %If it is negative, I will set NaN, so that it do not plot it
    alpha=C_1*p.*exp(-C_2*p./E);    %First townsend coefficient
    
    f_i=0.05; % [%] of ionization fraction when the avalanche ends [Iter2019]
    
    if alpha-1/L<=0 %no avalanche
        tau_radwall=NaN;
    else%avalanche
        tau_radwall= log(2*f_i*p/(K_B*T_prefill))./(43*E./p.*(alpha-1/L)); %41 in the cotient says LLoyd
    end

    
end

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Function for the colormap (given by Jose Rueda)@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   function [map] = Gamma_II(t)
%--------------------------------------------------------------------------
% This function creates the colormap that coincides with the
% Gamma_II_colormap of IDL. 
% INPUTS:
%    - t: Number of points to define the color scale (default = 256)
%--------------------------------------------------------------------------

% Initialise settings
if nargin < 1 
    t=256;
end

% Colors of the scale
T = [0,   0,   0            % dark
     0, 0,  255             % blue
     255, 0, 0              % red
     255, 255, 0            % yellow
     255, 255, 255]/255;   % white
 
 % Setting color intervals length
 x = [0
     70
     130
     200
     255]/255;
 
 %Interpolation between colors
 map = interp1(x,T,linspace(0,1,t));

   end