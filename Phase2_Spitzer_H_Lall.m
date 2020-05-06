%% SEVILLE SPHERICAL TOKAMAK STIMATION
%V2C2 (the winner). PHASE 1, REDUCED SIZE
%
%Daniel López Aires // danlopair@alum.us.es
%Updated to Juanjo Toledo, but modified, and optimized
%Includes, at the end, eddy plots and stresses stimations

clear 
clc
close all

% %print -depsc2 NOMBREPLOT.eps %this is for saving plots in .eps, without
% saving margins
tic

%#################################################
%######################TFM#########################
%#################################################

%%%PARAMETER TO BE VARIED FOR THE BREAKDOWN!!!!!!!!!!!!!!!!!!!!!
    %1)GAS TYPE

        %i)H    
        %Z_eff=Z_nucleus-Debye shielding. It can be calculated in
        %http://calistry.org/calculate/slaterRuleCalculator , giving:
             Gas_type='H'
             Z_eff=1                                                                    %H
             C_1=510                              %[m-1 Torr-1] This is the constant A, I have changed its name
             C_2=1.25e4                     %[V-1 m-1 Torr-1] This is the constant B, I have changed its name

        %ii)He  
%              Gas_type='H'
%              Z_eff=2 %He
%              C_1=300                                %[m-1 Torr-1] This is the constant A, I have changed its name
%              C_2=3.4e4                        %[V-1 m-1 Torr-1] This is the constant B, I have changed its name
            
        %iii)Ar  
            %Gas_type='Ar'
             %Z_eff=11.85 %Ar
             %C_1=?                                  %[m-1 Torr-1] This is the constant A, I have changed its name
             %C_2=?                                 %[V-1 m-1 Torr-1] This is the constant B, I have changed its name
            
    %2) Greenwald fraction
    
                Gr_fraction=0.7%0.15;                              %This is to scale the Gr_limit==> <=1
                                        %0.7 to match Elis calc
    %3) a_eff                                                        %[m] This defines the size of the field null region.
   
             %a_eff=0.05;                                   % little null region WHAT JJ USED, ST25D BASED
             a_eff=0.15;                                        % large null region

%%%%%%%%END PARAMETERS TO BE VARIED ON BREAKDWON!!!!!!!!!!!!!

%#################################################
%#####################TFM############
%#################################################


%%%%%%%%%%%%%%%%%  DEFINE DATA OUTPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%
%Taken from the God Scott

%Define figure extension
FigExt = '.png'; 		%'.png','.eps','.pdf'

%Define project and series names
ProjectName = 'SMART-P2-9L';			%Define global project name
%SeriesName = 'VaryTauP';		%Define parameter scan series name

%Create global output folders for saved data and figures
%ASCIIDir = 'RawData/'; mkdir(ASCIIDir);
%FigDir = 'Figures/'; mkdir(FigDir);			%NOT CURRENTLY USED

%Create simulation name based upon relevant run parameters
SimName = 'DefaultSimName';
disp([ 'SimName: ' SimName ]);
disp([ ' ' ]);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@CREATION OF THE TOKAMAK@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%%%SEVILLE TOKAMAK(ST)  outer profile%%%%%

%The heigh of vessel, PF,DiV, and Sol (Sol does not matter for the eq calculation
%of course, since it is turned off, but if will be important in RZIp) 
%All will be modified with
mod_Z_pos=-0.2;                                                   %[m], modification for Z>0
mod_Z_neg=+0.2;                                                 %[m], modification for Z<0


%This are values of the internal vessel. For the outer, the width and heigh
%must be considered

VesselRMinPoint=0.15;                                               %[m] R min position
VesselRMaxPoint=0.8;                                                %[m] R max position
VesselZMinPoint=-1+mod_Z_neg;                               %[m] Z min position 
VesselZMaxPoint=1+mod_Z_pos;                                %[m] Z max position 

%R,Z of the filament: for this, whe widht and height have to be considered,
%since fiesta needs major R,Z to cosntruct the filament, and we hace the
%inner limits of the vessel:.

ww_R=2.01e-2%1.5e-2;                                                 %[m] filament width
ww_Z=ww_R;         %ww_R originally                            %[m] filament heigh
    

ZMinCenter=VesselZMinPoint-ww_Z/2;
ZMaxCenter=VesselZMaxPoint+ww_Z/2;
RMinCenter=VesselRMinPoint-ww_R/2;
RMaxCenter=VesselRMaxPoint+ww_R/2;

%They do follow an order!
point1=[RMaxCenter ZMaxCenter ];
point2=[RMaxCenter ZMinCenter];
point3=[RMinCenter ZMinCenter];
point4=[RMinCenter ZMaxCenter];

%%%%%%  Make PF coils  %%%%%%%%%%%%%%
%
%
coil_temperature = 293;                                 %[K] Temperature of the coils. value set by ST25D

% Coil indicies and turns, 
iSol = 1;       %inductor coil
iDiv1 = 4;
iDiv2 = 5;
iPF1 = 2;
iPF2 = 3;

% Turns in R and Z of the coils
nZDiv1=6;
nRDiv1=4;
nZDiv2=6;
nRDiv2=4;
nZPF1=6;
nRPF1=4;
nZPF2=6;
nRPF2=4;

nDiv1=nZDiv1*nRDiv1;
nDiv2=nZDiv2*nRDiv2; 
nPF1=nZPF1*nRPF1; 
nPF2=nZPF2*nRPF2;  

%Dimensions of a turn
width_PF=0.042;                                             % [m] Width of a turn
height_PF=0.035;                                        %  [m] Height of a turn 

%Position of the coils (m)
R_PF1=0.9+ww_R/2; 
Z_PF1=0.5+mod_Z_pos+ww_Z/2; 
R_PF2=0.9+ww_R/2;
Z_PF2=0.8+mod_Z_pos+ww_Z/2;
R_Div1=0.15+ww_R/2; %0.25 originally (me) Updated to Scott
Z_Div1=1.05+mod_Z_pos+ww_Z/2;
R_Div2=0.45+ww_R/2; %0.55 originally (me) Updated to Scott
Z_Div2=1.05+mod_Z_pos+ww_Z/2;

% Make Solenoid
nSol=210 %800;                                      % number of turns of the solenoid (Agredano suggested)

%RSol=0.09-ww_R/2;                                  %[m] R position of the solenoid (Inner Solenoid)
 RSol=0.13-ww_R/2;                                    % [m] R position of the solenoid (Outer Solenoid) 
 
ZMax_Sol=RMaxCenter+ww_Z/2;                                   % [m] Max Z position
ZMin_Sol=-ZMax_Sol;                                                  % [m] Min Z position

turns=[];
turns(iSol) =nSol; %100 in my tfg
turns(iDiv1) = nDiv1; %8 in my tfg
turns(iDiv2) = nDiv2; %8 in my tfg
turns(iPF1) = nPF1; %24 in my tfg
turns(iPF2) = nPF2; %24 in my tfg

nPF = length(turns);                                %The number of total coils, counting poloidal and inductor

%Creation of the coils%%%%%%%
% the function createVESTPFCircuit (made by Carlos Soria) create two PF
% coils. One in (R, Z) and another in (R, -Z)

resistivity = copper_resistivity_at_temperature( coil_temperature );
density = 1; %default value in FIESTA, ST 25D

PF1  = createVestPFCircuit( 'PF1',R_PF1,Z_PF1,width_PF,height_PF,turns(iPF1),nZPF1,nRPF1,true, coil_temperature, resistivity, density);
PF2  = createVestPFCircuit( 'PF2',R_PF2,Z_PF2,width_PF,height_PF,turns(iPF2),nZPF2,nRPF2,true, coil_temperature, resistivity, density);
Div1 = createVestPFCircuit('Div1', R_Div1, Z_Div1, width_PF,height_PF, turns(iDiv1), nZDiv1,  nRDiv1, true, coil_temperature, resistivity, density); 
Div2 = createVestPFCircuit('Div2', R_Div2, Z_Div2, width_PF,height_PF, turns(iDiv2), nZDiv2,  nRDiv2, true, coil_temperature, resistivity, density);

%%%%%%Inductor coil (Sol)%%%%%%%%%%%%%%%

nfil_ind_coil=turns(iSol); %number of filaments of the inductor coil=number of turns (yes, check it if you want)
clear('coil_filaments');
Z_filament = linspace(ZMin_Sol,ZMax_Sol,nfil_ind_coil); 
for iFilament=1:nfil_ind_coil
    coil_filaments(iFilament) = fiesta_filament( RSol,Z_filament(iFilament), sqrt(70e-6),sqrt(70e-6) ); 
    %values by default, and R value defined so that the solenoid fit its
    %region
end

coil_1  = fiesta_coil( 'psh_coil', coil_filaments, 'Blue', resistivity, density );
Sol_circuit = fiesta_circuit( 'Sol', [1], [coil_1] );

%%%%Creation of the vessel object%%%%%%%%%%%%%%%

 %Lets create the lines manually
 
 n_fil_Z=round((ZMaxCenter-ZMinCenter+2*ww_Z)/ww_Z) %8;         %Number of filaments in the Z direction, from Zmin to Zmax
 n_fil_R=round((RMaxCenter-RMinCenter+2*ww_R)/ww_R) %8;        %Number of filaments in the R direction, from Rmin to Rmax

    %From point 1 to point 2, R=Rmax
    R_lin1_2=RMaxCenter*ones(1,n_fil_Z);
    Z_lin1_2=linspace(ZMaxCenter,ZMinCenter,n_fil_Z);
 
    %From point 2 to point 3, Z=Zmin
    R_lin2_3=linspace(RMaxCenter,RMinCenter,n_fil_R);
    Z_lin2_3=ZMinCenter*ones(1,n_fil_R);
 
   %From point 3 to point 4, R=Rmin
     R_lin3_4=RMinCenter*ones(1,n_fil_Z);
     Z_lin3_4=linspace(ZMinCenter,ZMaxCenter,n_fil_Z);

    %From point 4 to point 1, Z=Zmax
    R_lin4_1=linspace(RMinCenter,RMaxCenter,n_fil_R);
    Z_lin4_1=ZMaxCenter*ones(1,n_fil_R);
 
%       figure;
%       plot(R_lin1_2,Z_lin1_2,'ro');
%        xlabel('R (m)')
%        ylabel('Z (m)')
%       hold on
%        plot(R_lin2_3,Z_lin2_3,'b*');
%       plot(R_lin3_4,Z_lin3_4,'go');
%       plot(R_lin4_1,Z_lin4_1,'k*');
      
%Accumualtions of R,Z points:
xaccum=[R_lin1_2 R_lin2_3 R_lin3_4 R_lin4_1]';
yaccum=[Z_lin1_2 Z_lin2_3 Z_lin3_4 Z_lin4_1]';

    %accum_rep=[xaccum; yaccum];    %Check

%Remove duplicates
dup = (abs(diff(xaccum))+abs(diff(yaccum))) > 0;
xaccum = xaccum(dup);
yaccum = yaccum(dup);

% figure;
% plot(xaccum,yaccum,'r.')
% xlabel('R (m)')
% ylabel('Z (m)')

    %accum=[xaccum; yaccum];    %Check

%Creation of the vessel filaments
    %The dimensions of the filaments are a key factor to the eddy currents;
    %if the width (R) increases, the eddy increases, since you induced a
    %current density, so the higher the surface, the higher the total
    %current

for i=length(xaccum):-1:1
    vessel_filament(i) = fiesta_filament(xaccum(i),yaccum(i),ww_R,ww_Z,1,0,0); 
        %The fiesta_filament inputs are R,Z,2*r,2*z,1,0,0 R mayor radius, Z heigh, 
        %r minor radius, z minor z (2r width in R axis, 2z in Z axis of the filament)
end

%Creation of the vessel passives

passive = fiesta_passive('STVesselPas',vessel_filament,'g');
    
    %To simplify the calculus of the vessel
    %nmodes=300; %nmodes for the vessel, to simplify calculus. Originally, n=number of filaments
        %3 alters too much the eddys with respect to not choosing nmodes
        %5 too
        %10 reports an error on the second eq calc
        %100 altered
        %200 altered too
        %300 error, out of memory ==> :)))
    %passive=set(passive,'nmodes',nmodes)
    %in this you could define the vessel resistivity

   
%Finally, creation of the vessel object

vessel = fiesta_vessel( 'STVessel',passive);

%Plot of the cross section
    figure;
    plot(PF1);
    hold on;
    plot(PF2);
    plot(Div1);
    plot(Div2);
    plot(Sol_circuit);
    plot(vessel);
    axis equal;
    xlabel('R (m)')
    ylabel('Z (m)')
    title('Cross-section for the reduced size')
    %%%OPtionf for tfg
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    axis([0,1.1,-1.1,1.1]) 
    %print -depsc2 NOMBREPLOT.eps
    Filename = '_CrossSection';
    saveas(gcf, strcat(ProjectName,Filename,FigExt));

    
%%%%Borders

%%%%%%%%%%%%%%%%%%%%COILSET%%%%%%%%%%%%%%

coilset = fiesta_coilset('STcoilset',[Sol_circuit,PF1,PF2,Div1,Div2],false,xaccum',yaccum');

% @@@@@@@END CREATION OF THE TOKAMAK@@@@@@@@@@@@@@@


%%%%%%CURRENT PROFILE FOR THE COILSET

%For simplicity, and following Scott code, will declare the input profile
%at the beginnin, since this will be clearer!

%Notes:
    %-Provided that Ip>0, I_PF1,2 have to be<0, so it repel the plasma, and it do not 
        %touches the outer (R max) wall
    
%%%PF coil currents (At Equilibrium, time(5,6))                                                   [A]
I_PF1_Equil=-1.5e3;					%First guess, efit will found it. 
I_PF2_Equil=-.5e3;					%First guess, efit will found it. 
                                               %These first guesses have to
                                               %be accurate!!!! It affects
                                               %the results!!!!!!!!!
%I_Div1_Equil=NaN;					Div1 in series with Sol!
I_Div2_Equil=+3.3e3;			%4e3 		

%%%Sol current
 I_Sol_max =3.5e3%2.5e3                %[A] Max solenoid current, to achieve the desire Ip in RZIp.
ISol_equil=-I_Sol_max;%.05e3;                            % [A] the current of the Sol in the target eq calc
%If>0, Sol will attract the plasma, no lower PF and DIv are
%needed.


%%%Time intervals
nTime = 7;      	% Coil Waveform Timesteps	[-]
discharge_time=100;                                                          %[ms], time duration of the flat-top 
T_ramp_Sol=25;       %Seta as 25ms to match JJ               %[ms], min time to increase or decrease the Sol

t_add=2*T_ramp_Sol; %2 with short breakdown               %time for the aditional point>T_ramp_Sol

time = [-4*T_ramp_Sol -2*T_ramp_Sol 0 T_ramp_Sol t_add ...
      t_add+discharge_time+20 t_add+discharge_time+T_ramp_Sol+20]*1e-3;                   %[s]
 %%%%%% 
 
 %Points 3 and 4 are (0,I_Sol_start), (t,ISol_to). Lets find the equation
 %for that line, and create the point 5 so it has the same slope
 
 [coef]=polyfit([0 T_ramp_Sol],[I_Sol_max 0],1)
 
 %If t5=something, we get I_5 should be:
 I_5=coef(1)*t_add+coef(2); %=I_Sol_max, como debe ser. Por eso, uso I_max.!!
 
 %%%
 
%Construct Sol, PF/Div coil current waveforms vertices
											  %!Breakdown!					%!Efit Icoil!
        %Time   	     [1,    2,              3,           4,        5,                  6,               7];
ISol_Waveform =  [0,  I_Sol_max, I_Sol_max, 0, -I_Sol_max, -I_Sol_max*0.95,   0];
IPF1_Waveform =  [0,  NaN,        NaN,        NaN,           I_PF1_Equil,   I_PF1_Equil,   0];
IPF2_Waveform =  [0,  NaN,        NaN,        NaN,           I_PF2_Equil,   I_PF2_Equil,   0];
%IDiv1_Waveform = [0,  NaN,        NaN,        NaN,           I_Div1_Equil,  I_Div1_Equil, 0];
IDiv1_Waveform = ISol_Waveform; %Div1 in series with Sol
IDiv2_Waveform = [0,  NaN,        NaN,        NaN,           I_Div2_Equil,  I_Div2_Equil,  0];
%%%%%
CoilWaveforms = [ISol_Waveform; IPF1_Waveform; IPF2_Waveform; IDiv1_Waveform; IDiv2_Waveform]


%%  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@CONFIGURATION OF FIESTA@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%%% Fiesta grid %%%%%%%%%%%%%%%%%%%%
%The grid where FIESTA will run
%Modifying the grid will get different things

R_simulation_limits = [0.03 1];                                                    %R limits of the grid(has to contain vessel)
Z_simulation_limits = [-1.3 1.3]+[mod_Z_neg mod_Z_pos];         %Z limits of the grid(has to contain vessel)
Grid_size_R = 200;                                                                       %Grid points in R, ST25D based
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
Ip=100e3;                                    %[A] plasma current PHASE 1$%&/
RGeo=0.45;                                %[m] This is Rmajor in the excel. From the eq plot!
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
control = fiesta_control( 'diagnose',true, 'quiet',false, 'convergence', 1e-5, 'boundary_method',2 );

%%  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@TARGET EQUILIBRIA@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
config  = fiesta_configuration( 'STV2C2', Grid, coilset);

icoil_equil = fiesta_icoil(coilset);

%Assign equilibrium coil currents to icoil object [kA]
icoil_equil.Sol=ISol_equil;         %Solenoid Equilibrium Current at time(5,6)
icoil_equil.PF1=CoilWaveforms(2,6);	%PF1 Equilibrium Current at time(5,6)
icoil_equil.PF2=CoilWaveforms(3,6);	%PF2 Equilibrium Current at time(5,6)
%icoil_init.Div1=CoilWaveforms(4,6);	%Div1 Equilibrium Current at time(5,6)
icoil_equil.Div1=icoil_equil.Sol;         %Div1 in series with Sol
icoil_equil.Div2=CoilWaveforms(5,6);	%Div2 Equilibrium Current at time(5,6)

%To do the equilibria calc, we can use EFIT algorithm or not:

    %1) NO EFIT..........
    
    %equil = fiesta_equilibrium( 'STV2C2', config, Irod, jprofile, control, [],icoil_equil ); 


    %2)EFIT..........
    %this calculates PF and Divs current given plasma parameters
    %Discovered by Juanjo Toledo Garcia
    
    %[efit_config, signals, weights, index]=fiesta_efit_configuration(config, {'PF1','PF2'}, [0.44, 0, 0.44/1.85 1.8 0.2])
    [efit_config, signals, weights, index]=efit_shape_controller(config, {'PF2'}, [0.44, 0, 0.44/1.85 1.8 0.2])
    
    % The numbers you give are [Rgeo, Zgeo, a, kappa, delta], Rgeo,Zgeo,a are
    % mandatory.
    %I use the values of the standar shape, to get a similar equil

    equil=fiesta_equilibrium('ST', config, Irod, jprofile, control,efit_config, icoil_equil, signals, weights) %%EFIT!!!
    %It does the case in line 96!! The equil calc is in lin 124

    %Now we have to extract the new currents from the equil, provided that EFIT
    %changed some of them to satisfy the conditions requested:
    icoil_equil=get(equil,'icoil');                         %redefine icoil to save the new current values
    current_post_EFIT=get(icoil_equil,'currents');
    CoilWaveforms(iPF1,5:6) =current_post_EFIT(iPF1);
	CoilWaveforms(iPF2,5:6) =current_post_EFIT(iPF2);

%%%%%%END EFIT................

%%%%%%  Fiesta Plot  %%%%%%%%%%%%%%%%%

    figure; hold on; axis equal;
    plot(coilset);
    contour( get(equil,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
    contour( get(equil,'Psi'),get(equil,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
    plot(vessel);
    fileName = 'ST_target_equilibrium';
    legend(gca,'hide');
    %set(gca,'XLim',[0 1]);
    %set(gca,'YLim',[-1.5 1.5]);
    xlabel(gca,'R (m)');
    ylabel(gca,'Z (m)');
    title('Target equilibria for SST phase 2');
    % save_to_pdf( gcf, fileName );
    %%%OPtionf for tfg
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    axis([0,1.1,-1.1,1.1]) 
    %print -depsc2 NOMBREPLOT.eps
    
 %Plot from the demo in examples of FIESTA folder
        section_figure=section(equil); %THIS IS A PLOT
        
        figure;
        plot(equil)
        parametersshow(equil)           %this plots the parameters in the equil
        hold on
        plot(vessel)
        plot(coilset)
        title('Target equilibria phase 2')
            %Watch out, R0 in the plot is r0_mag, not r0_geom!!
        Filename = '_TargetEquilibrium';
        saveas(gcf, strcat(ProjectName,Filename,FigExt));
%Equilibrium parameters

parameters(equil); %its better in this way, because it shows units and coil currents. If you define a variable, it wont do that
param_equil=parameters(equil)                             %Will be used in the null sensors

%%%%%%  Vessel time constant  %%%%%%%%%%%
%This has to deal with the change of eddy currents in the vessel (See the
%article). However, this is not necessary (at least, yet), so I wont use
%it.
% [~,tau_vessel] = eig(curlyR(1:end-3,1:end-3)\curlyM(1:end-3,1:end-3));
% tau_vessel = max(diag(tau_vessel));
% disp([ 'tau_vessel=' num2str(tau_vessel*1e3) 'ms' ]);


%%% Make virtual sensors where we want breakdown  %%%%%%%
%This is to null the poloidal field(BP), to increase the connective length, and
%allow the plasma breakdown. Juanjo created it in a small region, setting
%the max and min intervals. What I will do, for simplicity, is create them
%using RGeo and Zgeo, The plasma central positions. 
%I will magnify the region, since it is reasonable that the plasma will be
%created in a greater region than what Juanjo uses

BP_virt_R = linspace(param_equil.r0_geom-a_eff,param_equil.r0_geom+a_eff,10);   %R values, 100 values
BP_virt_Z = linspace(ZGeo-a_eff,ZGeo+a_eff,10);             %Z values, 100 values

[BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
BP_virt_R = BP_virt_R(:)';
BP_virt_Z = BP_virt_Z(:)';

BP_virt_theta = zeros(1,length(BP_virt_R));
nSensors = length(BP_virt_theta); %100

BP_virt_names = {};
for iSensor=1:nSensors
    BP_virt_names{iSensor} = ['Radial Bp Virtual Sensor #' num2str(iSensor) ];
end

BP_virt_R = [BP_virt_R  BP_virt_R];
BP_virt_Z = [BP_virt_Z  BP_virt_Z];         %Both size 1*200. It is replicated, so element 101=element 1 

BP_virt_theta = [BP_virt_theta  BP_virt_theta+pi/2];        %size 1*200. The first 100 have 0, and the second has pi/2 

%Taken from a sensro Btheta function:
    %     theta=0.00 --> sensor is pointing in the R direction
    %     theta=pi/2 --> sensor is pointing in the Z direction

for iSensor=nSensors+1:2*nSensors
    BP_virt_names{iSensor} = ['Vertical Bp Virtual Sensor #' num2str(iSensor) ];
end

sensor_btheta = fiesta_sensor_btheta( 'sensor', BP_virt_R, BP_virt_Z,BP_virt_theta, BP_virt_names );

%(r,z) of the sensors

r_sensors=get(sensor_btheta,'r'); %size 1*200
z_sensors=get(sensor_btheta,'z'); %size 1*200
[R_sensors,Z_sensors]=meshgrid(r_sensors,z_sensors); %size 200*200

%Plot of the sensors
    figure; hold on; axis equal;
    plot(coilset);
    contour( get(equil,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
    contour( get(equil,'Psi'),get(equil,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
    plot(vessel);
    fileName = 'ST_target_equilibrium';
    legend(gca,'hide');
    plot(sensor_btheta);
    %set(gca,'XLim',[0 1]);
    %set(gca,'YLim',[-1.5 1.5]);
    xlabel(gca,'R (m)');
    ylabel(gca,'Z (m)');
    title('Target equilibria phase 2 with sensors to null Bpol');
    % save_to_pdf( gcf, fileName );
    %%%OPtionf for tfg
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    axis([0,1.1,-1.1,1.1]) 

%%%%END OF FIESTA EQ@@@@@@@@@@@@@@@@@@@@@@@@@@

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@RZIp@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% It has the virtual sensors

rzip_config = fiesta_rzip_configuration( 'RZIP', config, vessel, {sensor_btheta} );
[A, B, C, D, curlyM, curlyR, gamma, plasma_parameters, index, label_index, state] = response(rzip_config, equil, 'rp',plasma_resistance);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        %D1=C_temp(:,iPF1:end);          %Elements of C_temp(Cn) related to the PF and Div coils
        %I_PF_null = -pinv(D1) * (C1*ISol_Waveform(2));	%(PF1,PF2,Div1,Div2)

    %2)Div1 in serie with Sol! :
        %Since we want Div1 in serie with Sol, it should be reomved from the matrix
        %to do the calc:
        D1=[D1_PF1 D1_PF2 D1_Div2];                     %This is for Div1 in series with Sol, so 
                                                                    %must not include it in the calculation for the currents                                                                                                                 
        I_PF_null = -pinv(D1) * (C1*ISol_Waveform(2)+D1_Div1*ISol_Waveform(2));	%Div1 in serie with Sol
        I_PF_null=[I_PF_null(1); I_PF_null(2); ISol_Waveform(2); I_PF_null(3)]; %To not have problems
        
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
    title('I_{{input}} versus time')
    legend('Sol','PF1','PF2','Div1','Div2')
    %%%optinos for tfg
    set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
    set(gca, 'FontSize', 13, 'LineWidth', 0.75);                    %<- Set properties TFG
    Filename = '_InputCurrents';
    saveas(gcf, strcat(ProjectName,Filename,FigExt));



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


%% %%%%%RZIp run%%%%%%%%%%%%%
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

        %%%I_PF_output and plasma current
            figure;
            subplot(3,1,1)
            %plot(time_adaptive*1e3,I_PF_output/(1e3))
            %hold on
            plot(time_adaptive*1e3,Ip_output/(1e3))
            ylabel('I (kA)')
            title('Dynamic response SST phase 2')
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
            saveas(gcf, strcat(ProjectName,Filename,FigExt));
          %print -depsc2 NOMBREPLOT.eps

%%% END RZIP@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%% RE-DOING EQUILIBRIA CALC WITH EDDYS

    %Now that we havr compute the eddys, we could re do all the calc, to
    %get an equil with those eddys, with the new equil do the RZip again,
    %and, it the new eddys do not change much, could accept the result.
    %Yes, this is a iterative process, but with a huge consumption of
    %computer resources==> :)

    %To save some memory, will rewrite things:
      %coilset_noVV=coilset; %the old coilset, without VV
      coilsetVV=fiesta_loadassembly(coilset, vessel);
      configVV = fiesta_configuration( 'SMART with VV', Grid, coilsetVV);
     
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
%This is to null the poloidal field(BP), to increase the connective length, and
%allow the plasma breakdown. 

BP_virt_R = linspace(param_equil.r0_geom-a_eff,param_equil.r0_geom+a_eff,10);       %R values, 100 values
BP_virt_Z = linspace(ZGeo-a_eff,ZGeo+a_eff,10);                     %Z values, 100 values

[BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
BP_virt_R = BP_virt_R(:)';
BP_virt_Z = BP_virt_Z(:)';

BP_virt_theta = zeros(1,length(BP_virt_R));
nSensors = length(BP_virt_theta); %100

BP_virt_names = {};
for iSensor=1:nSensors
    BP_virt_names{iSensor} = ['Radial Bp Virtual Sensor #' num2str(iSensor) ];
end

BP_virt_R = [BP_virt_R  BP_virt_R];
BP_virt_Z = [BP_virt_Z  BP_virt_Z];         %Both size 1*200. It is replicated, so element 101=element 1 

BP_virt_theta = [BP_virt_theta  BP_virt_theta+pi/2];    %size 1*200. The first 100 have 0, and the second has pi/2 

for iSensor=nSensors+1:2*nSensors
    BP_virt_names{iSensor} = ['Vertical Bp Virtual Sensor #' num2str(iSensor) ];
end

sensor_btheta = fiesta_sensor_btheta( 'sensor', BP_virt_R, BP_virt_Z,BP_virt_theta, BP_virt_names );

%(r,z) of the sensors
r_sensors=get(sensor_btheta,'r'); %size 1*200
z_sensors=get(sensor_btheta,'z'); %size 1*200
[R_sensors,Z_sensors]=meshgrid(r_sensors,z_sensors); %size 200*200

%Plot of the sensors
    figure; hold on; axis equal;
    plot(coilset);
    contour( get(equil_eddy,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
    contour( get(equil_eddy,'Psi'),get(equil_eddy,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
    plot(vessel);
    fileName = 'ST_target_equilibrium';
    legend(gca,'hide');
    plot(sensor_btheta);
    %set(gca,'XLim',[0 1]);
    %set(gca,'YLim',[-1.5 1.5]);
    xlabel(gca,'R (m)');
    ylabel(gca,'Z (m)');
    title('Sensors to null Bpol with eddys');
    % save_to_pdf( gcf, fileName );
    %%%OPtionf for tfg
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    axis([0,1.1,-1.1,1.1]) 

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
Rfil = get(ftmp(:),'r'); %dim 1*number of filaments; number of filaments=nf
Zfil = get(ftmp(:),'z'); %dim 1*number of filaments
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

figure;
 scatter3(Rfil,Zfil,I_Passive_fil/(1e3),100,I_Passive_fil/(1e3),'filled')
 hold on
 plot(coilset)
  view(2) %2D view
 colorbar %colorbar
xlabel(' R (m)')
ylabel('Z (m)')
%zlabel('I (A)')
axis([0,1.03,-1.1,1.1]) %for the tfg EDDY Y FORCES
title('Eddy current in the vessel in kA')
set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
grid off
print -depsc2 NOMBREPLOT.eps

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

%%%Field
Br=get(equil,'Br'); %Here is it a fiesta_field,containig Br, in all points of the grid
Bz=get(equil,'Bz'); 
Bphi=get(equil,'Bphi'); 

Brdata = get(get(equil,'Br'),'data'); %Note this is 200*251, GridR*GridZ dimension
Bzdata = get(get(equil,'Bz'),'data');
Bphidata = get(get(equil,'Bphi'),'data');

%%Field values, sized Grid_size_R*Grid_size_Z, as expected.
%However, It is neccesary to reshape it in orden to do the interpolation:
Brdata = reshape( Brdata, length(zGrid), length(rGrid)); %251*200
Bzdata = reshape( Bzdata, length(zGrid), length(rGrid)); %251*200
Bphidata = reshape( Bphidata, length(zGrid), length(rGrid)); %251*200

%this are the value of the field in the grid points, but we want its value
%at the vessel. How do we calculate that? The easiest way (I wonder
%wether it can be done by another way) is to interpolate the data; Carlos came
%up with:

Br_interp = @(r,z) interpn(zGrid,rGrid,Brdata,z,r);
Bz_interp = @(r,z) interpn(zGrid,rGrid,Bzdata,z,r);
Bphi_interp = @(r,z) interpn(zGrid,rGrid,Bphidata,z,r);

%So the field at the vessel is:

Br_vessel=Br_interp(Rfil,Zfil);
Bz_vessel=Bz_interp(Rfil,Zfil);
Bphi_vessel=Bphi_interp(Rfil,Zfil);
B_vessel=[Br_vessel' Bphi_vessel' Bz_vessel']; 
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

%plot
figure; 
 scale_factor=2*10^5; %graphic needs to be scaled
quiver(xaccum,yaccum,stress_R/scale_factor,stress_Z/scale_factor,'color',[1 0 0],'AutoScale','off')
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

%%

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@BREAKDOWN CALC@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%%%Field%%%%%%%
%To compute the field, will do an equ calc, with no plasma, and the
%breakdown currents

    %I will use here the fiesta_load_assembly to include the eddys of the
    %vessel
    
coil_currents_break = zeros(1,nPF+get(vessel,'n'));                       %n=308, number of filaments
coil_currents_break(iSol) = I_Sol_max;
coil_currents_break(2:nPF) = I_PF_null';
coil_currents_break(iDiv1) =coil_currents_break(iSol); %Div1 in series with Sol

%To set the eddy currents for the vessel, will select the maximum eddy
%between a certain time interval. Will chose to determine the interval, the
%plasma current, of course. Since the eddy increases as Ip icnreases, will
%fix the limit on Ip, and select those values.

    %1) Select time>0
    
    index_time_Iplimit=time_adaptive>=0 & time_adaptive<=10e-3;             %indexes within the desired interval
    time_Iplimit=time_adaptive(index_time_Iplimit);

    [MAX,index_maxIPassive_VV]=max(I_Passive_VV(index_time_Iplimit));       %searching the max IPassive (sum over all the filaments)
    I_Passive_Iplimit=I_Passive(index_maxIPassive_VV,:);            %chosing the values of the filaments

    %now will choose the one with the maximum total eddy current

coil_currents_break(nPF+1:end)=I_Passive_Iplimit;                 %eddy currents, have to define each filament

icoilVV_break = fiesta_icoil( coilsetVV, coil_currents_break );
equil_optimised_null = fiesta_equilibrium( 'optimised null', configVV, Irod, icoilVV_break );
    
    %Plot
    figure;
    plot(equil_optimised_null);
    hold on;
    plot(coilset)
    fileName = '\Psi contour optimised_null';
    legend(gca,'hide');
    title('\Psi contours poloidal null phase 2 XL')


%%Field inside vessel (breakdown)

%On the forces section I get the field at each point of the grid. I also
%ahve the (R,Z) of the filaments of the vessel, so I should get (R,Z) of
%the points inside, and with them get the field.
%Since I know the limit, I can create a linspace with the values of R and Z
%of the number of points I desired, and then interpolate

n_pnts_inside=100; %100 is the ideal to have good plots of the fields, but the 
    %L int failures. Need to be 30 for the integration.

r_inside_VV=linspace(VesselRMinPoint,VesselRMaxPoint,n_pnts_inside); 
z_inside_VV=linspace(VesselZMinPoint,VesselZMaxPoint,n_pnts_inside);

%Have to do ameshgrid for the interpolation, so each R value has a Z value
[R_in,Z_in]=meshgrid(r_inside_VV,z_inside_VV);

%Comprobation!. To plot it, I need to do a meshgrid
%     figure;
%     plot(R_in,Z_in,'b.')
%     hold on;
%     plot(vessel)
%It seems right, so awesome.


%With the grid, now I can get the field easily, as I did on the stresses
%calc. I need the field during the breakdown phase, which means without
%plasma, in the optimized null setup. Doing similarly as in the stresses
%section. 
    
    %%%Field in the breakdown phase
    Br_break=get(equil_optimised_null,'Br'); 
    %Here is it a fiesta_field,containig Br, in all points of the grid
    Bz_break=get(equil_optimised_null,'Bz'); 
    Bphi_break=get(equil_optimised_null,'Bphi_vac'); 
    %Dont know why, but Bphi do not work, I have to use Bphi_vac. However, by
    %seeing its graph, it has from linear tinterpolation, BT=0.1 for
    %Rgeo=0.4993, which is virtually RGeo from the stimations, 0.45, so it
    %shouldnt be wrong (or so wrong)

    Br_breakdata = get(get(equil_optimised_null,'Br'),'data'); 
    %Note this is 200*251, GridR*GridZ dimension
    Bz_breakdata = get(get(equil_optimised_null,'Bz'),'data');
    Bphi_breakdata = get(get(equil_optimised_null,'Bphi_vac'),'data');

    %%Field values, sized Grid_size_R*Grid_size_Z, as expected.
    %However, It is neccesary to reshape it in orden to do the interpolation:
    Br_breakdata = reshape( Br_breakdata, length(zGrid), length(rGrid)); %251*200
    Bz_breakdata = reshape( Bz_breakdata, length(zGrid), length(rGrid)); %251*200
    Bphi_breakdata = reshape( Bphi_breakdata, length(zGrid), length(rGrid)); %251*200
    
%this are the value of the field in the grid points, but we want its value
%at the vessel. How do we calculate that? The easiest way (I wonder
%wether it can be done by another way) is to interpolate the data; Carlos came
%up with:

Br_break_interp = @(r,z) interpn(zGrid,rGrid,Br_breakdata,z,r,'mikama');
Bz_break_interp = @(r,z) interpn(zGrid,rGrid,Bz_breakdata,z,r,'mikama');
Bphi_break_interp = @(r,z) interpn(zGrid,rGrid,Bphi_breakdata,z,r,'mikama');

%I need pairs of points (R_in,Z_ins), and for this I use the mesh

Br_ins_vessel=Br_break_interp(R_in,Z_in);
Bz_ins_vessel=Bz_break_interp(R_in,Z_in);
Bphi_ins_vessel=Bphi_break_interp(R_in,Z_in);
%B_ins_vessel=[Br_ins_vessel' Bphi_ins_vessel' Bz_ins_vessel'];

%min(min(Bz_ins_vessel))
%min(min(Br_ins_vessel))
%This value have not been altered by the addition of the Earth field

Bpol_ins_vessel=sqrt(Br_ins_vessel.^2+Bz_ins_vessel.^2);

%Field in the sensor region %%%%%%%
%The same applies here, I know the R and Z values, so making a linspace:

r_inside_sensors=linspace(min(r_sensors),max(r_sensors),30);
z_inside_sensors=linspace(min(z_sensors),max(z_sensors),30);

%Meshgrid
[R_sensorsG,Z_sensorsG]=meshgrid(r_inside_sensors,z_inside_sensors);
%Comprobation!. To plot it, I need to do a meshgrid
%     figure;
%     plot(R_sensorsG,Z_sensorsG,'b.')
%     hold on;
%     plot(vessel)
%     plot(sensor_btheta)
%It seems right, so awesome.
%...........................................
%Field
Br_sensor=Br_break_interp(R_sensorsG,Z_sensorsG);
Bz_sensor=Bz_break_interp(R_sensorsG,Z_sensorsG);
Bphi_sensor=Bphi_break_interp(R_sensorsG,Z_sensorsG);

Bpol_sensor=sqrt(Br_sensor.^2+Bz_sensor.^2);

% min(min(Bpol_sensor))

%%%%PLOTS   

% % %%Quiver plot
% %  scale_factor=1e-1; %graphic needs to be scaled
% % figure; 
% % subplot(1,3,1)
% % quiver(R_in,Z_in,Br_ins_vessel/scale_factor, Bz_ins_vessel/scale_factor,'color',[1 0 0],'AutoScale','off')
% % hold on;
% % plot(vessel)
% %    view(2) %2D view
% %       plot(sensor_btheta)
% % xlabel('R (m)')
% % ylabel('Z (m)')
% % title('\vec{B_tokamak} inside vessel')
% % subplot(1,3,2)
% % quiver(R_in,Z_in,BrEarthVV/scale_factor, BzEarthVV/scale_factor,'color',[1 0 0],'AutoScale','off')
% % hold on;
% % plot(vessel)
% %    view(2) %2D view
% %       plot(sensor_btheta)
% % xlabel('R (m)')
% % ylabel('Z (m)')
% % title('\vec{B_Earth} inside vessel')
% % subplot(1,3,3)
% % quiver(R_in,Z_in,Br_ins_vesselE/scale_factor, Bz_ins_vesselE/scale_factor,'color',[1 0 0],'AutoScale','off')
% % hold on;
% % plot(vessel)
% %    view(2) %2D view
% %       plot(sensor_btheta)
% % xlabel('R (m)')
% % ylabel('Z (m)')
% % title('\vec{B_tokamak}+\vec{B_Earth} inside vessel')
% 

%  % Plot mod Br   
%     figure; 
%     surf(R_in,Z_in,log10(abs(Br_ins_vessel)),'FaceAlpha',0.5,'EdgeColor','none');
%     hold on;
%     plot(vessel)
%     colorbar %colorbar
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     %plot(sensor_btheta)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('log_{10}(Br) inside vessel ')
%        
%   % Plot Br   
%     figure; 
%     surf(R_in,Z_in,Br_ins_vessel,'FaceAlpha',0.5,'EdgeColor','none');
%     hold on;
%     plot(vessel)
%     colorbar %colorbar
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     %plot(sensor_btheta)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Br ') 
%     
% % Plot mod Bz
%     figure; 
%     surf(R_in,Z_in,log10(abs(Bz_ins_vessel)),'FaceAlpha',0.5,'EdgeColor','none');
%     hold on;
%     plot(vessel)
%     colorbar %colorbar
%      view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('log_{10}(Bz) inside vessel')
%   
% 
%  % Plot Bz   
%     figure; 
%     surf(R_in,Z_in,Bz_ins_vessel,'FaceAlpha',0.5,'EdgeColor','none');
%     hold on;
%     plot(vessel)
%     colorbar %colorbar
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     %plot(sensor_btheta)
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Bz ') 
    
 %Plot Bpol y Bphi   
    figure; 
    subplot(1,2,1)
    surf(R_in,Z_in,log10(abs(Bpol_ins_vessel)),'EdgeColor','none');
    shading('interp') %this is to make the transition between values continuous,
    %instedad of discontinuously between pixels
    hold on;
    plot(vessel)
    colorbar %colorbar
    view(2) %2D view
    plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
     xlabel('R (m)')
    ylabel('Z (m)')
    title('log_{10}(Bpol) phase2')
    subplot(1,2,2)
    surf(R_in,Z_in,Bphi_ins_vessel,'EdgeColor','none');
    shading('interp') %this is to make the transition between values continuous,
    %instedad of discontinuously between pixels    
    hold on;
    plot(vessel)
    colorbar %colorbar
    view(2) %2D view
    plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
    xlabel('R (m)')
    ylabel('Z (m)')
    title('Bphi')

%% %%L CALC FORMULAE%%%%%%%%%%%%

%The field in the whole vessel will be used, since it is the field that is
%used in the field line integrator.

%Minimun with veesel interp 
Bpolmin=min(min(Bpol_ins_vessel))                                               %minimun of poloidal field
[Bpolmin_index_row Bpolmin_index_column]=find(Bpol_ins_vessel==Bpolmin);            %indexes  

        %%%%%THIS IS TO COMPARE
%     %Minimun with sensors interp
% Bpolmin=min(min(Bpol_sensor)) %minimun of poloidal field
% [Bpolmin_index_row Bpolmin_index_column]=find(Bpol_sensor==Bpolmin);
%     %indexes
        %%%%%THIS IS TO COMPARE

        
%And the toroidal field at the center of the null region is:
Bphi_centerNull=Bphi_break_interp((min(r_sensors)+max(r_sensors))/2,...
    (min(z_sensors)+max(z_sensors))/2)

%The connective length without averaging is
L_no_aver=0.25*a_eff*Bphi_centerNull/Bpolmin                            %[m]

%A more precise way to calculate L is by averaging in the poloidal field.
%As stated by TCV thesis and the newer ITER article, Bpol is obtained by
%averaging on the surface of the null region. The toroidal field is not
%averaged. To do this, if I use the sensor field, it will be much much more
%easier, so I will do it.

Bpolmin_av=mean(mean(Bpol_sensor));              %to compute the mean inside the sensor region
L_aver=0.25*a_eff*Bphi_centerNull/Bpolmin_av                %[m] L with the average pol field

%% L CALC INTEGRATION%%%%%%%

%Now that we have the field, we want to comput the field lines, whose eq is
%cross(B,dr)=0 ==> Br/dr=Bphi/(r dphi)=Bz/dz =cte=tt
%I have to create a function. I dont want the constant tt, and I can remove
%it by diving the 3 eq Br=ctt*dr, Bphi=ctt*r*dphi, Bz=ctt*dz between each other
%One of the coordinates have to be the independent one, and phi seems the
%best option, since is >0 while the rest are ><0, and phi increase (or
%decrease) when you follow a field line. (Carlos Soria)

phi_values=1000;

n_iter=1000                                         %integer, Number of iterations!!!!

phiSpan=linspace(0,2*pi*n_iter,phi_values*n_iter);    %the range of values of phi to do the integration. remember phi is
                                        %the independent variable. To do more than 1 loop, we multiply
                                        %by an integer, n_iter

L_max=1000;                                                 %max value for the integration; when L achieves
                                                            %this value, the integration stops. ST have around 50m.

odefun= @(phi, rzL) Field_LineIntegrator(phi,rzL,Br_break_interp,...
    Bz_break_interp,Bphi_break_interp);
event_colission_wall=@(phi,rzL)Colission_wall(phi,rzL,VesselRMaxPoint,...
    VesselRMinPoint,VesselZMaxPoint,VesselZMinPoint,Br_break_interp,...
    Bz_break_interp,Bphi_break_interp,L_max); 
options = odeset('OutputFcn',@ode_progress_bar,'Events',event_colission_wall,'AbsTol',1e-10,'RelTol',1e-6); 
                                    %I include a fiesta funciton to show the progress of the ode

%Int eh Lazarus paper 1998, they compute connective length by avergaing on
%9 lines, the line with Bpol min, and the 8 surroundings. I will do the
%same. No need for loop, to avoid problems. I'll do it manually xD

%As suggested by Manolo, i will calculate L in thw whole grid inside the
%vessel. For doing this, provided that i already ahve the inner points,
%I can do it easily

%I redefine the grid to compute the connection length, for less computer
%demands (time)

n_pnts_insideL=15 %10%25 for a_eff=0.15;            %100 is the ideal to have good plots of the fields, but the L int failures. 

r_inside_VVL=linspace(VesselRMinPoint,VesselRMaxPoint,n_pnts_insideL); 
z_inside_VVL=linspace(VesselZMinPoint,VesselZMaxPoint,n_pnts_insideL);

    %Plot
    [r_ins_VVL,z_ins_VVL]=meshgrid(r_inside_VVL,z_inside_VVL);
    figure;
    plot(r_ins_VVL,z_ins_VVL,'r.')
    hold on
    plot(vessel)
    xlabel('R (m)')
    ylabel('Z (m)')

%Points inside, without the extremal points. This will be used in the ode45
r_inside_VV_noLimits=r_inside_VVL(2:end-1);
z_inside_VV_noLimits=z_inside_VVL(2:end-1);

r0_z0_L0=[0 0 0]; %the third column is always zero, since L starts at 0

tic

for i=1:length(r_inside_VV_noLimits)

    for j=1:length(z_inside_VV_noLimits)
        
        points=[r_inside_VV_noLimits(i) z_inside_VV_noLimits(j)];
        
        r0_z0_L0=[ r0_z0_L0; points 0];
        
    end
    
end

%I have the additional point 0 0 0 form the begining, that i can remove
%easily with
r0_z0_L0=r0_z0_L0(2:end,:);

%So, gotta solve 10^2*10^2=10^4 eq xD. Since I am only interested in L,
%could do a for loop simply saving the L value. Will do that:

for i=1:length(r0_z0_L0)
    fprintf('Iter %d de %d',i,length(r0_z0_L0))
     [phi_fieldline, rzL_fieldline]=ode45(odefun,phiSpan,r0_z0_L0(i,:),options);        %ode15s Carlos
    
    %To save the last values of R,Z,L
    RZL_lastTEST(i,1)=rzL_fieldline(end,1);
    RZL_lastTEST(i,2)=rzL_fieldline(end,2);
    RZL_lastTEST(i,3)=rzL_fieldline(end,3);                                                     %this contains L
    
end

%To store start points that do not collide: first I get the index of both R
%and Z, but together, since they do not collide if oth R and Z are greater
%than the min value, and lower than the greatest value
    
    RZ_store_index=RZL_lastTEST(:,1)<VesselRMaxPoint &...
        RZL_lastTEST(:,1)>VesselRMinPoint & RZL_lastTEST(:,2)<VesselZMaxPoint &...
        RZL_lastTEST(:,2)>VesselZMinPoint;

    RZ_no_collide=[r0_z0_L0(RZ_store_index,1) r0_z0_L0(RZ_store_index,2)];
    
    
%Plot    

%     figure;
%     ss=scatter3(r0_z0_L0(:,1),r0_z0_L0(:,2),RZL_lastTEST(:,3),400,L_int,'filled','MarkerFaceAlpha',1)
%     hold on
%     plot(RZ_no_collide(:,1),RZ_no_collide(:,2),'m*','MarkerSize',5)
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     hh=plot(vessel);
%     set(get(get(hh(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     legend('L (m)','Not colliding starting points','Field null region')  
%     view(2) %2D view
%     colorbar %colorbar
%     xlabel(' R (m)')
%     ylabel('Z (m)')
%     %axis([0,1.03,-1.1,1.1]) 
%     title('L (m) depending on the starting point for the integration')
%     set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG


%%Plot of one of the line

figure;
plot3(r0_z0_L0(end,1),0,r0_z0_L0(end,2),'k*','LineWidth',3)
hold on;
plot3(vessel)
plot3(coilset)
plot3(rzL_fieldline(:,1).*cos(phi_fieldline),rzL_fieldline(:,1).*sin(phi_fieldline),...
    rzL_fieldline(:,2),'r','LineWidth',3)
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');  
hold on
plot3(rzL_fieldline(length(rzL_fieldline),1).*cos(phi_fieldline(length(rzL_fieldline)))...
    ,rzL_fieldline(length(rzL_fieldline),1).*sin(phi_fieldline(length(rzL_fieldline))),...
    rzL_fieldline(length(rzL_fieldline),2),'g*','LineWidth',3)
title('Field line integration for the calculation of the connective length 3')
set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%legend('Starting point (Point with less Bpol)','Field line',...
%'End point (wall collision')



   %Plot shading
    figure; 
    omega_II_colormap                                                       %colormap used in the group
    tri = delaunay(r0_z0_L0(:,1),r0_z0_L0(:,2));
    plot(r0_z0_L0(:,1),r0_z0_L0(:,2),'.')
    %
    % How many triangles are there?
    [r,c] = size(tri);
    disp(r)
    % Plot it with TRISURF
    trisurf(tri, r0_z0_L0(:,1), r0_z0_L0(:,2),RZL_lastTEST(:,3),'FaceAlpha',1,'EdgeColor','none');
    %view(2)
    colorbar
    hold on
    shading('interp')                       %this is to make the transition between values continuous,
                                                        %instedad of discontinuously between pixels
    plot3(RZ_no_collide(:,1),RZ_no_collide(:,2),...
        max(RZL_lastTEST(:,3))*ones(length(RZ_no_collide(:,1)),1),...
        'm*','MarkerSize',7)
    plot3([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
    [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],...
    ones(1,5)*max(RZL_lastTEST(:,3)),'r.--')
    hh=plot(vessel);
    set(get(get(hh(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('L (m)','Not colliding starting points','Field null region')  
    view(2) %2D view
    colorbar
    %caxis([min(RZ_no_collide(:,3)),max(RZ_no_collide(:,3))])
    xlabel(' R (m)')
    ylabel('Z (m)')
    %axis([0,1.03,-1.1,1.1]) 
    title('L (m) depending on the starting point for the integration')
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        %Note that I convert 2D plots into 3D to plot it at the plane z=max
        %value of L, because the L plot is a surf, so it has height. 
        
        
%%%%%%%%%5
%The event function is at the end!!
    

time_Ode=toc;     
%% PASCHEN CURVE PLOT
 
 p=linspace(1e-5,1e-2,10000);                           %[Tor]pressure of the prefill gas. size>100 because if
                                                                            %not, It do not work properly
 Emin= @(L,p) C_2*p./log(C_1*p*L);
 
 figure;
    loglog(p,Emin(10,p))
    hold on
    loglog(p,Emin(100,p),'r')
    loglog(p,Emin(1000,p),'g')
    %
    loglog(p,3*ones(1,length(p)))
    loglog(p,Emin(50,p))
    xlabel('Prefill pressure (Torr)')
    ylabel('E_{min} (V/m)')
    legend('L=10m','L=100m','L=1000m','VEST E (no L data)','Globus-M, L=50m, V_{loop}=4.5-8V')
    title(sprintf('Paschen curve, Gas=%s',Gas_type)); %d for numbers

%%%Estimation of the breakdown time

%Imagine the values choosen are
E_choosen=5 %V/m
p_choosen=10^-5 %Torr (multiplo de Pa)
K_B=1.380649*10^-23 %J K
T_prefill=298*K_B; %K, temperature of the prefill gas

tau_bd= @(percent) log(percent/100*p_choosen/T_prefill)/...
    (43*E_choosen/p_choosen*(C_1*p_choosen*exp(-C_2*p_choosen/E_choosen)*(1-0.3)));
        %has choosen L to be 0.3alpha

%tau_bd= @(percent) log(percent/100)/43;

percentage=linspace(2,100);
figure;
plot(percentage, tau_bd(percentage)*1e3)
xlabel('Percentage of ionized prefill gas (%)')
ylabel('Breakdown duration (ms)')
title('Breakdown duration (starts at t=0, finish at the y axis time')
    
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


%%
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Functions for the breakdown integration@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    %1) Field line integrator function
        %this solves the field line eq, using phi as the independent
        %variable, so the derviatvies are dr/dphi and dz/dphi
    
    
    function [results]=Field_LineIntegrator(phi,rzL,Br_break_interp,Bz_break_interp,Bphi_break_interp)
    %rzL=[r z L]
   
    %First, the field needs to be evaluated at the point (r,phi,z):
    
    Br_eval=Br_break_interp(rzL(1),rzL(2));
    Bphi_eval=Bphi_break_interp(rzL(1),rzL(2));
    Bz_eval=Bz_break_interp(rzL(1),rzL(2));
    
    %With the field, the eq to solve is:
    
    dr_dphi=rzL(1)*Br_eval/Bphi_eval;
    dz_dphi=rzL(1)*Bz_eval/Bphi_eval;
    length=rzL(1)*sqrt(1+(Br_eval/Bphi_eval)^2+(Bz_eval/Bphi_eval)^2);
    
    results=zeros(3,1); %column vector to group the results
    results(1)=dr_dphi;
    results(2)=dz_dphi;
    results(3)=length;
    
    end

    %2)EVENT FUNCTION FOR THE ODE

   %Only one function can be introduced into ode, so this event has to have
   %all the conditions
   
   function [rz_value isterminal direction]=Colission_wall(phi,rzL,...
       VesselRMaxPoint,VesselRMinPoint,VesselZMaxPoint,VesselZMinPoint,...
       Br_break_interp,Bz_break_interp,Bphi_break_interp,L_max) 
   
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
   
    Br_eval=Br_break_interp(rzL(1),rzL(2));
    Bphi_eval=Bphi_break_interp(rzL(1),rzL(2));
    Bz_eval=Bz_break_interp(rzL(1),rzL(2));
   
   L_lim=rzL(3)-L_max;                                                  %[m] Maximum L value=20
   
   rz_value=[up_colission down_colission out_colission in_colission L_lim];
   
    %Have checked that if I dont use the minR condition, it will impige 
    %in the upper wall, which was was happens at the beggining, when
    %I dont have the inner wall condition
   end

 