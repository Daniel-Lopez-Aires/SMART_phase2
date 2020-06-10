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

%#################################################
%######################TFM#########################
%#################################################

%%%PARAMETER TO BE VARIED FOR THE BREAKDOWN!!!!!!!!!!!!!!!!!!!!!
    %1)GAS TYPE

        %i)H    
        %Z_eff=Z_nucleus-Debye shielding. It can be calculated in
        %http://calistry.org/calculate/slaterRuleCalculator , giving:
             Gas_type='H'
             Z_eff=1                                                          %H
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
    
                Gr_fraction=0.7;%0.15;                              %This is to scale the Gr_limit==> <=1
    
    %3) a_eff                                                        %[m] This defines the size of the field null region.
   
             %a_eff=0.05;                                   % little null region WHAT JJ USED, ST25D BASED
             a_eff=0.15;                                        % large null region
            %a_eff=0.3
%%%%%%%%END PARAMETERS TO BE VARIED ON BREAKDWON!!!!!!!!!!!!!

%#################################################
%#####################TFM############
%#################################################


%%%%%%%%%%%%%%%%%  DEFINE DATA OUTPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%
%Taken from the God Scott

%Define figure extension
FigExt = '.png'; 		%'.png','.eps','.pdf'

%Define project and series names
ProjectName = 'SMART-P1-9L';			%Define global project name
%SeriesName = 'VaryTauP';		%Define parameter scan series name

%Create global output folders for saved data and figures
ASCIIDir = 'RawData/'; mkdir(ASCIIDir);
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


%Define Vessel Wall Thickness
VWall_Inboard=0.004;	%Inboard Wall Thickness		[m]
VWall_Outboard=0.008;	%Outboard Wall Thickness	[m]
VWall_Upper=0.015;		%Top Wall Thickness			[m]
VWall_Lower=0.015;		%Bottom Wall Thickness		[m]
Cross_filament=0.00015;           %[m^2] cross section of the VV's filaments

%This are values of the internal vessel. For the outer, the width and heigh
%must be considered

VesselRMinPoint=0.15+VWall_Inboard;;              %[m] R min position
VesselRMaxPoint=0.8;                                                %[m] R max position
VesselZMinPoint=-.8;                               %[m] Z min position 
VesselZMaxPoint=.8;                                %[m] Z max position 

ZMinCenter=VesselZMinPoint-VWall_Lower/2;
ZMaxCenter=VesselZMaxPoint+VWall_Upper/2; %+
RMinCenter=VesselRMinPoint-VWall_Inboard/2;
RMaxCenter=VesselRMaxPoint+VWall_Outboard/2; %+

%They do follow an order!
point1=[RMaxCenter ZMaxCenter ];                     %Top Right vertice
point2=[RMaxCenter ZMinCenter];                     %Bottom Right
point3=[RMinCenter ZMinCenter];                     %Bottom Left
point4=[RMinCenter ZMaxCenter];                     %Top Left

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
nZDiv1=6; nRDiv1=4;
nZDiv2=6; nRDiv2=4;
nZPF1=6; nRPF1=4;
nZPF2=6; nRPF2=4;

nDiv1=nZDiv1*nRDiv1;
nDiv2=nZDiv2*nRDiv2; 
nPF1=nZPF1*nRPF1; 
nPF2=nZPF2*nRPF2;  

%Dimensions of the PF coils
width_PF=75*1e-3;  %[m] Width of the PF coil (m)  0.111chonky not updated   (initially 0.042m) (diff 0.069 - to maintain relative edge)
height_PF=50*1e-3; %[m] Height of a the PF coil (m) 0.074 chonky updated	(Initially 0.035m) (diff 0.039 - to maintain relative edge)
        %This neglects the separation between conductors in real desgin!!
        
%Position of the coils (m)
R_PF1=0.9375;  Z_PF1=0.3075;  %[m] Position of PF1 %(0.90,0.30) before 4/6 update
R_PF2=0.9375;  Z_PF2=0.6075;  %[m] Position of PF2 	%(0.90,0.60) before 4/6 update
R_Div1=0.2355;  Z_Div1=0.890;  %[m] Position of Div1	%(0.2,0.86+0.039) before 4/6 update
R_Div2=0.4475;  Z_Div2=0.890;  %[m] Position of Div2  %(0.45,0.86+0.039) before 4/6 update

% Make Solenoid
nSol=210; %800;                                      % number of turns of the solenoid (Agredano suggested)

RSolInner = 0.115; RSolOuter=0.145;   % Inner and Outer solenoid radii    [m]
 RSol=(RSolInner+RSolOuter)/2; % Central radius of solenoid (0.13) [m]
 ZMax_Sol=ZMaxCenter+(VWall_Lower/2);                              % [m] Max Z position
 ZMin_Sol=-ZMax_Sol;                                                  % [m] Min Z position
 
%turns 
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

SolWidth = RSolOuter-RSolInner;                 %[m] - Previously sqrt(70e-6)
SolHeight = (2*ZMax_Sol)/(length(Z_filament));   %[m] - Previously sqrt(70e-6)

for iFilament=1:nfil_ind_coil
    coil_filaments(iFilament) = fiesta_filament( RSol,Z_filament(iFilament), SolWidth,SolHeight ); 
    %values by default, and R value defined so that the solenoid fit its
    %region
end

coil_1  = fiesta_coil( 'psh_coil', coil_filaments, 'Blue', resistivity, density );
Sol_circuit = fiesta_circuit( 'Sol', [1], [coil_1] );

%%%%Creation of the vessel object%%%%%%%%%%%%%%%

 %Lets create the lines manually
 
 %Each wall ahve a certain width
 
 %%%1->2) From Zmax to Zmin, with Rmax= point 1 to 2
 
 ww_R=VWall_Outboard;
 ww_Z=Cross_filament/ww_R;
 
 n_fil_Z=round((ZMaxCenter-ZMinCenter+2*ww_Z)/ww_Z); %8;         %Number of filaments in the Z direction, from Zmin to Zmax
 n_fil_R=round((RMaxCenter-RMinCenter+2*ww_R)/ww_R); %8;        %Number of filaments in the R direction, from Rmin to Rmax

    %From point 1 to point 2, R=Rmax
    R_lin1_2=RMaxCenter*ones(1,n_fil_Z);
    Z_lin1_2=linspace(ZMaxCenter,ZMinCenter,n_fil_Z);

    %Array with widths
    wwR_1_2=ww_R*ones(1,n_fil_Z);
    wwZ_1_2=ww_Z*ones(1,n_fil_Z);
    
 %%%2->3) From Rmax to Rmin, with Zmin= point 2 to 3
 
 ww_Z=VWall_Lower;
 ww_R=Cross_filament/ww_Z;
 
 n_fil_Z=round((ZMaxCenter-ZMinCenter+2*ww_Z)/ww_Z); %8;         %Number of filaments in the Z direction, from Zmin to Zmax
 n_fil_R=round((RMaxCenter-RMinCenter+2*ww_R)/ww_R); %8;        %Number of filaments in the R direction, from Rmin to Rmax
 
    %From point 2 to point 3, Z=Zmin
    R_lin2_3=linspace(RMaxCenter,RMinCenter,n_fil_R);
    Z_lin2_3=ZMinCenter*ones(1,n_fil_R);    

    %Array with widths
    wwR_2_3=ww_R*ones(1,n_fil_R);
    wwZ_2_3=ww_Z*ones(1,n_fil_R);  
    
 %%%3->4) From Zmin to Zmax, with Rmin= point 3 to 4
 
 ww_R=VWall_Inboard;
 ww_Z=Cross_filament/ww_R;
 
 n_fil_Z=round((ZMaxCenter-ZMinCenter+2*ww_Z)/ww_Z); %8;         %Number of filaments in the Z direction, from Zmin to Zmax
 n_fil_R=round((RMaxCenter-RMinCenter+2*ww_R)/ww_R); %8;        %Number of filaments in the R direction, from Rmin to Rmax
    
    
   %From point 3 to point 4, R=Rmin
     R_lin3_4=RMinCenter*ones(1,n_fil_Z);
     Z_lin3_4=linspace(ZMinCenter,ZMaxCenter,n_fil_Z);
    %Array with widths
    wwR_3_4=ww_R*ones(1,n_fil_Z);
    wwZ_3_4=ww_Z*ones(1,n_fil_Z);
 
    %%%4->1) From Rmin to Rmax, with Zmax= point 4 to 1
 
 ww_Z=VWall_Upper;
 ww_R=Cross_filament/ww_Z;
 
 n_fil_Z=round((ZMaxCenter-ZMinCenter+2*ww_Z)/ww_Z); %8;         %Number of filaments in the Z direction, from Zmin to Zmax
 n_fil_R=round((RMaxCenter-RMinCenter+2*ww_R)/ww_R); %8;        %Number of filaments in the R direction, from Rmin to Rmax    
     
    %From point 4 to point 1, Z=Zmax
    R_lin4_1=linspace(RMinCenter,RMaxCenter,n_fil_R);
    Z_lin4_1=ZMaxCenter*ones(1,n_fil_R);
    
    %Array with widths
    wwR_4_1=ww_R*ones(1,n_fil_R);
    wwZ_4_1=ww_Z*ones(1,n_fil_R); 
    
%       figure;
%       plot(R_lin1_2,Z_lin1_2,'ro');
%        xlabel('R (m)')
%        ylabel('Z (m)')
%       hold on
%        plot(R_lin2_3,Z_lin2_3,'b*');
%       plot(R_lin3_4,Z_lin3_4,'go');
%       plot(R_lin4_1,Z_lin4_1,'k*');
      
%Accumualtions of R,Z points and widths:
xaccum=[R_lin1_2 R_lin2_3 R_lin3_4 R_lin4_1]';
yaccum=[Z_lin1_2 Z_lin2_3 Z_lin3_4 Z_lin4_1]';

ww_R=[wwR_1_2 wwR_2_3 wwR_3_4 wwR_4_1]';
ww_Z=[wwZ_1_2 wwZ_2_3 wwZ_3_4 wwZ_4_1]';

    %accum_rep=[xaccum; yaccum];    %Check

%Remove duplicates
dup = (abs(diff(xaccum))+abs(diff(yaccum))) > 0;
xaccum = xaccum(dup);
yaccum = yaccum(dup);

ww_R=ww_R(dup);
ww_Z=ww_Z(dup);

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
    vessel_filament(i) = fiesta_filament(xaccum(i),yaccum(i),ww_R(i),ww_Z(i),1,0,0); 
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

%%%%%%%%%%%%%%%%%%%%COILSET%%%%%%%%%%%%%%

coilset = fiesta_coilset('STcoilset',[Sol_circuit,PF1,PF2,Div1,Div2],false,xaccum',yaccum');


%Plot of the cross section
%     figure;
%     set(gca, 'DataAspectRatio', [1,1,1], 'NextPlot', 'add')
%     c=plot(coilset);
%     set(c, 'EdgeColor', 'k')
%     hold on
%     c=plot(vessel);
%     set(c, 'EdgeColor', 'k')    
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     title('Cross-section')
%     %%%OPtionf for tfg
%     %set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%     %axis([0,1.1,-1.1,1.1]) 
%     %print -depsc2 NOMBREPLOT.eps
%     Filename = '_CrossSection';
%     saveas(gcf, strcat(ProjectName,Filename,FigExt));

% @@@@@@@END CREATION OF THE TOKAMAK@@@@@@@@@@@@@@@


%%%%%%CURRENT PROFILE FOR THE COILSET

%For simplicity, and following Scott code, will declare the input profile
%at the beginnin, since this will be clearer!

%Notes:
    %-Provided that Ip>0, I_PF1,2 have to be<0, so it repel the plasma, and it do not 
        %touches the outer (R max) wall
    
%%%PF coil currents (At Equilibrium, time(5,6))                                                   [A]
I_PF1_Equil=-.8e3;	%-.9 good	%-1.5e3			%First guess, efit will found it. 
I_PF2_Equil=-1e3;	%-1e3 good	%-.5e4			%First guess, efit will found it. 
                                               %These first guesses have to
                                               %be accurate!!!! It affects
                                               %the results!!!!!!!!!
%I_Div1_Equil=NaN;					Div1 in series with Sol!
I_Div2_Equil=+3.3e3;	%3.2e3 good %3.3e3		%4e3 		

%%%Sol current
 I_Sol_max =2.7e3;%2.5e3                %[A] Max solenoid current, to achieve the desire Ip in RZIp.
ISol_equil=-I_Sol_max;%.05e3;                            % [A] the current of the Sol in the target eq calc
%If>0, Sol will attract the plasma, no lower PF and DIv are
%needed.



%%%Time intervals
nTime = 7;      	% Coil Waveform Timesteps	[-]
discharge_time=100; %[ms], time duration of the flat-top 

T_ramp_Sol_loop=[25 10 5 2.5];

tic %to measure time for the loop

for loopM=1:length(T_ramp_Sol_loop) 

T_ramp_Sol=T_ramp_Sol_loop(loopM)              %[ms], min time to increase or decrease the Sol

t_add=2*T_ramp_Sol; %2 with short breakdown               %time for the aditional point>T_ramp_Sol

time = [-100 -50 0 T_ramp_Sol t_add ...
      t_add+discharge_time+10 t_add+discharge_time+50]*1e-3;                   %[s]
 %%%%%% 
 
 %Points 3 and 4 are (0,I_Sol_start), (t,ISol_to). Lets find the equation
 %for that line, and create the point 5 so it has the same slope
 
 [coef]=polyfit([0 T_ramp_Sol],[I_Sol_max 0],1);
 
 %If t5=something, we get I_5 should be:
 I_5=coef(1)*t_add+coef(2); %=I_Sol_max, como debe ser. Por eso, uso I_max.!!
 
 %%%
 
%Construct Sol, PF/Div coil current waveforms vertices
											  %!Breakdown!					%!Efit Icoil!
        %Time   	     [1,    2,              3,           4,        5,                  6,               7];
ISol_Waveform =  [0,  I_Sol_max, I_Sol_max, 0, -I_Sol_max, -I_Sol_max*1,   0];
IPF1_Waveform =  [0,  NaN,        NaN,        NaN,           I_PF1_Equil,   I_PF1_Equil,   0];
IPF2_Waveform =  [0,  NaN,        NaN,        NaN,           I_PF2_Equil,   I_PF2_Equil,   0];
%IDiv1_Waveform = [0,  NaN,        NaN,        NaN,           I_Div1_Equil,  I_Div1_Equil, 0];
IDiv1_Waveform = ISol_Waveform; %Div1 in series with Sol
IDiv2_Waveform = [0,  NaN,        NaN,        NaN,           I_Div2_Equil,  I_Div2_Equil,  0];
%%%%%
CoilWaveforms = [ISol_Waveform; IPF1_Waveform; IPF2_Waveform; IDiv1_Waveform; IDiv2_Waveform];


%Loop voltage calc

slope= (ISol_Waveform(5)-ISol_Waveform(3))/(time(5)-time(3));          %[A/t] slope of the ramp down of the Sol
V_loop=abs(mu0*turns(iSol)/(ZMax_Sol-ZMin_Sol)*pi*RSol^2*slope)     %[V] 1 loop voltage induced by Sol only

%VEST value: Vloop about 3V     GlobusM: about 2V

%Vloop=0.24V for 0.25s T ramp Sol.  Too low
%Vloop=0.6V for 0.10s T ramp Sol. Too low still 
%Vloop=1.2V for 0.05(5ms) T ramp Sol, a bit low yet, but time is too low

E=V_loop/(2*pi*0.45)            %[V/m]


% %% Loop for Paschen plot
% 
%     %2) Paschen curves!
%     C_1=[510 300]; %C_1 constant [ m^-1 Tor^-1] for H and He
%     C_2=[1.25e4 3.4e4];  %C_2 constant [V m^-1 Tor^-1]
%     Gas_type=["H_2","He"];
%     p=linspace(1e-6,1e-3,100000);                           %[Tor]pressure of the prefill gas. size>1000 because if
%                                                                             %not, It do not work properly
%     Emin= @(L,p,C1,C2) C2*p./log(C1*p*L);
% 
%     T_ramp_Sol=[25 10 5 2.5];  
% 
% for tt=1:length(T_ramp_Sol)
%     T_ramp=T_ramp_Sol(tt);
%     t_add=2*T_ramp; %2 with short breakdown               %time for the aditional point>T_ramp_Sol
% 
%     time = [-100 -50 0 T_ramp t_add ...
%       t_add+discharge_time+10 t_add+discharge_time+50]*1e-3;  
%     %Loop voltage calc
% 
%     slope= (ISol_Waveform(5)-ISol_Waveform(3))/(time(5)-time(3));          %[A/t] slope of the ramp down of the Sol
%     V_loop(tt)=mu0*turns(iSol)/(ZMax_Sol-ZMin_Sol)*pi*RSol^2*slope;     %[V] 1 loop voltage induced by Sol only
%     E(tt)=abs(V_loop(tt)/(2*pi*0.45));  
%     
%  
% 
% end
%         figure;
%         %subplot(1,2,1)
%         loglog(p,Emin(10,p,C_1(1),C_2(1)),'*-')
%         hold on
%         loglog(p,Emin(50,p,C_1(1),C_2(1)),'*-')
%         loglog(p,Emin(100,p,C_1(1),C_2(1)),'*-')
%         loglog(p,Emin(500,p,C_1(1),C_2(1)),'*-')
%         %VEST
%         loglog(p,3/(2*pi*0.36)*ones(1,length(p)),'b--','LineWidth', 0.95) 
%         loglog([2e-5 2e-5],[10^-1 10^5],'b--','LineWidth', 0.95)
%         loglog([3e-5 3e-5],[10^-1 10^5],'b-.','LineWidth', 0.95)
%         %GlobusM
%         loglog(p,4.5/(2*pi*0.36)*ones(1,length(p)),'k--','LineWidth', 0.95)
%         loglog(p,8/(2*pi*0.36)*ones(1,length(p)),'k--','LineWidth', 0.95)
%         loglog([3e-5 3e-5],[10^-1 10^5],'k--','LineWidth', 0.95)
%         loglog([6e-5 6e-5],[10^-1 10^5],'k--','LineWidth', 0.95)
%         %
%         loglog(p,abs(E(1))*ones(1,length(p)),'r-','LineWidth', 0.95)           
%         loglog(p,abs(E(2))*ones(1,length(p)),'r-','LineWidth', 1.15)   
%         loglog(p,abs(E(3))*ones(1,length(p)),'r-','LineWidth', 1.35)  
%         loglog(p,abs(E(4))*ones(1,length(p)),'r-','LineWidth', 1.55)  
%         xlabel('Prefill pressure (Torr)')
%         ylabel('E_{min} (V/m)')
%         legend('L=10m','L=50m (Globus)','L=100m','L=500m','','VEST','','','GlobusM','','','T_{ramp}=50ms (original)','T_{ramp}=20ms','T_{ramp}=10ms','T_{ramp}=5ms')
%         title(sprintf('Paschen curve, Gas=%s',Gas_type(1))); %d for numbers
%         set(gca, 'FontSize', 13); %<- Set properties TFG

%%  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@CONFIGURATION OF FIESTA@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%%% Fiesta grid %%%%%%%%%%%%%%%%%%%%
%The grid where FIESTA will run
%Modifying the grid will get different things

R_simulation_limits = [0.03 1];                                                    %R limits of the grid(has to contain vessel)
Z_simulation_limits = [-1.1 1.1];         %Z limits of the grid(has to contain vessel)
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
Ip=100e3;                                    %[A] plasma current
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
icoil_equil.Sol=ISol_equil;         %Solenoid Equilibrium Current at time(5,6)
icoil_equil.PF1=CoilWaveforms(2,6);	%PF1 Equilibrium Current at time(5,6)
icoil_equil.PF2=CoilWaveforms(3,6);	%PF2 Equilibrium Current at time(5,6)
%icoil_init.Div1=CoilWaveforms(4,6);	%Div1 Equilibrium Current at time(5,6)
icoil_equil.Div1=icoil_equil.Sol;         %Div1 in series with Sol
icoil_equil.Div2=CoilWaveforms(5,6);	%Div2 Equilibrium Current at time(5,6)

%To do the equilibria calc, we can use EFIT algorithm or not:

    %1) NO EFIT..........
    
    %equil = fiesta_equilibrium( 'STV2C2', config, Irod, jprofile, control, [],icoil ); 


    %2)EFIT..........
    %this calculates PF and Divs current given plasma parameters
    %Discovered by Juanjo Toledo Garcia
    
    %[efit_config, signals, weights, index]=fiesta_efit_configuration(config, {'PF1','PF2'}, [0.44, 0, 0.44/1.85 1.8 0.2])
    [efit_config, signals, weights, index]=efit_shape_controller(config, {'PF1','PF2'}, [0.44, 0, 0.44/1.85 1.8 0.2]);
    
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

%     figure; hold on; axis equal;
%     plot(coilset);
%     contour( get(equil,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
%     contour( get(equil,'Psi'),get(equil,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
%     plot(vessel);
%     fileName = 'ST_target_equilibrium';
%     legend(gca,'hide');
%     %set(gca,'XLim',[0 1]);
%     %set(gca,'YLim',[-1.5 1.5]);
%     xlabel(gca,'R (m)');
%     ylabel(gca,'Z (m)');
%     title('Target equilibria phase 1');
%     % save_to_pdf( gcf, fileName );
%     %%%OPtionf for tfg
%     set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%     axis([0,1.1,-1.1,1.1]) 
%     %print -depsc2 NOMBREPLOT.eps
    
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
            %Watch out, R0 in the plot is r0_mag, not r0_geom!!
        %Filename = '_TargetEquilibrium';
        %saveas(gcf, strcat(ProjectName,Filename,FigExt));
        
        clc %to delete all teh warnings that appears

        %Equilibrium parameters

parameters(equil); %its better in this way, because it shows units and coil currents. If you define a variable, it wont do that
param_equil=parameters(equil);                             %Will be used in the null sensors


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

% %All the vessel
% BP_virt_R = linspace(VesselRMinPoint,VesselRMaxPoint,10);   %R values, 100 values
% BP_virt_Z = linspace(VesselZMinPoint,VesselZMaxPoint,10);             %Z values, 100 values
% 
% a_eff=VesselRMaxPoint-VesselRMinPoint
%%%%%%%%%%5

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

global R_sensor Z_sensor
[R_sensor,Z_sensor]=meshgrid(r_sensors,z_sensors); %size 200*200

%Plot of the sensors
%     figure; hold on; axis equal;
%     plot(coilset);
%     contour( get(equil,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
%     contour( get(equil,'Psi'),get(equil,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
%     plot(vessel);
%     fileName = 'ST_target_equilibrium';
%     legend(gca,'hide');
%     plot(sensor_btheta);
%     %set(gca,'XLim',[0 1]);
%     %set(gca,'YLim',[-1.5 1.5]);
%     xlabel(gca,'R (m)');
%     ylabel(gca,'Z (m)');
%     title('Target equilibria phase 1 with sensors to null Bpol');
%     % save_to_pdf( gcf, fileName );
%     %%%OPtionf for tfg
%     set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%     axis([0,1.1,-1.1,1.1]) 
     
%         figure;
%         plot(equil)        
%         hold on
%         plot(vessel)
%         plot(coilset)
%         parametersshow(equil)   %this plots the parameters in the equil
%         title('Target equilibria')
%         plot(sensor_btheta);
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
%     figure;
%     plot( time*1e3, I_PF_input/(1e3),'*-' );
%     xlabel('time (ms)')
%     ylabel('I_{{input}} (kA)')
%     title('I_{{input}} versus time')
%     legend('Sol','PF1','PF2','Div1','Div2')
%     %%%optinos for tfg
%     set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%     set(gca, 'FontSize', 13, 'LineWidth', 0.75);                    %<- Set properties TFG
    %Filename = '_InputCurrents';
    %saveas(gcf, strcat(ProjectName,Filename,FigExt));



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


%% RZIP AS FIESTA%%%%%%EXPERIMENTAL%%%%

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
%             figure;
%             subplot(3,1,1)
%             %plot(time_adaptive*1e3,I_PF_output/(1e3))
%             %hold on
%             plot(time_adaptive*1e3,Ip_output/(1e3))
%             ylabel('I (kA)')
%             title('Dynamic response SMART phase 1')
%             set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%             %set(gca,'YLim',[-5 35]);
%             set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%         %%%vOLTAGE
%             subplot(3,1,2)
%             plot(time_adaptive*1e3,V_PF_output/(1e3))
%             hold on
%             plot(time_adaptive*1e3,Vp_output/(1e3))
%             ylabel('V (kV)')
%             legend('Sol','PF2','PF3','Div1','Div2','Plasma')
%             set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%             %set(gca,'YLim',[-1.500 1.500]);
%             set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%         %%%I_PASSIVE VERUS TIME
%             subplot(3,1,3)
%             plot(time_adaptive*1e3,I_Passive_VV/(1e3))
%             xlabel(' time (ms)')
%             ylabel('I_{passive} (kA)')
%             set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%             %set(gca,'YLim',[-800 800]);
%             set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%             
%             Filename = '_RZIpOutputs';
%             saveas(gcf, strcat(ProjectName,Filename,FigExt));
%           %print -depsc2 NOMBREPLOT.eps
% 
%             figure;           
%             plot(time_adaptive*1e3,Ip_output/(1e3))
%             ylabel('I_p (kA)')
%             xlabel(' time (ms)')
%             title('I_p SMART phase 1')
%             set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%             %set(gca,'YLim',[-5 35]);
%             set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
% 
%             figure;
%             plot(time_adaptive*1e3,I_Passive_VV/(1e3))
%             xlabel(' time (ms)')
%             ylabel('I_{passive} (kA)')
%             title('I_{VV}')
%             set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%             set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
            %Experimental plots to see Vp in detail
            
            %save_10ms=[time_adaptive Ip_output I_Passive_VV]
            %save('10ms_time_Ip_IVV','save_10ms')
        
            %% PLOT THINGS FOR SLIDES JUNE
%             %Ip plot
%             figure;
%             plot(save_5ms(:,1)*1e3,save_5ms(:,2)*1e-3, 'LineWidth', 0.95)
%             hold on
%             plot(save_10ms(:,1)*1e3,save_10ms(:,2)*1e-3, 'LineWidth', 0.95)
%             plot(save_25ms(:,1)*1e3,save_25ms(:,2)*1e-3, 'LineWidth', 0.95)
%             plot(Save_50ms(:,1)*1e3,Save_50ms(:,2)*1e-3, 'LineWidth', 0.95)
%             xlabel('t (ms)')
%             ylabel('I_p (kA)')
%             title('I_p as a function of the ramp down time ')
%             legend('5ms','10ms','20ms','50ms (original)')
%             set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%             
%             %IVV
%                         figure;
%             plot(save_5ms(:,1)*1e3,save_5ms(:,3)*1e-3, 'LineWidth', 0.95)
%             hold on
%             plot(save_10ms(:,1)*1e3,save_10ms(:,3)*1e-3, 'LineWidth', 0.95)
%             plot(save_25ms(:,1)*1e3,save_25ms(:,3)*1e-3, 'LineWidth', 0.95)
%             plot(Save_50ms(:,1)*1e3,Save_50ms(:,3)*1e-3, 'LineWidth', 0.95)
%             xlabel('t (ms)')
%             ylabel('I_{VV} (kA)')
%             title('I_{VV} as a function of the ramp down time ')
%             legend('5ms','10ms','20ms','50ms (original)')
%             set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
            %%       
%             
%             figure;
%             subplot(1,4,1)
%             plot(time_adaptive*1e3,I_PF_output(:,1)/(1e3))
%             xlabel(' t (ms)')            
%             ylabel(' I CS (kA)')
%             title('I Sol')
%             subplot(1,4,2)
%             plot(time_adaptive*1e3,Ip_output/(1e3))
%             ylabel('I_p (kA)')
%             xlabel(' t (ms)')
%             title('I_p phase 1')
%             subplot(1,4,3)
%             plot(time_adaptive*1e3,Vp_output)
%             xlabel('t (ms)')
%             ylabel('Vp (V)')            
%             title('Vp')
%             subplot(1,4,4)
%             plot(time_adaptive*1e3,V_PF_output(:,1))
%             xlabel('t (ms)')
%             ylabel('V Sol (V)')
%             title('V CS')
            
            %DO NOT UNDERSTAND WHY THE F Vp IS NOT Vloop AT LEAST AT THE 
            %BEGINNING!!!!!!!!!!!!!
%%% END RZIP@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


%     %Inside uFinal(x)
%     figure;
%     subplot(2,4,1)
%     plot(time_adaptive*1e3,uFinal(:,end)/(1e3))  
%     xlabel('t(ms)')
%     ylabel('Ip (kA)')
%     title('uFinal(:,end)')
%         subplot(2,4,2)
%     plot(time_adaptive*1e3,uFinal(:,[end-5:end-1]))  
%     xlabel('t(ms)')
%     ylabel('¿?')
%     title('uFinal(:,[end-5:end-1]) ¿?')
%             subplot(2,4,3)
%     plot(time_adaptive*1e3,uFinal(:,[1:end-6])/(1e3))  
%     xlabel('t(ms)')
%     ylabel('I_VV (kA)')
%     title('uFinal(:,[1:end-6])')
%                 subplot(2,4,4)
%     plot(time_adaptive*1e3,sum(uFinal(:,[1:end-6])/(1e3),2))  
%     xlabel('t(ms)')
%     ylabel('I_VV(total) (kA)')
%     title('sum(uFinal(:,[1:end-6]),2)')
%             subplot(2,4,5)
%     plot(time_adaptive*1e3,uFinal(:,[end-5:end-1])./uFinal(:,end))  
%     xlabel('t(ms)')
%     ylabel('I PF (kA)')
%     title('uFinal(:,[end-5:end-1])/uFinal(:,end)')
%             subplot(2,4,7)
%                 plot(time_adaptive*1e3,V_PF_output/(1e3))
%             hold on
%             plot(time_adaptive*1e3,Vp_output/(1e3))
%             ylabel('V (kV)')
%             legend('Sol','PF2','PF3','Div1','Div2','Plasma')
%             title('Volages')
%             %set(gca,'YLim',[-1.500 1.500]);
%             subplot(2,4,8)
%             plot(time_adaptive*1e3,I_PF_output/(1e3))
%             title('I_PF (kA)')
%             
%             %Weirds elements of uFinal in detail
%             figure;
%             subplot(2,5,1)
%             plot(time_adaptive*1e3,uFinal(:,end-5))  
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-5)')
%             subplot(2,5,2)
%             plot(time_adaptive*1e3,uFinal(:,end-4))  
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-4)')
%             subplot(2,5,3)
%             plot(time_adaptive*1e3,uFinal(:,end-3))  
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-3)')
%             subplot(2,5,4)
%             plot(time_adaptive*1e3,uFinal(:,end-2))  
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-2)')           
%             subplot(2,5,5)
%             plot(time_adaptive*1e3,uFinal(:,end-1))  
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-1)')             
%             subplot(2,5,6)
%             plot(time_adaptive*1e3,uFinal(:,end-5)./uFinal(:,end))  
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-5)/uFinal(:,end)')    
%             subplot(2,5,7)
%             plot(time_adaptive*1e3,uFinal(:,end-4)./uFinal(:,end))    
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-4)/uFinal(:,end)')    
%             subplot(2,5,8)
%             plot(time_adaptive*1e3,uFinal(:,end-3)./uFinal(:,end))    
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-3)/uFinal(:,end)')    
%             subplot(2,5,9)
%             plot(time_adaptive*1e3,uFinal(:,end-2)./uFinal(:,end))    
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-2)/uFinal(:,end)')              
%             subplot(2,5,10)
%             plot(time_adaptive*1e3,uFinal(:,end-1)./uFinal(:,end))    
%             xlabel('t(ms)')
%            ylabel('¿?')
%            title('uFinal(:,end-1)/uFinal(:,end)')    
                
    %This mean that the last column of uFinal is Ip, and the previous one
    %are coils current and structure's current, so no RIp,ZIp :((
    
    
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

r_inside_VV=linspace(VesselRMinPoint,VesselRMaxPoint,n_pnts_inside); 
z_inside_VV=linspace(VesselZMinPoint,VesselZMaxPoint,n_pnts_inside);

%Have to do ameshgrid for the interpolation, so each R value has a Z value
global R_in Z_in
[R_in,Z_in]=meshgrid(r_inside_VV,z_inside_VV);

%Some fiesta things are needed for the equil calc:
coilsetVV=fiesta_loadassembly(coilset, vessel);
configVV = fiesta_configuration( 'SMART with VV', Grid, coilsetVV);
  

%%%%Breakdown non optimised(only Sol)%%%%%%%%
 coil_currents_Break_non = zeros(1,nPF);
 coil_currents_Break_non(iSol) = I_Sol_max; %current at the top previous to ramp down
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%So, what I will do is a loop with the time from RZIp, to define the
%current values. The time steps will be
time_loop=[0]*1e-3 %[s]
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
    end
    
    %With the currents, the icoil object can be created, and with it the
    %equil
    icoilVV_break = fiesta_icoil( coilsetVV, coil_currents_break );
    equil_break = fiesta_equilibrium( 'optimised null (eddys in)', configVV, Irod, icoilVV_break );

%     %Plot
%     figure;
%     plot(equil_break);
%     hold on;
%     plot(coilset)
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
%         xlabel('R (m)')
%     fileName = '\psi contour optimised_null';
%     legend(gca,'hide');
%     title('\psi contours at t=0ms phase 1')
%     set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
 
%Manual plot
psi_null=get(equil_break,'Psi_vac'); %fiesta field


%
psi_null=get(psi_null,'data','2D'); %this retunrs 2D data, so no need to reshape!
                    %251*200

%  % PLOT WHOLE VV
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
%         
psi_null_interpn=  @(r,z) interpn(zGrid,rGrid,psi_null,z,r,'mikama');        
psi_null_ins_VV=psi_null_interpn(R_in,Z_in);    
    
%     figure;
%     contour(R_in,Z_in,log10(psi_null_ins_VV),1000)
%     hold on
%     plot(vessel)
%     plot(coilset)
%     colormap(Gamma_II)
%     c=colorbar; %colorbar
%         ylabel(c, 'log10(\psi)');
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%       [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
%      xlabel('R (m)')
%      ylabel('Z (m)')
%      title(sprintf('log10(psi)  at t=%d ms (iter %d/%d) 1000c',time_loop(loop)*1e3,loop,length(time_loop)))
%           
 %Extraction of the fields (Earth inside)!
 [FieldsBreak, FieldsBreakNoEarth]=fields(equil_break);
    
    %Plots!
    % % %%Quiver plot
%  scale_factor=1; %graphic needs to be scaled
% figure; 
% %subplot(1,3,1)
% quiver(R_in,Z_in,Br_ins_vessel/scale_factor, Bz_ins_vessel/scale_factor,'color',[1 0 0],'AutoScale','on','AutoScaleFactor', 10)
% hold on;
% plot(vessel)
%    view(2) %2D view
%       plot(sensor_btheta)
% xlabel('R (m)')
% ylabel('Z (m)')
% title('\vec{B} inside vessel')
% subplot(1,3,2)
% quiver(R_in,Z_in,BrEarthVV/scale_factor, BzEarthVV/scale_factor,'color',[1 0 0],'AutoScale','off')
% hold on;
% plot(vessel)
%    view(2) %2D view
%       plot(sensor_btheta)
% xlabel('R (m)')
% ylabel('Z (m)')
% title('\vec{B_Earth} inside vessel')
% subplot(1,3,3)
% quiver(R_in,Z_in,Br_ins_vesselE/scale_factor, Bz_ins_vesselE/scale_factor,'color',[1 0 0],'AutoScale','off')
% hold on;
% plot(vessel)
%    view(2) %2D view
%       plot(sensor_btheta)
% xlabel('R (m)')
% ylabel('Z (m)')
% title('\vec{B_tokamak}+\vec{B_Earth} inside vessel')


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
%     title('log10(Br) inside vessel ')
%        
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
%     title('log10(Bz) inside vessel')
%   
% 

    %Plot Bpol y quiver
%      figure; 
%     contourf(R_in,Z_in,log10(abs(Bpol_ins_vessel)));
%     shading('interp') %this is to make the transition between values continuous,
%     %instedad of discontinuously between pixels
%     hold on;
%     plot(vessel)
%     colorbar %colorbar
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%      [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%     quiver(R_in,Z_in,Br_ins_vessel, Bz_ins_vessel,'color',[1 0 0],'AutoScale','off')
%      xlabel('R (m)')
%     ylabel('Z (m)')
%     title(sprintf('log10(Bpol) and Bpol quiver at t %d ms (iter %d/%d (P2)',time_loop(loop)*1e3,loop,length(time_loop)))


%     %Field of the plasma with the circular wire approx
%     
%     for k=1:length(R_in)
%     
%     for l=1:length(Z_in)
%         B_p_pol_wire(k,l)=Campo_mg_espira_Cilind(R_in(k,l),0,Z_in(k,l),param_equil.r0_geom,Ip_loop(loop),10)*[1,0,0]'+...
%             Campo_mg_espira_Cilind(R_in(k,l),0,Z_in(k,l),param_equil.r0_geom,Ip_loop(loop),10)*[0,0,1]';  %Pol field of the wire model
%                                                                                                     %of the plasma
%         B_p_phi_wire(k,l)=Campo_mg_espira_Cilind(R_in(k,l),0,Z_in(k,l),param_equil.r0_geom,Ip_loop(loop),10)*[0,1,0]'; %Toroidal field
%         
%         %Check
%         B_p_r_wire(k,l)=Campo_mg_espira_Cilind(R_in(k,l),0,Z_in(k,l),param_equil.r0_geom,Ip_loop(loop),10)*[1,0,0]'; %Br
%         B_p_z_wire(k,l)=Campo_mg_espira_Cilind(R_in(k,l),0,Z_in(k,l),param_equil.r0_geom,Ip_loop(loop),10)*[0,0,1]'; %Bz
%     end
%     end

    %Plots plasma field as a circular wirw
       
%         subplot(1,2,2)
%         contourf(R_in,Z_in,log10(abs(B_p_pol_wire)));
%         shading('interp') %this is to make the transition between values continuous,
%         %instedad of discontinuously between pixels
%         colormap(Gamma_II)
%         hold on;
%         plot(vessel)
%         colorbar %colorbar
%         view(2) %2D view
%         plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%         [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
%         xlabel('R (m)')
%         ylabel('Z (m)')
%         title(sprintf('log10(B_pol plasma)  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))
%         
        %Bphi structure
%         figure;
%         contourf(R_in,Z_in,log10(FieldsBreak.VV.Bphi));
%         shading('interp') %this is to make the transition between values continuous,
%         %instedad of discontinuously between pixels  
%         colormap(Gamma_II)
%         hold on;
%         plot(vessel)
%         colorbar %colorbar
%         view(2) %2D view
%         plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%         [  min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%         xlabel('R (m)')
%         ylabel('Z (m)')
%         title(sprintf('log10(B_phi [T]) at t %d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))
%     
        %B pol alone
         figure; 
        %contourf(R_in,Z_in,log10(abs(Bpol_ins_vessel)),'EdgeColor','none');
        %surf(R_in,Z_in,log10(abs(Bpol_ins_vessel)),'EdgeColor','none');
        contour(R_in,Z_in,log10(FieldsBreak.VV.Bpol),100);
        shading('interp') %this is to make the transition between values continuous,
        %instedad of discontinuously between pixels
        hold on;
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        %colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'log10(Bpol[T])');
        colorbar %colorbar
        view(2) %2D view
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
        [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        xlabel('R (m)')
        ylabel('Z (m)')
        %title(sprintf('B_{pol}  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))
        title(sprintf('B_{pol} at t=%dms for %dms',time_loop(loop)*1e3,T_ramp_Sol))

    
%         figure; 
%         subplot(2,2,1)
%         contourf(R_in,Z_in,log10(abs(B_p_r_wire)));
%         shading('interp') %this is to make the transition between values continuous,
%        %instedad of discontinuously between pixels
%         hold on;
%         plot(vessel)
%         colorbar %colorbar
%         view(2) %2D view
%         plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%         [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%         xlabel('R (m)')
%         ylabel('Z (m)')
%         title(sprintf('log10(Br plasma (wire))  at t %d ms (iter %d/%d (P2)',time_loop(loop)*1e3,loop,length(time_loop)))
%         %
%         subplot(2,2,2) 
%         contourf(R_in,Z_in,log10(abs(B_p_z_wire)));
%         shading('interp') %this is to make the transition between values continuous,
%         %instedad of discontinuously between pixels
%         hold on;
%         plot(vessel)
%         colorbar %colorbar
%         view(2) %2D view
%         plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%         [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%         xlabel('R (m)')
%         ylabel('Z (m)')
%         title(sprintf('log10(Bz plasma (wire))  at t %d ms (iter %d/%d (P2)',time_loop(loop)*1e3,loop,length(time_loop)))
%         %
%         subplot(2,2,3)
%         contourf(R_in,Z_in,log10(abs(B_p_pol_wire)));
%         shading('interp') %this is to make the transition between values continuous,
%         %instedad of discontinuously between pixels
%         hold on;
%         plot(vessel)
%         colorbar %colorbar
%         view(2) %2D view
%         plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%         [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%         xlabel('R (m)')
%         ylabel('Z (m)')
%         title(sprintf('log10(Bpol plasma (wire))  at t %d ms (iter %d/%d (P2)',time_loop(loop)*1e3,loop,length(time_loop)))
%         %
%         subplot(2,2,4)
%         contourf(R_in,Z_in,log10(abs(B_p_phi_wire)));
%         shading('interp') %this is to make the transition between values continuous,
%         %instedad of discontinuously between pixels
%         hold on;
%         plot(vessel)
%         colorbar %colorbar
%         view(2) %2D view
%         plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%         [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'r.--')
%         xlabel('R (m)')
%         ylabel('Z (m)')
%         title(sprintf('log10(Bphi plasma (wire))  at t %d ms (iter %d/%d (P2)',time_loop(loop)*1e3,loop,length(time_loop)))

   %%L CALC FORMULAE%%%%%%%%%%%%

    %The field in the whole vessel will be used, since it is the field that is
    %used in the field line integrator.

    %Minimun with veesel interp 
    Bpolmin=min(min(FieldsBreak.sensor.Bpol));                                         %minimun of poloidal field
    [Bpolmin_index_row Bpolmin_index_column]=find(FieldsBreak.sensor.Bpol==Bpolmin);     %indexes  

        %%%%%THIS IS TO COMPARE
    %     %Minimun with sensors interp
    % Bpolmin=min(min(Bpol_sensor)) %minimun of poloidal field
    % [Bpolmin_index_row Bpolmin_index_column]=find(Bpol_sensor==Bpolmin);
    %     %indexes
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
    L_aver(loop)=0.25*a_eff*Bphi_centerNull/Bpolmin_av                %[m] L with the average pol field

    %REMEMBER THAT THIS IS COMPUTED INSIDE NULL REGION, WHILE THE LOWER
    %BPOL IS NOT IN THAT REGION DUE TO EDDYS (AND EARTH)!!!!!!!!!!!!!!!!!!
    %%NEED TO THINK ON THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

    %%L calc by field line integration###########################

           %Int eh Lazarus paper 1998, they compute connective length by avergaing on
        %9 lines, the line with Bpol min, and the 8 surroundings. I will do the
        %same. No need for loop, to avoid problems. I'll do it manually xD.
        %However, here I will do the calc inside all the VV
      

        %Grid
        %I redefine the grid to compute the connection length, for less computer
        %demands (time)

        n_pnts_insideL=30 %15            %100 is the ideal to have good plots of the fields, but the L int failures. 

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
        r_inside_VV_noLimits=r_inside_VVL; %r_inside_VVL(2:end-1);
        z_inside_VV_noLimits=r_inside_VVL; %r_inside_VVL(2:end-1);
                    %   NOW IT CONTAINS ALL THE POINTS(REFINED GRID)
        %Mesh to future plots
        [r_insVV_noLimit,z_insVV_noLimit]=meshgrid(r_inside_VV_noLimits,z_inside_VV_noLimits);


        r0_z0_L0_Phi0_U0=[0 0 0 0 0]; 
            %note it has to be r z fro using the same event function

        for i=1:length(r_inside_VV_noLimits)

            for j=1:length(z_inside_VV_noLimits)
        
                points=[r_inside_VV_noLimits(i) z_inside_VV_noLimits(j) 0];  %r z phi         
                r0_z0_L0_Phi0_U0=[ r0_z0_L0_Phi0_U0; points 0 0];    %U(0)=0 (arbitrary)    
            end
    
        end

        %I have the additional point 0 0 0 form the begining, that i can remove
        %easily with
        r0_z0_L0_Phi0_U0=r0_z0_L0_Phi0_U0(2:end,:);
    
    %%%%%%%            
            
        Lp_values=1000; 
        n_iterLp=3000;         %1000                                %integer, Number of iterations!!!!
        LpSpan=linspace(0,10*n_iterLp,Lp_values*n_iterLp);    %the range of values of phi
        L_max=10000;                                                 %max value for the integration; when L achieves
                                                            %this value, the integration stops. ST have around 50m.

        odefun= @(Lp, rzLphiU) Field_LineIntegrator_Lp(Lp,rzLphiU,FieldsBreak.interpn.Br,...
            FieldsBreak.interpn.Bz,FieldsBreak.interpn.Bphi);
        event_colission_wall=@(Lp,rzLphiU) Colission_wall(Lp,rzLphiU,VesselRMaxPoint,...
            VesselRMinPoint,VesselZMaxPoint,VesselZMinPoint,FieldsBreak.interpn.Br,...
            FieldsBreak.interpn.Bz,FieldsBreak.interpn.Bphi,L_max); 
        options = odeset('OutputFcn',@ode_progress_bar,'Events',event_colission_wall,'AbsTol',1e-10,'RelTol',1e-6); 
                                    %I include a fiesta funciton to show the progress of the ode

        tic              %to know the time
%         %%%%%%%%SINGLE FIELD LINE TRACER
% 
%         %Single integrator and plotter of lines
%             %need to find i for the chosen R,Z value in r0_z0_L0_U0.
%             %I= 85 for a line inside, 49 for a max L outside, 152 for the
%             %outward arm (Z>0). 135 for the outward Z<0 line. 64 for the upper
%             %arm
%         
%             i=1 %looked in the vector
%             [Lp_fieldline, rzLphiU_fieldline]=ode45(odefun,LpSpan,r0_z0_L0_Phi0_U0(i,:),options);        %ode15s Carlos
%     
%             %To save the last values of R,Z,L
%             RZLPhiU_end(i,1)=rzLphiU_fieldline(end,1);
%             RZLPhiU_end(i,2)=rzLphiU_fieldline(end,2);
%             RZLPhiU_end(i,3)=rzLphiU_fieldline(end,3);                          
%             RZLPhiU_end(i,4)=rzLphiU_fieldline(end,4);   
%             RZLPhiU_end(i,5)=rzLphiU_fieldline(end,5);
%         % %%Plot of one of the line
%         % 
%         figure;
%         plot3(r0_z0_L0_Phi0_U0(i,1),r0_z0_L0_Phi0_U0(i,3),r0_z0_L0_Phi0_U0(i,2),'k*','LineWidth',3)
%         hold on;
%         plot3(vessel)
%         plot3(coilset)
%         plot3(rzLphiU_fieldline(:,1).*cos(rzLphiU_fieldline(:,3)),rzLphiU_fieldline(:,1).*sin(rzLphiU_fieldline(:,3)),...
%             rzLphiU_fieldline(:,2),'r','LineWidth',3)
%         xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');  
%         hold on
%         plot3(rzLphiU_fieldline(end,1).*cos(rzLphiU_fieldline(end,3)),rzLphiU_fieldline(end,1).*sin(rzLphiU_fieldline(end,3)),...
%             rzLphiU_fieldline(end,2),'g*','LineWidth',3)
%         title('Field line integration Lp (single)')
%         set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%         %legend('Starting point (Point with less Bpol)','Field line',...

        %%%%%%%%END One line tracer%%%%%%%%%%%%%5
    
    %So, gotta solve 10^2*10^2=10^4 eq xD. Since I am only interested in L,
    %could do a for loop simply saving the L value. Will do that:
    
        for i=1:length(r0_z0_L0_Phi0_U0)
            fprintf('Iter %d de %d',i,length(r0_z0_L0_Phi0_U0))
            [Lp_fieldline, rzLphiU_fieldline]=ode45(odefun,LpSpan,r0_z0_L0_Phi0_U0(i,:),options);        %ode15s Carlos
    
            %To save the last values of R,Z,L
            RZLPhiU_end(i,1)=rzLphiU_fieldline(end,1);
            RZLPhiU_end(i,2)=rzLphiU_fieldline(end,2);
            RZLPhiU_end(i,3)=rzLphiU_fieldline(end,3);                          
            RZLPhiU_end(i,4)=rzLphiU_fieldline(end,4);   
            RZLPhiU_end(i,5)=rzLphiU_fieldline(end,5);
        end
       time_int_Lp=toc           %time of the ode
        %To store start points that do not collide: first I get the index of both R
        %and Z, but together, since they do not collide if oth R and Z are greater
        %than the min value, and lower than the greatest value
    
            RZ_store_indexLp=RZLPhiU_end(:,1)<VesselRMaxPoint & ...
                RZLPhiU_end(:,1)>VesselRMinPoint & RZLPhiU_end(:,2)<VesselZMaxPoint &...
                RZLPhiU_end(:,2)>VesselZMinPoint; %100*5==> error, has to ve vector,, not matrix!!!
                        
            RZ_no_collideLp=[r0_z0_L0_Phi0_U0(RZ_store_indexLp,1) r0_z0_L0_Phi0_U0(RZ_store_indexLp,2)];    
       
        %WHEN INCLUDING ALL THE POINTS, WEIRD THINGS HAPPENS WITH THE
        %POINTS AT THE BORDER OF THE GRID, SO THAT THE NON COLLIDE VECTOR
        %COLLIDE SOME OF THOSE POINTS
       
        %Plot contour
        L_int=reshape(RZLPhiU_end(:,3),size(r_insVV_noLimit,1),size(r_insVV_noLimit,2));
        
        figure;
        contourf(r_insVV_noLimit,z_insVV_noLimit,L_int,10)
        %surf(r_insVV_noLimit,z_insVV_noLimit,L_int,'EdgeColor','none'), shading('interp')
        hold on
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'g.--')
        plot(RZ_no_collideLp(:,1),RZ_no_collideLp(:,2),'m*')
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'L(m)');
        xlabel('R (m)')
        ylabel('Z (m)')
        legend('L','Field null region')
        %title(sprintf('L  at t=%d ms (iter %d/%d) LP',time_loop(loop)*1e3,loop,length(time_loop)))   
        title(sprintf('L at t=%dms for %dms',time_loop(loop)*1e3,T_ramp_Sol))
    
        %Plot contour 'potential'
        U_int=reshape(RZLPhiU_end(:,5),size(r_insVV_noLimit,1),size(r_insVV_noLimit,2));
        figure;
        contourf(r_insVV_noLimit,z_insVV_noLimit,U_int,10)
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
        xlabel('R (m)')
        ylabel('Z (m)')
        %title(sprintf('Pseudo potential  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        title(sprintf('Pseudo-potential at t=%dms for %dms',time_loop(loop)*1e3,T_ramp_Sol))
        %Plot Lloyd criteria
        Lloyd=E.*FieldsBreak.VV.Bphi./FieldsBreak.VV.Bpol;
        
        figure;
        contour(R_in,Z_in,Lloyd,50,'ShowText','on')
        %surf(r_insVV_noLimit,z_insVV_noLimit,U_int), shading('interp')
        hold on
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        %colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'E*Bphi/Bpol (V/m)');
        view(2) %2D view
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
        xlabel('R (m)')
        ylabel('Z (m)')
        %title(sprintf('Lloyd criteria  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        title(sprintf('Lloyd criteria at t=%dms for %dms',time_loop(loop)*1e3,T_ramp_Sol))
end    %end loop time (not in use right now)!!!!

        L_emp(loopM)=L_aver(loop);              %[m]to store Lempirical!!!
        I_PassiveVV(loopM)=IPassive_loop; %[kA]to store total eddy at t=0ms
end  %end loop time ramp sol

time_loopRamp=toc %time loop for ramp time




% %%
% figure;
% subplot(1,2,2)
% plot(time_loop*1e3,L_aver,'r*-','MarkerSize',5)
% xlabel('time (ms)')
% ylabel('L_{emp} (m)')
% title('L (empirical) versus time')
% set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
% subplot(1,2,1)
% plot(time_loop*1e3,Ip_loop*1e-3,'b*-','MarkerSize',5)
% xlabel('time (ms)')
% ylabel('I_{p} (kA)')
% title('I_p versus time')
% set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
% 
% 
% time_plot=(0:1:5)*1e-3; %s
% I_plot=interpn(time_adaptive,Ip_output,time_plot);
% 
% figure;
% plot(time_plot*1e3,I_plot*1e-3)
% xlabel('time (ms)')
% ylabel('I_{p} (kA)')
% title('I_p versus time')
% set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
% 
% time_plot=(0:1:20)*1e-3; %s
% I_plot=interpn(time_adaptive,Ip_output,time_plot);
% 
% figure;
% plot(time_plot*1e3,I_plot*1e-3)
% xlabel('time (ms)')
% ylabel('I_{p} (kA)')
% title('I_p versus time')
% set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG

%%%%%%TEST FIELD LINE METHOD

	figure;
    contour(R_in,Z_in,log10(sqrt(1+FieldsBreak.VV.Bphi.^2./FieldsBreak.VV.Bpol.^2)),500);
    shading('interp') %this is to make the transition between values continuous,
    %instedad of discontinuously between pixels
    %colormap(Gamma_II)
    hold on;
    plot(vessel)
    c=colorbar; %colorbar
      ylabel(c, 'ds/dLp');
    plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
     [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
     xlabel('R (m)')
    ylabel('Z (m)')
     title(sprintf('log10(ds/dLp)  500c at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))
              
     
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
     
 %Lloyd criteria plot
 
 
%  figure;
%  contour(R_in,Z_in,1./(2*pi*R_in).*FieldsBreak.VV.Bphi./FieldsBreak.VV.Bpol,50,'ShowText','on')
%  %surf(R_in,Z_in,1./(2*pi*R_in).*FieldsBreak.VV.Bphi./FieldsBreak.VV.Bpol), shading('interp')
%  colormap(Gamma_II)
%     hold on;
%     plot(vessel)
%     plot(coilset)
%     c=colorbar; %colorbar
%       ylabel(c, '(E_{phi}B_{phi}/B_{pol})/V_{loop} (1/m)');
%     view(2) %2D view
%     plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
%      [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'k.--')
%      xlabel('R (m)')
%     ylabel('Z (m)')
%      title(sprintf('(E_phiB_phi/B_pol)/V_loop (iter %d/%d)',loop,length(time_loop)))
%      axis([0.1 0.9 0 0.9])       
        

%% ADDITIONAL FEATURES (EDDYS, FORCES, INTRODUCING EDDYS ON EQUIL (TRY), ETC


% %% RE-DOING EQUILIBRIA CALC WITH EDDYS
% 
%     %Now that we havr compute the eddys, we could re do all the calc, to
%     %get an equil with those eddys, with the new equil do the RZip again,
%     %and, it the new eddys do not change much, could accept the result.
%     %Yes, this is a iterative process, but with a huge consumption of
%     %computer resources==> :)
%      
%      %coil currents for the new equilibria calc
%      
%     coil_currents_eddys = zeros(1,nPF+get(vessel,'n'));                       %n=308, number of filaments
%     coil_currents_eddys(iSol) =	CoilWaveforms(iSol,5);             %Do not change, in principle
%     coil_currents_eddys(iPF1) =CoilWaveforms(iPF1,5);  
%     coil_currents_eddys(iPF2) =CoilWaveforms(iPF2,5);  
%     coil_currents_eddys(iDiv1) =CoilWaveforms(iDiv1,5);  
%     coil_currents_eddys(iDiv2) =CoilWaveforms(iDiv2,5);            %May need to be changed
% 
%     %To select the eddy currents, provided that at the equil time the eddys are amx, can
%     %just look at the max:
%     
%     [max_eddy,index_max_eddy]=max(I_Passive_VV);                %the total eddy, the max value
%     coil_currents_eddys(nPF+1:end)=I_Passive(index_max_eddy,:);       %eddy currents, have to define each filament
%     
%     icoil_eddy = fiesta_icoil( coilsetVV, coil_currents_eddys );
%     
%     %%%1) EFIT
%     [efit_configVV, signalsVV, weightsVV, indexVV]=efit_shape_controller(configVV, {'PF1','PF2'}, [0.44, 0, 0.44/1.85 1.8 0.2])
%     % The numbers you give are [Rgeo, Zgeo, a, kappa, delta], Rgeo,Zgeo,a are
%     % mandatory.
%     %I use the values of the standar shape, to get a similar equil
% 
%     equil_eddy=fiesta_equilibrium('Target+eddys', configVV, Irod, jprofile, control,efit_configVV, icoil_eddy, signalsVV, weightsVV) %%EFIT!!!
%     %It does the case in line 96!! The equil calc is in lin 124
% 
%     %Now we have to extract the new currents from the equil, provided that EFIT
%     %changed some of them to satisfy the conditions requested:
%     icoil_eddy=get(equil_eddy,'icoil');
%     current_post_EFIT=get(icoil_eddy,'currents');
%     coil_currents_eddys(iPF1) =current_post_EFIT(iPF1);
%     coil_currents_eddys(iPF2) =current_post_EFIT(iPF2);
%     coil_currents_eddys(iDiv1) =current_post_EFIT(iDiv1);
%     coil_currents_eddys(iDiv2) =current_post_EFIT(iDiv2);
% %No need of redefine the Sol current of course. Actually, Div 1 and 2 are
% %not neccessary, since I am not changing them
% 
% 
% %     %%2) NO EFIT
% %             equil = fiesta_equilibrium( 'STV2C2', config, Irod, jprofile, control, [],icoil );
% 
%                 %This is not an option, since it does not
%                 %convnerge with the values without eddys!!!!!!!!!!
%                              
%   %Plot of the equil with the eddys!
%         %section_figure=section(equil); %THIS IS A PLOT
%         figure;
%         plot(equil_eddy)
%         parametersshow(equil_eddy) %this plots the parameters in the equil
%         hold on
%         plot(vessel)
%         plot(coilset)
%         title('Target equilibria eddys included')


%%% Make virtual sensors where we want breakdown  %%%%%%%
%This is to null the poloidal field(BP), to increase the connective length, and
%allow the plasma breakdown. 

% BP_virt_R = linspace(param_equil.r0_geom-a_eff,param_equil.r0_geom+a_eff,10);       %R values, 100 values
% BP_virt_Z = linspace(ZGeo-a_eff,ZGeo+a_eff,10);                     %Z values, 100 values
% 
% [BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
% BP_virt_R = BP_virt_R(:)';
% BP_virt_Z = BP_virt_Z(:)';
% 
% BP_virt_theta = zeros(1,length(BP_virt_R));
% nSensors = length(BP_virt_theta); %100
% 
% BP_virt_names = {};
% for iSensor=1:nSensors
%     BP_virt_names{iSensor} = ['Radial Bp Virtual Sensor #' num2str(iSensor) ];
% end
% 
% BP_virt_R = [BP_virt_R  BP_virt_R];
% BP_virt_Z = [BP_virt_Z  BP_virt_Z];         %Both size 1*200. It is replicated, so element 101=element 1 
% 
% BP_virt_theta = [BP_virt_theta  BP_virt_theta+pi/2];    %size 1*200. The first 100 have 0, and the second has pi/2 
% 
% for iSensor=nSensors+1:2*nSensors
%     BP_virt_names{iSensor} = ['Vertical Bp Virtual Sensor #' num2str(iSensor) ];
% end
% 
% sensor_btheta = fiesta_sensor_btheta( 'sensor', BP_virt_R, BP_virt_Z,BP_virt_theta, BP_virt_names );
% 
% %(r,z) of the sensors
% r_sensors=get(sensor_btheta,'r'); %size 1*200
% z_sensors=get(sensor_btheta,'z'); %size 1*200
% [R_sensors,Z_sensors]=meshgrid(r_sensors,z_sensors); %size 200*200
% 
% %Plot of the sensors
%     figure; hold on; axis equal;
%     plot(coilset);
%     contour( get(equil_eddy,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
%     contour( get(equil_eddy,'Psi'),get(equil_eddy,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
%     plot(vessel);
%     fileName = 'ST_target_equilibrium';
%     legend(gca,'hide');
%     plot(sensor_btheta);
%     %set(gca,'XLim',[0 1]);
%     %set(gca,'YLim',[-1.5 1.5]);
%     xlabel(gca,'R (m)');
%     ylabel(gca,'Z (m)');
%     title('Sensors to null Bpol with eddys');
%     % save_to_pdf( gcf, fileName );
%     %%%OPtionf for tfg
%     set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
%     axis([0,1.1,-1.1,1.1]) 

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


%% PASCHEN CURVE PLOT
 

C_1=[510 300]; %C_1 constant [ m^-1 Tor^-1] for H and He
C_2=[1.25e4 3.4e4];  %C_2 constant [V m^-1 Tor^-1]
Gas_type=["H_2","He"];
 p=linspace(1e-6,1e-3,100000);                           %[Tor]pressure of the prefill gas. size>1000 because if
                                                                            %not, It do not work properly
 Emin= @(L,p,C1,C2) C2*p./log(C1*p*L);
 
 %Plot H
 figure;
    %subplot(1,2,1)
    loglog(p,Emin(10,p,C_1(1),C_2(1)))
    hold on
    loglog(p,Emin(40,p,C_1(1),C_2(1)))
    loglog(p,Emin(100,p,C_1(1),C_2(1)))
    loglog(p,abs(E)*ones(1,length(p)),'k-')
    xlabel('Prefill pressure (Torr)')
    ylabel('E_{min} (V/m)')
    legend('L=10m','L=50m (GlobusM)','L=100m','Actual E')
    title(sprintf('Paschen curve, Gas=%s',Gas_type(1))); %d for numbers
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    %
    
    subplot(1,2,2)
    loglog(p,Emin(5,p,C_1(2),C_2(2)))
    hold on
    %
    loglog(p,Emin(10,p,C_1(2),C_2(2)))
    hold on
    loglog(p,Emin(40,p,C_1(2),C_2(2)))
    loglog(p,Emin(400,p,C_1(2),C_2(2)))
    loglog(p,Emin(1000,p,C_1(2),C_2(2)))
    loglog(p,abs(E)*ones(size(p)),'k-')
    xlabel('Prefill pressure (Torr)')
    ylabel('E_{min} (V/m)')
    legend('L=10m','L=40m','L=400m (GlobusM)','L=1000m','Actual E')
    title(sprintf('Paschen curve, Gas=%s',Gas_type(2))); %d for numbers
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
 
 %%%Estimation of the time when the avalancha has ended (begin
 %%%burn-through)

%Imagine the values choosen are
%E_choosen=5 %V/m
%p_choosen=10^-5 %Torr (multiplo de Pa)
K_B=1.380649*10^(-23); %J K^-1
 T_prefill=373; %K, temperature of the prefill gas

%     %Plot versus p, dado E
%     L_c=20 %L choosen
%     E_c=3 %[V/m] E choosen
%     p=linspace(1e-6,1e-2,10000);  %[Tor]
%     figure;
%     semilogx(p, tau_bd(p,E_c,L_c)*1e3)
%     xlabel('Prefill pressure (Tor)')
%     ylabel('Time when D_\alpha peaks (ms)')
%     title('Breakdown duration (starts at t=0, finish at the y axis time) L=70m')
% 
%     %Plot versus E, dado Ep
%     p_c=5e-4 %[Tor] p choosen
%     E=linspace(0.2,10);  %[Tor]
%     figure;
%     plot(E, tau_bd(p_c,E,L_c)*1e3)
%     xlabel('E (V/m)')
%     ylabel('Time when D_\alpha peaks (ms)')
%     title('Breakdown duration (starts at t=0, finish at the y axis time) L=70m')
%     
    %Plot3 varying E,p
    Ep=linspace(0.1,5e1,1000);  %[Tor]
    p=linspace(1e-5,1e-3,1000);  %[Tor]
    
    %Have to a meshgrid with p and E so that all the values of E are combined
    %with the p, and viceversa, which is was I want to see.
    [p,Ep]=meshgrid(p,Ep);
    
    %To plot this, will do as a I did with L inside the vessel. Since there
    %are problems with negative values, with substitute them with NaN. Have
    %to evaluate the function, which want only a single value of E and p
    %for not crashing or giving wrong things
    
    %Will do a loop to plot it for different L values
    clear tau
    tic
    L_plot=[5 10 70 90]             %[m]
    
    for i_L=1:length(L_plot)
    
        for i_p=1:length(p)
        
            for i_E=1:length(Ep)
            
                tau(i_p,i_E)=tau_bd(Ep(i_p,i_E),p(i_p,i_E),L_plot(i_L),C_1,C_2); %[s] time 
            end
        end
      time_Ep_loop=toc  
        %Plot
        figure;
        contourf(p,Ep,log10(tau*1e3));
        %contourf(p,E,tau);
        shading('interp') %this is to make the transition between values continuous,
                                    %instedad of discontinuously between pixels
       colormap(Gamma_II)
        hold on;
        colorbar %colorbar
        view(2) %2D view
        xlabel('p (Tor)')
        ylabel('E (V/m)')
        set(gca,'Xscale','log') %para poner el eje x en log   
        set(gca,'Yscale','log') %para poner el eje x en log 
        title(sprintf('log10(time radiation wal[ms]) for L=%d m',L_plot(i_L)))
        %title(sprintf('time radiation wal[ms]) for L=%d m',L_plot(i_L)))
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG

    end

  %Also, assuming the ionization fraction of 5%=0.05, the plasma current when
  %the line radiation peaks could be founded:
  
  E_ej=3%[V/m] example of field
  j_rad_wall=2*1.6*10^(-19)*1/(K_B*T_prefill)*E_ej*0.05 % [A/m^2] plasma current
                               %density when line rad peakd. Note that p cancels out
%asuming a surface of 2pi R, can find the current:
I_rad_wall= j_rad_wall*2*pi*0.2 %Pasma current [A], Assuming circular plasma of R=0.2m
    
 
 figure;
     p=linspace(1e-5,1e-3,10000);  %[Tor]
    loglog(p,Emin(5,p))
    hold on
    %
    loglog(p,Emin(10,p))
    loglog(p,Emin(70,p))
    loglog(p,Emin(90,p))
    loglog(p,Emin(50,p))
    xlabel('Prefill pressure (Torr)')
    ylabel('E_{min} (V/m)')
    legend('L=5m','L=10m','L=70m','L=90m','Globus-M, L=50m, V_{loop}=4.5-8V')
    title(sprintf('Paschen curve, Gas=%s',Gas_type)); %d for numbers
    set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
    axis([10^-5 10^-3 10^-1 50])
    
        
  %Plot of the mean free path%%%%%
  %The mean free path is lambda=1/C_1*p. Can easily plot it for the
  %interesting gases, H, He and Ar.
  
  lambda= @(C1,p) 1./(C1.*p); 
  figure;
  plot(p,lambda(C_1(1),p))
  hold on;
  plot(p,lambda(C_1(2),p))
  xlabel('p (Tor)')
  ylabel('\lambda')
  set(gca,'Xscale','log') %para poner el eje x en log   
  title('Mean free path \lambda versus p')
  legend(Gas_type(1),Gas_type(2))
  set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
  
 
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




%%
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Functions for the breakdown integration@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  
     %2) Field line integrator function with Lp
        %this solves the field line eq, using Lp (poloidal length) 
            %as an independent variable
        %way more faster than with phi (5 min when Lmax=3000m, 15 inside points)
    
    
    function [results]=Field_LineIntegrator_Lp(Lp,rzLphiU,Br_interpn,Bz_interpn,Bphi_interpn)
    %rzphiLU=[r z phi L U]
    %Lp= poloidal length (have to write capital L so it not apperas as
    %internsity I)
    
    %First, the field needs to be evaluated at the point (r,phi,z):
    
    Br_eval=Br_interpn(rzLphiU(1),rzLphiU(2));
    Bphi_eval=Bphi_interpn(rzLphiU(1),rzLphiU(2));
    Bz_eval=Bz_interpn(rzLphiU(1),rzLphiU(2)); 
    Bpol_eval=sqrt(Br_eval^2+Bz_eval^2);
    
    %With the field, the eq to solve is:
    
    dr_dLp=Br_eval/Bpol_eval;
    dphi_dLp=rzLphiU(1)*Bphi_eval/Bpol_eval;
    dz_dLp=Bz_eval/Bpol_eval;
    length=sqrt(1+(Bphi_eval/Bpol_eval)^2);
    U_Vloop=1/(2*pi*rzLphiU(1)); %pseudo potential U/V_loop
    
    results=zeros(5,1); %column vector to group the results
    results(1)=dr_dLp;
    results(4)=dphi_dLp;
    results(2)=dz_dLp;
    results(3)=length;
    results(5)=U_Vloop;
    
    end
    
    %3)EVENT FUNCTION FOR THE ODE

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

%%

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Function for the fields@@@@
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

%2) Earth Field

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

%%

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Function for tau_bd@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function [tau_radwall]=tau_bd(E,p,L,C_1,C_2)
    %This is the stimation for the time when the peak of line radiation occurs 
    %(Lloyd1991) 
    
    %If it is negative, I will set NaN, so that it do not plot it
    alpha=C_1*p.*exp(-C_2*p./E);    %First townsend coefficient
    
    if alpha-1/L<=0 %no avalanche
        tau_radwall=NaN;
    else%avalanche
        tau_radwall= 41./(43*E./p.*(alpha-1/L)); 
    end

    
end
%%
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Function for the colormap (Jose Rueda)@@@@
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
 