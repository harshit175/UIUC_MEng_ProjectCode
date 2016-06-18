clear all

global n1 n2 n3 n m multfactor

rc=0.00e-3; %center radius
rin=0.04e-3; %wire radius

rinter=.17e-3; %microchannel radius
rout=.5e-3;   %outer radius, rout-rinter is the PDMS thickness
%rout==rinter implies there is no PDMS

% rc=60e-3;
% rin=60e-3; %inner radius
% %rinter=0.17e-3; %interface radius
% rinter=60.01e-3;
% rout=60.5e-3;   %outer radius


Tambient=30; %this is the initial temperature
Tb=200;  %this is ignition temperature

alphainitial = 0.00;  %this is the initial cure

%conductivity, (c suffix added to k, to distinguish from k)
kc=[.19 .27 401]; %first entry is polymer (Frulloni), second is PDMS, third is wire

%wire-polymer-PDMS (order in which the three occur, from left to right)
%n3-n1-n2

%number of radial elements in wire are denoted by n3
%number of radial elements in microchannel are denoted by n1
%number of radial elements in PDMS are denoted by n2

%polymer Frulloni
rho1=1.19e3;
cp1=3118;

%polymer TMPTA
%rho1=1.060e3;
%cp1=1700;

%PDMS properties
rho2=1.03e3;
cp2=1100;

% stainless steel
% k=16
% rho3=8.06e3;
% cp3=490;

%cadmium
% kc(3)=92;
% rho3=8.65e3;
% cp3=230;

%copper properties
rho3=8.92e3;
cp3=386;

Hr=550e3; %heat of reaction in Frulloni paper

%Hr = 360e3; %Heat of reaction of TMPTA


%density and enthalpy, lumped together as rhocp
rhocp=[rho1*cp1 rho2*cp2 rho3*cp3]; 

Q=0e13;      %resistive heating W/m^3 ( 1e13 leads to quick heating)
q=008e6;     %heat supplied to start the process (classic value is 8e6) W/m^2
%q=0;

%if q is not zero than the bottom boundary temperature (Tb) specification
%becomes redundant 

%hz=15e-3; %height in z direction;
hz=.2e-3;

%elements along radial direction
n1=double(int64((rinter-rin)*1e3*200)); %elements in microchannel
n3=double(int64((rin-rc)*1e3*200));     %elements in conducting wire

gpratio=1.1; %the geometric progression ratio used for radial length of PDMS elements
%if gpratio is specified as 1, n2 has to be specified in the solver file.

beta=0.5; %used in implicit scheme
curelimit=0.95;
%elements along z direction

m=hz*1e3*200;

%time step in seconds;
dt=4e-6;


te=1000  %time step size

%the Temperature and cure fields are recorded at certain time intervals,
%denoted by recordstep
recordstep=10;


%Data input ends


