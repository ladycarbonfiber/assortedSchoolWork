%ECE 4445 Design Project Appendix
%Tom Wells



%Read in Z_off
%data set one is from driver in free space
%Z_off = importdata('12_off.txt');

%Read in Z_on
%data set two is from the driver in a box of 1ft^3
%used to calculate VAS
%Z_on = importdata('12_on.txt');

data_obs_off = Z_off.data;
fq_obs_off = data_obs_off(:,1);
mag_obs_off = data_obs_off(:,2);
loglog(fq_obs_off, mag_obs_off);
fq = fq_obs_off;
xlabel('Frequency, Hz');
ylabel('Magnitude');
title('Measured Driver Impedence');
figure

RE = 5.16; %given
RES = 42.28047 - RE; %peak value observed
R1 = sqrt(RE*(RES+RE)); %geometric mean
%we exploit the fact that due to the low fq we can neglect the inductance
%and treat the transfer function as symetric around the resonance
fL = 21.36206; %observed at near R1;
fU = 50.56323; %%upper observed @ near R1;
fs = 33.73340; %fq @ resonance
ws = 2*pi*fs; %omega simplifies equation below
%the mechanical quality factor is calculated based on simplification of 
% the text's equation 11.112.
QMS = fs *sqrt((RE+RES)/RE)/(fU-fL);
%follows simply from QMS
QES = QMS*(RE/RES);
QTS = QMS *(RE/(RE+RES));

%Box parameters used to calculate VAS
fs_box = 61.95850; %shifted right but still before inductive regime
RES_box = 42.38421 - RE; %observed resonance impedance 
ws_box = 2*pi*fs_box;
R1_box = sqrt(RE*(RES_box+RE));
fL_box = 48.07739; %observed
fU_box = 75.87695;
QMS_box = fs_box*sqrt((RE+RES_box)/RE)/(fU_box-fL_box);
%calculate QMS for box mounted speaker in the same manner as above
QES_box = QMS_box *(RE/RES_box);
VT = .0283168; %given  1ft^3 as m^3
VAS = VT*((fs_box/fs)*(QES_box/QES) -1);%given


a = 5*.0254; %inches to m
SD = a^2*pi;%effective piston area of driver
%knowing VAS we can now calculate CMS and MMS
p0 = 1.17; %air density factor
c = 345;%speed of sound 
%Curcuit parameters
CMS = VAS/(p0*c^2*SD^2);%= ~12.5uF m/N
CAS =VAS/(p0*c^2); % = 516nF
MMS = 1/((ws^2)*CMS);%=~ .1107 kg
MAS = 1/((ws^2)*CAS);%= ~43

Bl = sqrt((RE/QES)*sqrt(MMS/CMS)); %64.93T-m
MA1_infinite_baffal = ((8*p0)/(3*pi^2*a)); %=1.24 kg/m^4
RAS = (1/QMS)*sqrt(MMS/CMS);%=~3.2 N s /m
MMD = (MAS - 2*MA1_infinite_baffal)*SD^2;
MAD = (MAS - 2*MA1_infinite_baffal);
RA1 = (.4410 *p0*c)/SD;%=3.5k
RA2 = (p0*c)/SD;%=8k
CA1 = 5.94*a^3/(p0*c^2);%=87nF



data_obs_on = Z_on.data;
fq_obs_on = data_obs_on(:,1);
mag_obs_on = data_obs_on(:,2);
%figure
%loglog(fq_obs_on, mag_obs_on);

%two fqs are polled to calculate n and LEE
%Z1, Z2 picked arbritraily to be in the lossless and lossy regimes
Z1 = 48.15852;
Z2 = 70.77856;
%frequencies corresponding
w1 = 2*pi*9586.98047;
w2 = 2*pi*17612.38379;
%linear regression mapped to loglog scale
n_raw = log((1/Z1)/(1/Z2))/log(w2/w1);%~.63
LEE_raw = cos(n_raw*(pi/2))/((1/Z1)*w1^n_raw);%~.0247
LE_raw = 1/((sin(n*pi/2)/(.5*(w1^n_raw)*LEE_raw))*w1*(1/Z1));%~.02 Im not sure if I did that right
le_raw = LE_raw*j*2*pi.*fq;%lossless
lee_raw = LEE_raw*((j*2*pi.*fq).^n_raw);%lossy
inductive_part_raw = (le_raw .*lee_raw)./(le_raw+lee_raw);
%tweaked values to get the graph to line up
n=.66;
LE =.01;%lossless part standard inductor as LE*jw
LEE = .036; %lossy part modeled as LEE*jw^n
le = LE*j*2*pi.*fq;%lossless
lee = LEE*((j*2*pi.*fq).^n);%lossy
inductive_part = (le .*lee)./(le+lee);

s = 2*pi*1i .* fq;
Z_n = (1/(QMS)).*(s./ws);%(1/(QMS)).*(s./ws);
Z_d = (((1/ws).*s).^2) + (1/(QMS)).*(s./ws) + 1;
Z_mot = RES .* (Z_n ./ Z_d);
%Z_T = tf([1/(QMS*ws),0],[1/ws^2 1/(QMS*ws) 1]);
loglog(fq, abs(Z_mot));
xlabel('Frequency, Hz');
ylabel('Magnitude');
title('Driver Motional Impedance');

model_T_beforeTweak = abs((RE + inductive_part_raw + Z_mot));
model_T = abs((RE + inductive_part + Z_mot));

%------------------------------------------------------------------------
%-------------------Table of Small Signal Parameters---------------------
%----fs = 33.73
%----QTS = 0.4036
%----QMS = 3.3068
%----QES = 0.4597
%----VAS = 2.5381
%------------------------------------------------------------------------
figure
loglog(fq_obs_off,mag_obs_off);%refernce curve
hold on
loglog(fq,model_T_beforeTweak,'g');
loglog(fq, model_T,'r');

legend('Measured','Calculated', 'Fitted');
xlabel('Frequency, Hz');
ylabel('Magnitude');
title('Modeled Driver Impedence Comparisons');


%ZEL plots
figure
%parallel combination
loglog(fq, abs(inductive_part));
%lossless
hold on
loglog(fq, abs(le),'r');
%lossy
loglog(fq, abs(lee),'g');
%measured
ZEL = mag_obs_off - Z_mot - RE;
%ZEL(ZEL<.5) =.01; %clean up noise
loglog(fq, ZEL,'m');
legend('Parallel Combination','Lossless', 'Lossy','Measured');
xlabel('Frequency, Hz');
ylabel('Magnitude');
title('Modeled Inductor Impedence Comparisons');

%Circuit Comparison
%figure
%zMot2 = ((Bl^2)/SD^2) ./ (MAD + RAS + (1./(CAS * 2* pi * j .*fq)) + 2*MA1_infinite_baffal*2*pi*j .*fq);
%model_T2 = abs((RE + inductive_part + zMot2));
%loglog(fq, model_T, 'r')
%hold on
%loglog(fq, model_T2)
