% this model evaluates the perceived positions and the number of sources
% inferred for multisensory stimuli presented in different spatial 
% configurations.
% IvMean, IaMean set up the efficacy of the first auditory and visual inputs
% IvMean2, IaMean2 set up the efficacy of the second auditory and visual inputs
% AV_Displacement identifies the AV distance
% this file has been used to generate data files utilized to create
% manuscript's figures n.4-11-12-13


clear, close all, clc

s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);

%% Model's parameters

% duration of the simulation
iter = 1000;
dt = 0.2;
t = [0:iter]*dt;
L = length(t);

% number of nerual elements in each region
dim_Areas = 180;
k = [1:1:dim_Areas];

% Parameters of the sigmoidal relationship defining single elements
G = 1;
phi = 20;
pend = 0.3;

tau_x=3;        % auditory time constant
tau_y=15;       % visual time constant
tau=1;          % causal inference layer time constant

%% CORTICAL LAYERS are defined

% VISUAL CORTICAL LAYER
Nv=dim_Areas;           % number of neural elements in this area
Iv = zeros(Nv,L);       % external input to the visual region
xv = zeros(Nv,L);       % activity of the visual region  
% ACOUSTIC CORTICAL LAYER
Na=dim_Areas;           % number of neural elements in this area
Ia = zeros(Na,L);       % external input to the auditory region
xa = zeros(Na,L);       % activity of the auditory region  
% CAUSAL INFERENCE AREA
Nm=dim_Areas;           % number of neural elements in this area
xm=zeros(Nm,L);         % activity of the causal inference region
xm_T=zeros(Nm,L);       % activity above the threshold of the causal inference region

% detection threshold in the causal inference layer
thr = 0.15;

%% FEEDFORWARD Projections are generated - % Eq.13 in the manuscript
% Auditory -> Multisensory
Wa=zeros(Na,Na);
Wa_SD=0.5;          % standard deviation
Wa0=18;             % W max
ka=[1:1:Na];
for i = 1:Na,
    
    distance_A=abs(ka-i);
    control_dist_A=double((abs(ka-i)>(Na/2))); 
    D2_A=((distance_A-sign(control_dist_A)*Na).^2);
    
    Wa(i,:)=Wa0*exp(-D2_A/2/Wa_SD/Wa_SD);
end 
% Visual -> Multisensory
Wv=zeros(Nv,Nv);
Wv_SD=0.5;          % standard deviation
Wv0=18;             % W max
kv=[1:1:Nv];
for i = 1:Nv,
    
    distance_V=abs(kv-i);
    control_dist_V=double((abs(kv-i)>(Nv/2))); 
    D2_V=((distance_V-sign(control_dist_V)*Nv).^2);
    
    Wv(i,:)=Wv0*exp(-D2_V/2/Wv_SD/Wv_SD);
end 
%% CROSSMODAL CONNECTIONS are generated - % Eq.11 in the manuscript
Wav = zeros(Na,Na);
Wva = zeros(Nv,Nv);
WAV0 = 1.4;         % W max
Wav_SD = 5;         % standard deviation
WVA0 = 1.4;         % W max
Wva_SD = 5;         % standard deviation
W = zeros(Nm,Nm);

for i = 1:dim_Areas,
    distance_M = abs(k-i);
    control_dist_M = double((abs(k-i)>(dim_Areas/2))); 
    
    D2 = ((distance_M-sign(control_dist_M)*dim_Areas).^2);

    Wav(i,:) = WAV0*exp(-D2/2/Wav_SD/Wav_SD); 
    Wva(i,:) = WVA0*exp(-D2/2/Wva_SD/Wva_SD);
end

%% INTRA_AREA Mexican Hat Synapses are generated - % Eqs.5-6 in the manuscript
% causal inference Area
LLm=zeros(Nm,Nm);
LLm_SD_in=10;           % standard deviation inhibitory                                                    
LLm_SD_ex=2;            % standard deviation excitatory
kL=3/2.6;
LLm0_in=2.6;            % L max inhibitory 
LLm0_ex=kL*LLm0_in;     % L max excitatory

k=[1:1:Nm];
for i = 1:Nm,
    distance_LLm=abs(k-i);
    control_dist_LLm=double((abs(k-i)>(Nm/2))); 
    
    D2_LLm=((distance_LLm-sign(control_dist_LLm)*Nm).^2);
    LLm(i,:)=LLm0_ex*exp(-D2_LLm/2/LLm_SD_ex/LLm_SD_ex)-LLm0_in*exp(-D2_LLm/2/LLm_SD_in/LLm_SD_in); 
    LLm(i,i)=0;    
    
end

% Auditory Area
LLa=zeros(Na,Na);
k_SD=3/120;
LLa_SD_in=120;              % standard deviation inhibitory 
LLa_SD_ex=k_SD*LLa_SD_in;   % standard deviation excitatory       
k_L=5/4;
LLa0_in=4;                  % L max inhibitory
LLa0_ex=k_L*LLa0_in;        % L max excitatory

k=[1:1:Na];
for i = 1:Na,
    distance_LLa=abs(k-i);
    control_dist_LLa=double((abs(k-i)>(Na/2))); 
    
    D2_LLa=((distance_LLa-sign(control_dist_LLa)*Na).^2);
    LLa(i,:)=LLa0_ex*exp(-D2_LLa/2/LLa_SD_ex/LLa_SD_ex)-LLa0_in*exp(-D2_LLa/2/LLa_SD_in/LLa_SD_in); 
    LLa(i,i)=0;    
    
end

% Visual Area
LLv=zeros(Nv,Nv);
k_SD=3/120;
LLv_SD_in=120;              % standard deviation inhibitory
LLv_SD_ex=k_SD*LLv_SD_in;   % standard deviation excitatory     
k_L=5/4;
LLv0_in=4;%                 % L max inhibitory       
LLv0_ex=k_L*LLv0_in;%       % L max excitatory

k=[1:1:Nv];
for i = 1:Nv,
    distance_LLv=abs(k-i);
    control_dist_LLv=double((abs(k-i)>(Nv/2))); 
    
    D2_LLv=((distance_LLv-sign(control_dist_LLv)*Nv).^2);
    LLv(i,:)=LLv0_ex*exp(-D2_LLv/2/LLv_SD_ex/LLv_SD_ex)-LLv0_in*exp(-D2_LLv/2/LLv_SD_in/LLv_SD_in); 
    LLv(i,i)=0;    
    
end

%% AUDITORY AND VISUAL INPUTS are defined
% Visual and Auditory Inputs
sy = 4;            % Visual Standard Deviation SD
sx = 32;           % Auditory Standard Deviation SD
% Insert intensity of Input 1
IvMean = 0;%
IaMean = 28;%
% Insert intensity of Input 2
IvMean2 = 27;
IaMean2 = 0;

% Added noise to the visual and auditory inputs
I_Noise_V = -(IvMean*0.4)+(2*IvMean*0.4).*rand(dim_Areas,1);
I_Noise_A = -(IaMean*0.4)+(2*IaMean*0.4).*rand(dim_Areas,1);

E0y = IvMean;
E0x = IaMean;
E0y2 = IvMean2;
E0x2 = IaMean2;

% Inputs duration
T_A=500;
T_V=500;

% Distance in degrees between the inputs
AV_Displacement = -10;
% Set the Position of Input 1
position=90;
position_v = position; 
position_a = position;
% Set the Position of Input 2
position_v2 = position_v+AV_Displacement; 
position_a2 = position_a+AV_Displacement;

%% AUDITORY AND VISUAL INPUTS are generated in this section - Eqs.8-9 in the manuscript

distance_Ia=abs(k-position_a);
control_dist_Ia=double((abs(k-position_a)>(dim_Areas/2)));
D2=((distance_Ia-sign(control_dist_Ia)*dim_Areas).^2);
% Gaussian function defining the First Auditory Input
Ix0=E0x*exp(-D2/2/sx/sx);       

distance_Ia2=abs(k-position_a2);
control_dist_Ia2=double((abs(k-position_a2)>(dim_Areas/2)));
D2_2=((distance_Ia2-sign(control_dist_Ia2)*dim_Areas).^2);
% Gaussian function defining the Second Auditory Input
Ix02=E0x2*exp(-D2_2/2/sx/sx);

distance_Iv=abs(k-position_v);
control_dist_Iv=double((abs(k-position_v)>(dim_Areas/2)));
D2=((distance_Iv-sign(control_dist_Iv)*dim_Areas).^2);
% Gaussian function defining the First Visual Input
Iy0=E0y*exp(-D2/2/sy/sy);

distance_Iv2=abs(k-position_v2);
control_dist_Iv2=double((abs(k-position_v2)>(dim_Areas/2)));
D2_2=((distance_Iv2-sign(control_dist_Iv2)*dim_Areas).^2);
% Gaussian function defining the Second Visual Input
Iy02=E0y2*exp(-D2_2/2/sy/sy);

A_Onset = 0/dt;
A_Duration = T_A/dt; 

V_Onset = 0/dt;
V_Duration = T_V/dt; 

Diff_AV = (V_Duration + V_Onset) - (A_Duration + A_Onset);
rest = iter/dt - V_Duration - V_Onset;

% Inputs evolution over time during the simulation
Input_step_x = [zeros(1,A_Onset) ones(1,A_Duration) zeros(1,rest+Diff_AV)];
Input_step_y = [zeros(1,V_Onset) ones(1,V_Duration) zeros(1,rest)];

%% Main Body of the model

ux=zeros(Na,L);     % overall input to the auditory area
uy=zeros(Nv,L);     % overall input to the visual area
uVxm=zeros(Nm,L);   % visual input to the causal inference area
uAxm=zeros(Nm,L);   % auditory input to the causal inference area
Sm=zeros(Nm,L);     % Intra_area competition in the causal inference area
um=zeros(Nm,L);     % overall input to the causal inference area

for k=1:L-1,
    dt*k;
    % Eq.7 in the manuscript
    Iv(:,k)=(Iy0*Input_step_y(k))'+(Iy02*Input_step_y(k))'+I_Noise_V(:,1);
    Ia(:,k)=(Ix0*Input_step_x(k))'+(Ix02*Input_step_x(k))'+I_Noise_A(:,1);
    
    % INPUTS TO LAYERS
    % VISUAL LAYER
    % Eq.4 in the manuscript
    input_lateral_v = LLv*xv(:,k);
    % Eq.10 in the manuscript
    input_crossmodal_v = Wva*xa(:,k);
    % AUDITORY LAYER
    % Eq.4 in the manuscript
    input_lateral_a = LLa*xa(:,k);
    % Eq.10 in the manuscript
    input_crossmodal_a = Wav*xv(:,k);
    % CAUSAL INFERENCE LAYER
    % Eq.12 in the manuscript
    Im_v = Wv*xv(:,k);
    Im_a = Wa*xa(:,k);
    % Eq.4 in the manuscript
    Sm = LLm*xm(:,k);
    

    % UNISENSORY LAYERS
    % Eq.3 in the manuscript
    ux(:,k) = Ia(:,k) + input_crossmodal_a + input_lateral_a;
    % Eqs.1-2 in the manuscript
    xa(:,k+1) = xa(:,k) + dt*((1/tau_x)*(-xa(:,k)+1./(1+exp(-(ux(:,k)-phi)*pend))));
    
    % Eq.3 in the manuscript
    uy(:,k) = Iv(:,k) + input_crossmodal_v + input_lateral_v;
    % Eqs.1-2 in the manuscript
    xv(:,k+1) = xv(:,k) + dt*((1/tau_y)*(-xv(:,k)+1./(1+exp(-(uy(:,k)-phi)*pend))));
    
    uVxm(:,k) = Im_v; 
    uAxm(:,k) = Im_a; 
    
    % MULTISENSORY LAYER
    % Eq.3 in the manuscript
    um(:,k)= uVxm(:,k) + uAxm(:,k) + Sm;
    % Eqs.1-2 in the manuscript
    xm(:,k+1)=xm(:,k)+dt*((1/tau)*(-xm(:,k)+1./(1+exp(-(um(:,k)-phi)*pend))));    
    
    xm_T(:,k+1)=(xm(:,k+1)-thr).*(xm(:,k+1)>thr);

end

%% BARYCENTER Method in the Unisensory Areas to Evaluate the Auditory and Visual Percepts

if position_a < 90
    abscissa_x = [ [1:1:position_a+90] [position_a-90:1:-1] ];
end
if position_a > 90
        abscissa_x = [ [181:1:position_a+90] [position_a-89:1:180] ];
end
if position_a == 90
    abscissa_x = 1:1:180;
end

if position_v < 90
    abscissa_y = [ [1:1:position_v+90] [position_v-90:1:-1] ];
end
if position_v > 90
        abscissa_y = [ [181:1:position_v+90] [position_v-89:1:180] ];
end
if position_v == 90
    abscissa_y = 1:1:180;
end

xa_T=xa(:,L);
xv_T=xv(:,L);
% Perceived Position Auditory Stimulus is evalutated here
A_Percept = sum(xa_T.*abscissa_x')/sum(xa_T)
% Perceived Position Visual Stimulus is evalutated here
V_Percept = sum(xv_T.*abscissa_y')/sum(xv_T)

% Auditory and Visual Bias
A_Bias = (((A_Percept-position_a))/abs(position_a-position_v))*100;
V_Bias = (((V_Percept-position_v))/abs(position_v-position_a))*100;

A2_Bias = (((A_Percept-position_a2))/abs(position_a2-position_v2))*100;
V2_Bias = (((V_Percept-position_v2))/abs(position_v2-position_a2))*100;

%% Causal Inference
M(1,:) = xm(:,L);
M_thr(1,:) = (M(1,:)-thr).*(M(1,:)>thr);
M_index = find(M_thr(1,:)>0);                      
if sum(M_thr(1,min(M_index):max(M_index))<=0) == 0 && not(isempty(M_index)),
    
    % the network identifies a common cause
    S = 1;
    
elseif sum(M_thr(1,min(M_index):max(M_index))<=0) > 0,
    
    % the network identifies independent causes
    S = 2;
    
else
    
    % the network does not identify any input source
    S = 0;
    
end
S

%% BARYCENTER Method in the Causal Inference Area to Evaluate the Position of the Input Sources
if position < 90
    abscissa_m = [ [1:1:position+90] [position-90:1:-1] ];
end
if position > 90
        abscissa_m = [ [181:1:position+90] [position-89:1:180] ];
end
if position == 90
    abscissa_m = 1:1:180;
end

index_S1 = 0;
index_S2 = 0;
iS2 = 1;
iS1 = 1;
if S == 1,
    position_source1_S1 = sum(M_thr(1,:).*abscissa_m)/sum(M_thr(1,:)); 
    position_source2_S1 = position_source1_S1;                          
elseif S == 2,
    for j = 1 : length(M_index),
        if (j == 1) || ((M_index(j-1)+1 == M_index(j)) && (index_S1(1) == 0)),
            
            index_S2(iS2) = M_index(j);
            iS2 = iS2 + 1;
        elseif M_index(j) > 0,
            
            index_S1(iS1) = M_index(j);
            iS1 = iS1 + 1;
        end
    end
    position_source1_S2 = sum(M_thr(1,index_S1).*index_S1(1:length(index_S1)))/sum(M_thr(1,index_S1)); 
    position_source2_S2 = sum(M_thr(1,index_S2).*index_S2(1:length(index_S2)))/sum(M_thr(1,index_S2));
end


%% Figure
figure
for l=1:5:L-1
    subplot(3,1,1)
    imagesc(xm(:,l)')
    caxis([0 1])
    title('causal inference area')
    
    subplot(3,1,2)
    imagesc(xv(:,l)')
    caxis([0 1])
    title('visual layer')

    subplot(3,1,3)
    imagesc(xa(:,l)')
    caxis([0 1])
    title('auditory layer')
   
    pause(0.01)
end

figure
for l=1:5:L-1
    subplot(3,1,1)
    plot(xm(:,l)')
    ylim([0 1])
    xlim([1 dim_Areas])
    title('causal inference area')
    
    subplot(3,1,2)
    plot(xv(:,l)')
    ylim([0 1])
    xlim([1 dim_Areas])
    title('visual layer')

    subplot(3,1,3)
    plot(xa(:,l)')
    ylim([0 1])
    xlim([1 dim_Areas])
    title('auditory layer')
   
    pause(0.01)
end
