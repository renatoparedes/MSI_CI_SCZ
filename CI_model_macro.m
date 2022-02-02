% this macro evaluates the perceived positions and the number of sources
% inferred for multisensory stimuli presented in different spatial 
% configurations.
% n_iter sets up the number of simulations run with the same spatial
% configuration of the external stimuli
% min_AV_distance and max_AV_distance identify the range of the AV
% distances
% this file has been used to generate data files utilized to create
% manuscript's figures n.2-3-5-6-7-8-9-10

clear, close all, clc

s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);

% Model's parameters

% duration of the simulation
iter = 1500;
dt = 0.2;
t = [0:iter]*dt;
L = length(t);
n_iter = 1000;      % number of repetitions for each stimulus configuration

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


%% FEEDFORWARD Projections are generated
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

%% CROSSMODAL CONNECTIONS are generated
Wav = zeros(Na,Na);
Wva = zeros(Nv,Nv);
WAV0 = 1.4;%
WAV0_SD = 5;
WVA0 = 1.4;%
WVA0_SD = 5;
k = [1:1:dim_Areas];
for i = 1:dim_Areas,
    distance = abs(k-i);
    control_dist = double((abs(k-i)>(dim_Areas/2))); 
    D2 = ((distance-sign(control_dist)*dim_Areas).^2);
    Wav(i,:) = WAV0*exp(-D2/2/WAV0_SD/WAV0_SD); 
    Wva(i,:) = WVA0*exp(-D2/2/WVA0_SD/WVA0_SD);
end

%% Intra_Area Synapses - Competitive Mechanism - Mexican Hat Synapses are generated
% Multisensory Area
LLm=zeros(Nm,Nm);
LLm_SD_in=10;                                                    
LLm_SD_ex=2;
kL=3/2.6;
LLm0_in=2.6; 
LLm0_ex=kL*LLm0_in;
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
LLa_SD_in=120; 
LLa_SD_ex=k_SD*LLa_SD_in;       
k_L=5/4;
LLa0_in=4;         
LLa0_ex=k_L*LLa0_in;  
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
LLv_SD_in=120;
LLv_SD_ex=k_SD*LLv_SD_in;     
k_L=5/4;
LLv0_in=4;%         
LLv0_ex=k_L*LLv0_in;%
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
sy = 4;             % Visual Standard Deviation SD
sx = 32;            % Auditory Standard Deviation SD
% Insert intensity of Inputs
IvMean = 27;        % Visual Input
IaMean = 28;        % Auditory Input
% Added noise to the visual and auditory inputs
I_Noise_V = -(IvMean*0.4)+(2*IvMean*0.4).*rand(dim_Areas,n_iter);
I_Noise_A = -(IaMean*0.4)+(2*IaMean*0.4).*rand(dim_Areas,n_iter);
% Inputs duration
T_A=500;
T_V=500;

% Set the Position of the Visual Input
position=50;
% set the analyzed displacements between the auditory and visual stimuli (here min AV_distance = 2°, max AV_distance = 20°)
min_AV_distance = 2;
max_AV_distance = 20;
AV_distance_step = 2;

for d = min_AV_distance:AV_distance_step:max_AV_distance, 
    
    A_Displacement = d;     % Displacement of the auditory input
    
    S1 = 0;                 % counter common source evaluation (how many times the model identifies C=1)
    S2 = 0;                 % counter independent sources evaluation (how many times the model identifies C=2)
    
    % Positions of external stimuli
    position_a = position + A_Displacement;
    position_v = position; 
    
    position_source1_S1 = zeros(1,n_iter);
    position_source2_S1 = zeros(1,n_iter);
    position_source1_S2 = zeros(1,n_iter);
    position_source2_S2 = zeros(1,n_iter);
    
    M = zeros(n_iter,dim_Areas);
    M_thr = zeros(n_iter,dim_Areas);
    V = zeros(n_iter,dim_Areas);
    A = zeros(n_iter,dim_Areas);
    
    % Variables used to generate manuscript's figures
    A_Bias_S1=[];
    A_Bias_S2=[];
    V_Bias_S1=[];
    V_Bias_S2=[];
    V_Bias_TOT=[];
    A_Bias_TOT=[];
    A_Bias_S1_deg=[];
    A_Bias_S2_deg=[];
    V_Bias_S1_deg=[];
    V_Bias_S2_deg=[];
    V_Bias_TOT_deg=[];
    A_Bias_TOT_deg=[];
    
    for ii = 1:n_iter,
        
        %% main corp - single simulation repeated n_iter times for each AV distance analyzed
        E0x = IaMean;  
        E0y = IvMean;
        k=[1:1:dim_Areas];
        
        distance=abs(k-position_a);
        control_dist=double((abs(k-position_a)>(dim_Areas/2)));
        D2=((distance-sign(control_dist)*dim_Areas).^2);
        % Gaussian function defining the Auditory Input
        Ix0=E0x*exp(-D2/2/sx/sx);
        
        distance=abs(k-position_v);
        control_dist=double((abs(k-position_v)>(dim_Areas/2)));
        D2=((distance-sign(control_dist)*dim_Areas).^2);
        % Gaussian function defining the Visual Input
        Iy0=E0y*exp(-D2/2/sy/sy);
        
        A_Onset = 0/dt;             % set the Onset of the auditory stimulus
        A_Duration = T_A/dt;        % set the duration of the auditory stimulus
        
        V_Onset = 0/dt;             % set the Onset of the visual stimulus
        V_Duration = T_V/dt;        % set the duration of the visual stimulus
        
        Diff_AV = (V_Duration + V_Onset) - (A_Duration + A_Onset);
        rest = iter/dt - V_Duration - V_Onset;
        
        Input_step_x = [zeros(1,A_Onset) ones(1,A_Duration) zeros(1,rest+Diff_AV)];
        Input_step_y = [zeros(1,V_Onset) ones(1,V_Duration) zeros(1,rest)];
        
        ux=zeros(Na,L);     % overall input to the auditory area
        uy=zeros(Nv,L);     % overall input to the visual area
        uVxm=zeros(Nm,L);   % visual input to the causal inference area
        uAxm=zeros(Nm,L);   % auditory input to the causal inference area
        Sm=zeros(Nm,L);     % Intra_area competition in the causal inference area
        
        for k=1:L-1,
            
            Iv(:,k)=(Iy0*Input_step_y(k))'+I_Noise_V(:,ii);
            Ia(:,k)=(Ix0*Input_step_x(k))'+I_Noise_A(:,ii);
            
            % INPUTS TO LAYERS
            % visual
            input_lateral_v = LLv*xv(:,k);
            input_crossmodal_v = Wva*xa(:,k);
            % auditory
            input_lateral_a = LLa*xa(:,k);
            input_crossmodal_a = Wav*xv(:,k);
            % multisensory
            Im_v = Wv*xv(:,k);
            Im_a = Wa*xa(:,k);
            Sm = LLm*xm(:,k);
            
            % UNISENSORY LAYERS
            ux(:,k) = Ia(:,k) + input_crossmodal_a;
            xa(:,k+1) = xa(:,k) + dt*((1/tau_x)*(-xa(:,k)+1./(1+exp(-(ux(:,k)+input_lateral_a-phi)*pend))));

            uy(:,k) = Iv(:,k) + input_crossmodal_v;
            xv(:,k+1) = xv(:,k) + dt*((1/tau_y)*(-xv(:,k)+1./(1+exp(-(uy(:,k)+input_lateral_v-phi)*pend))));
            
            uVxm(:,k) = Im_v; 
            uAxm(:,k) = Im_a; 
            
            % MULTISENSORY LAYER
            Slm(:,k)=Sm;
            xm(:,k+1)=xm(:,k)+dt*((1/tau)*(-xm(:,k)+1./(1+exp(-(uAxm(:,k)+uVxm(:,k)+Slm(:,k)-phi)*pend))));%+I_Noise_M(:,ii)
            
            xm_T(:,k+1)=(xm(:,k+1)-thr).*(xm(:,k+1)>thr);
            
        end
        %%%%%% end main corp %%%%%%
        
        %%
        M(ii,:) = xm(:,L);
        A(ii,:) = xa(:,L);
        V(ii,:) = xv(:,L);
        %% BARYCENTER Method in the unisensory Areas
        
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
        % Perceived Position Auditory Stimulus is evalutated here
        A_Percept = sum(xa(:,L).*abscissa_x')/sum(xa(:,L));
        % Perceived Position Visual Stimulus is evalutated here
        V_Percept = sum(xv(:,L).*abscissa_y')/sum(xv(:,L));
        
        %% Causal Inference
        if position < 90
            abscissa_m = [ [1:1:position+90] [position-90:1:-1] ];
        end
        if position > 90
            abscissa_m = [ [181:1:position+90] [position-89:1:180] ];
        end
        if position == 90
            abscissa_m = 1:1:180;
        end
        M_thr(ii,:) = (M(ii,:)-thr).*(M(ii,:)>thr); 
        M_index = find(M_thr(ii,:)>0);
        
        % control if the network perceives at least one source
        is_source_null=sum(M_thr(ii,min(M_index):max(M_index)));            
        
        if sum(M_thr(ii,min(M_index):max(M_index))<=0) == 0 && is_source_null>0,
            
            % the network identifies a common cause
            S(ii) = 1;
            S1 = S1+1;      % counter common source evaluation (how many times the model identifies C=1)
        
        elseif sum(M_thr(ii,min(M_index):max(M_index))<=0) > 0 && is_source_null>0,
            
            % the network identifies independent causes
            S(ii) = 2;
            S2 = S2+1;      % counter independent sources evaluation (how many times the model identifies C=2)
        
        else
            
            % the network does not identify any input source
            S(ii) = 0;
        
        end
        
        index_S1 = 0;
        index_S2 = 0;
        iS2 = 1;
        iS1 = 1;
        
        if S(ii) == 1,      % Perception Bias computed in case of C=1
            
            position_source1_S1(ii) = sum(M_thr(ii,:).*abscissa_m)/sum(M_thr(ii,:));
            position_source2_S1(ii) = position_source1_S1(ii);
            if position_a~=position_v,
                
                A_Bias_S1(S1) = (((A_Percept-position_a))/abs((position_a-position_v)))*100; % Percentage Auditory Bias
                V_Bias_S1(S1) = (((V_Percept-position_v))/abs((position_v-position_a)))*100; % Percentage Visual Bias
                
            end
            
            A_Bias_S1_deg(S1) = abs(A_Percept-position_a); % Auditory Bias in degrees
            V_Bias_S1_deg(S1) = abs(V_Percept-position_v); % Visual Bias in degrees
            
        elseif S(ii) == 2,  % Perception Bias computed in case of C=2
            
            for j = 1 : length(M_index),
                if (j == 1) || ((M_index(j-1)+1 == M_index(j)) && (index_S1(1) == 0)),
                    % source n.2
                    index_S2(iS2) = M_index(j);
                    iS2 = iS2 + 1;
                elseif M_index(j) > 0,
                    % source n.1
                    index_S1(iS1) = M_index(j);
                    iS1 = iS1 + 1;
                end
            end
            position_source1_S2(ii) = sum(M_thr(ii,index_S1).*index_S1(1:length(index_S1)))/sum(M_thr(ii,index_S1));
            position_source2_S2(ii) = sum(M_thr(ii,index_S2).*index_S2(1:length(index_S2)))/sum(M_thr(ii,index_S2));
            if position_a~=position_v,
                
                A_Bias_S2(S2) = (((A_Percept-position_a))/abs((position_a-position_v)))*100;  % Percentage Auditory Bias  
                V_Bias_S2(S2) = (((V_Percept-position_v))/abs((position_v-position_a)))*100;  % Percentage Visual Bias
                
            end
            
            A_Bias_S2_deg(S2) = abs(A_Percept-position_a);  % Auditory Bias in degrees
            V_Bias_S2_deg(S2) = abs(V_Percept-position_v);  % Visual Bias in degrees
            
        end
    end
    
    % end of n_iter simulations for a specific AV distance
    % In the following we evaluate mean values of the auditory and visual perception Bias (in general, subdivided in the two cases C1 and C2) for
    % the current spatial configuration of the multisensory stimuli
    
    if position_a==position_v,      % In case of AV distance = 0 the Percentage Bias cannot be computed, only the Bias in degrees is reported 
        if S2==0,
            
            % case of only C=1 identifications
            A_Bias_TOT=0;
            V_Bias_TOT=0;
            A_Bias_TOT_deg=[A_Bias_S1_deg];
            V_Bias_TOT_deg=[V_Bias_S1_deg];
            
            Mean_position_source1_S1 = sum(position_source1_S1,2)/S1;
            Mean_position_source1_S2 = 0;
            Mean_position_source2_S1 = sum(position_source2_S1,2)/S1;
            Mean_position_source2_S2 = 0;
            
            A_Bias_S1 = 0;
            V_Bias_S1 = 0;
            A_Bias_S2 = 0;
            V_Bias_S2 = 0;
            A_Bias_S1_M = 0;
            A_Bias_S1_SEM = 0;
            A_Bias_S2_M = 0;
            A_Bias_S2_SEM = 0;
            V_Bias_S1_M = 0;
            V_Bias_S1_SEM = 0;
            V_Bias_S2_M = 0;
            V_Bias_S2_SEM = 0;
            
        elseif S2>0 && S1>0, 
            
            % case with C=1 and C=2 identifications
            A_Bias_TOT=0;
            V_Bias_TOT=0;
            A_Bias_TOT_deg=[A_Bias_S1_deg A_Bias_S2_deg];
            V_Bias_TOT_deg=[V_Bias_S1_deg V_Bias_S2_deg];
            
            Mean_position_source1_S1 = sum(position_source1_S1,2)/S1;
            Mean_position_source1_S2 = sum(position_source1_S2,2)/S2;
            Mean_position_source2_S1 = sum(position_source2_S1,2)/S1;
            Mean_position_source2_S2 = sum(position_source2_S2,2)/S2;
            
            A_Bias_S1 = 0;
            V_Bias_S1 = 0;
            A_Bias_S2 = 0;
            V_Bias_S2 = 0;
            A_Bias_S1_M = 0;
            A_Bias_S1_SEM = 0;
            A_Bias_S2_M = 0;
            A_Bias_S2_SEM = 0;
            V_Bias_S1_M = 0;
            V_Bias_S1_SEM = 0;
            V_Bias_S2_M = 0;
            V_Bias_S2_SEM = 0;
            
        else % S2>0 && S1==0,  
            
            % case with only C=2 identifications
            A_Bias_TOT=0;
            V_Bias_TOT=0;
            A_Bias_TOT_deg=[A_Bias_S2_deg];
            V_Bias_TOT_deg=[V_Bias_S2_deg];
            
            Mean_position_source1_S1 = 0;
            Mean_position_source1_S2 = sum(position_source1_S2,2)/S2;
            Mean_position_source2_S1 = 0;
            Mean_position_source2_S2 = sum(position_source2_S2,2)/S2;
            
            A_Bias_S1 = 0;
            V_Bias_S1 = 0;
            A_Bias_S2 = 0;
            V_Bias_S2 = 0;
            A_Bias_S1_M = 0;
            A_Bias_S1_SEM = 0;
            A_Bias_S2_M = 0;
            A_Bias_S2_SEM = 0;
            V_Bias_S1_M = 0;
            V_Bias_S1_SEM = 0;
            V_Bias_S2_M = 0;
            V_Bias_S2_SEM = 0;
            
        end
    elseif position_a~=position_v,  % In case of AV distance > 0 the Percentage Bias is expressed with its general mean value, independent from the inferred causal structure
                                    % but also separately for the two conditions: C=1 C=2
        if S2==0,
            
            % case of only C=1 identifications
            A_Bias_TOT=[A_Bias_S1];
            V_Bias_TOT=[V_Bias_S1];
            A_Bias_TOT_deg=[A_Bias_S1_deg];
            V_Bias_TOT_deg=[V_Bias_S1_deg];
            
            Mean_position_source1_S1 = sum(position_source1_S1,2)/S1;
            Mean_position_source1_S2 = 0;
            Mean_position_source2_S1 = sum(position_source2_S1,2)/S1;
            Mean_position_source2_S2 = 0;
            
            A_Bias_S2 = 0;
            V_Bias_S2 = 0;
            
            A_Bias_S1_M = mean(A_Bias_S1);
            A_Bias_S1_SEM = std(A_Bias_S1)/sqrt(S1);
            A_Bias_S2_M = 0;
            A_Bias_S2_SEM = 0;
            V_Bias_S1_M = mean(V_Bias_S1);
            V_Bias_S1_SEM = std(V_Bias_S1)/sqrt(S1);
            V_Bias_S2_M = 0;
            V_Bias_S2_SEM = 0;
            
        elseif S2>0 && S1>0, %
            
            % case with C=1 and C=2 identifications
            A_Bias_TOT=[A_Bias_S1 A_Bias_S2];
            V_Bias_TOT=[V_Bias_S1 V_Bias_S2];
            A_Bias_TOT_deg=[A_Bias_S1_deg A_Bias_S2_deg];
            V_Bias_TOT_deg=[V_Bias_S1_deg V_Bias_S2_deg];
            
            Mean_position_source1_S1 = sum(position_source1_S1,2)/S1;
            Mean_position_source1_S2 = sum(position_source1_S2,2)/S2;
            Mean_position_source2_S1 = sum(position_source2_S1,2)/S1;
            Mean_position_source2_S2 = sum(position_source2_S2,2)/S2;
            
            A_Bias_S1_M = mean(A_Bias_S1);
            A_Bias_S1_SEM = std(A_Bias_S1)/sqrt(S1);
            A_Bias_S2_M = mean(A_Bias_S2);
            A_Bias_S2_SEM = std(A_Bias_S2)/sqrt(S2);
            V_Bias_S1_M = mean(V_Bias_S1);
            V_Bias_S1_SEM = std(V_Bias_S1)/sqrt(S1);
            V_Bias_S2_M = mean(V_Bias_S2);
            V_Bias_S2_SEM = std(V_Bias_S2)/sqrt(S2);
            
        else % S2>0 && S1==0,
            
            % case with only C=2 identifications
            A_Bias_TOT=[A_Bias_S2];
            V_Bias_TOT=[V_Bias_S2];
            A_Bias_TOT_deg=[A_Bias_S2_deg];
            V_Bias_TOT_deg=[V_Bias_S2_deg];
            
            Mean_position_source1_S1 = 0;
            Mean_position_source1_S2 = sum(position_source1_S2,2)/S2;
            Mean_position_source2_S1 = 0;
            Mean_position_source2_S2 = sum(position_source2_S2,2)/S2;
            
            A_Bias_S1_M = 0;
            A_Bias_S1_SEM = 0;
            A_Bias_S2_M = mean(A_Bias_S2);
            A_Bias_S2_SEM = std(A_Bias_S2)/sqrt(S2);
            V_Bias_S1_M = 0;
            V_Bias_S1_SEM = 0;
            V_Bias_S2_M = mean(V_Bias_S2);
            V_Bias_S2_SEM = std(V_Bias_S2)/sqrt(S2);
        end
    end
    
    S_TOT=S1+S2;
    A_Bias_TOT_M = mean(A_Bias_TOT);
    A_Bias_TOT_SEM = std(A_Bias_TOT)/sqrt(S_TOT);
    A_Bias_TOT_deg_M = mean(A_Bias_TOT_deg);
    A_Bias_TOT_deg_SEM = std(A_Bias_TOT_deg)/sqrt(S_TOT);
    
    V_Bias_TOT_M = mean(V_Bias_TOT);
    V_Bias_TOT_SEM = std(V_Bias_TOT)/sqrt(S_TOT);
    V_Bias_TOT_deg_M = mean(V_Bias_TOT_deg);
    V_Bias_TOT_deg_SEM = std(V_Bias_TOT_deg)/sqrt(S_TOT);
    

    eval(['save Data_InputSources_D' num2str(A_Displacement) ' M A V S S1 S2 position_source1_S1 position_source1_S2 position_source2_S1 position_source2_S2 A_Bias_TOT_M A_Bias_TOT_SEM V_Bias_TOT_M V_Bias_TOT_SEM Mean_position_source1_S1 Mean_position_source2_S1 A_Bias_S1 V_Bias_S1 A_Bias_S1_M A_Bias_S1_SEM V_Bias_S1_M V_Bias_S1_SEM Mean_position_source1_S2 Mean_position_source2_S2 A_Bias_S2_M A_Bias_S2_SEM V_Bias_S2_M V_Bias_S2_SEM A_Bias_S2 V_Bias_S2 A_Bias_S2_deg A_Bias_S1_deg V_Bias_S2_deg V_Bias_S1_deg A_Bias_TOT_deg V_Bias_TOT_deg A_Bias_TOT_deg_M V_Bias_TOT_deg_M A_Bias_TOT_deg_SEM V_Bias_TOT_deg_SEM'])
end

