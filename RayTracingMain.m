%% Ray Tracing Engine Parameters    
tic
optimizationMode = 0;   % if 1, then reduces the Rx mesh size to measurement locations only
plotMode = 1;           % Ray Tracing engine plot mode
demoMode = 1;           % Ray Tracing Engine Plot Mode

losFlag = 1; 
reflectionFlag = 1;                     % whether or not calculate First reflections
secondReflectionFlag = 1;
reflectExaggerationFac = 1e0;          % Must be 1 unless to emphasize reflection for demonstration purposes
                                    % whether or not calculate LoS
disableIncidentAngle = 0;               % 1 Disables the incident angle calculation, if disableIncidentAngle= 1
solidIncidentAngle = 45;                % if disableIncidentAngle =1, then assign this which overwrites all the incident angles! This is unnecessary feature
polarizationSwap = 1;               % (See notes in HOW THIS WORSK)  % 1, Applies TE to walls and TM to ceiling. 0, applies TM to the walls and TE to the ceiling

imageRSSIScale = 5;         % increase this if number of meshes nodes are small
grayScaleImage = 0;

freq = 59.6e9;  % frequency in hertz
% freq = 2.45e9;  % frequency in hertz
lightVel = 3e8;
lambda = lightVel./freq;
vel = freq*lambda;
refDistance = 1;            % Reference distance from Tx
FPSLRefLoss = 0;

antennaGainRes = 40;
antennaEffiLoss = -11.5;         % dB antenna efficiency, missmatch and loss all together

ceilingEnable = 0; % Allowing to define ceiling and floor
groundLevel = 0;
ceilingLevel = 4060*1000  % Height of the ceiling

% For LOS dis = length of ceiling
timeLOS = ceilingLevel/vel;


mesh_.xNodeNum = 100;   % Keep the x and y mesh size the same, increase the size for better resolution and especially if you're increasing the frequency
mesh_.yNodeNum = 100;
mesh_.zNodeNum = 1;

%% Antenna Gain pattern calculation

[TxAntennaGainAE] = AntennaTemp (antennaGainRes,demoMode) + antennaEffiLoss;  % TxAntennaGainAE needs to be in dB
RxAntennaGainAE = TxAntennaGainAE;

Tx.xyz = [ 200, 200, 1]; % Location of the transmitter (s)
Tx.power =  67; % power of the transmitter dB(m) % start from 47 at 2.4 GHz
% start from 2 at 60 GHz and decrease by 20 for duct
% start from 2 at 60 GHz and decrease by 2 for duct

% Defining the boundary of the analysis (something like a boundary condition) 
boundary = [0,406
            0,406
            0,4060];    

%% CLOCK WISE WALL DEFINITION        
[wallxyz1, wallxyz2, wallxyz3, wallxyz4,wallX,wallY,wallZ] = CSV23D_V1(demoMode,groundLevel,ceilingLevel,Tx.xyz);

wall.xyz1 = wallxyz1;
wall.xyz2 = wallxyz2;
wall.xyz3 = wallxyz3;
wall.xyz4 = wallxyz4;

wall.X = wallX;
wall.Y = wallY;
wall.Z = wallZ;

if ceilingEnable == 1
    ceillFloor.xyz1 = [0,0,ceilingLevel
                       0,0,groundLevel
                        ];

    ceillFloor.xyz2 = [0,406,ceilingLevel
                       0,406,groundLevel
                        ];

    ceillFloor.xyz3 = [406,406,ceilingLevel
                        406,406,groundLevel
                        ];

    ceillFloor.xyz4 = [406,0,ceilingLevel
                        406,0,groundLevel
                        ];
else
    ceillFloor.xyz1 = [];
    ceillFloor.xyz2 = [];
    ceillFloor.xyz3 = [];
    ceillFloor.xyz4 = [];
end
wall.relativePerm = 6*ones(size(wall.xyz1,1)+size(ceillFloor.xyz1,1),1);
%% Adding Ceillilng and Floor to the structure
for i = 1:size(ceillFloor.xyz1,1)
    wall.xyz1 = [wall.xyz1;ceillFloor.xyz1(i,:)];
    wall.xyz2 = [wall.xyz2;ceillFloor.xyz2(i,:)];
    wall.xyz3 = [wall.xyz3;ceillFloor.xyz3(i,:)];
    wall.xyz4 = [wall.xyz4;ceillFloor.xyz4(i,:)];
end
RayTracingEng_V01

toc

L = [1 2 4 8];
RwP = [-34.1212 -37.7700 -42.4794 -47.2173];
RP = [-36.8592 -40.5080 -45.2173 -50.9552];
RMP = [-94.0744 -97.6685 -101.3379 -105.0773];
PP1 = [-82 -83 -82  -83];
PP2 = [-80 -80 -80 -80];
fspl = [-54 -60 -66 -72];
figure 
plot(L, RwP, 'r', L, RP, 'y',L, PP1, 'g',L, PP2, 'b',L, fspl, 'k');