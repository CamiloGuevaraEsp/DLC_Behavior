
close all; 
clear all;

%once you have analyzed a video using DLC, export that file as a csv from
%the command line using "deeplabcut.analyze_videos_converth5_to_csv(h5path,
%videotype='.avi')"

%%%%%%%%%%%%%START parameters to modify:
DLC_csv = 'fc2_save_2022-01-17-163110-_2Xspeed_frameDLC_resnet50_flyonaballRevisedDec2shuffle1_248000.csv';
file_path = 'A:\Theresa\Downloads';
outname = strrep( DLC_csv, '_2Xspeed_frameDLC_resnet50_flyonaballRevisedDec2shuffle1_248000', 'DLCmvmt');
%%%%%%%%%%%%%END parameters to modify:


cd(file_path);
%varNames = {'FrameNumber', 'FrontNearFoot', 'FrontFarFoot', 'MiddleNearFoot', 'MiddleFarFoot', ...
%    'BackNearFoot', 'BackFarFoot', 'tipofabdomen', 'tipofProboscis', 'CenterofFecalMatter', 'CenterofEgg' };

M = readtable(DLC_csv);
%M.Properties.VariableNames = varNames; % names of columns
M = table2array(M);

FrameNumber = [];
FrontNearFoot = [];
FrontFarFoot = [];
MiddleNearFoot = [];
MiddleFarFoot = [];
BackNearFoot = [];
BackFarFoot = [];
tipofabdomen = [];
tipofProboscis = [];
CenterofFecalMatter = [];
CenterofEgg = [];
mvmtscore = []; %will create movement score by adding together the euclidian 
% distances of the bottom of each leg and the abdomen score... these seem
% to have the least amount of uncertainty from the model and should be
% enough information to score movement

for i = 2:size(M,1)
    FrontNearFoot_currentX = M(i,2);
    FrontNearFoot_currentY = M(i,3);
    FrontNearFoot_currentP = M(i,4);
    
    FrontFarFoot_currentX = M(i,5);
    FrontFarFoot_currentY = M(i,6);
    FrontFarFoot_currentP = M(i,7);

    MiddleNearFoot_currentX = M(i,8);
    MiddleNearFoot_currentY = M(i,9);
    MiddleNearFoot_currentP = M(i,10);
    
    MiddleFarFoot_currentX = M(i,11);
    MiddleFarFoot_currentY = M(i,12);
    MiddleFarFoot_currentP = M(i,13);
    
    BackNearFoot_currentX = M(i,14);
    BackNearFoot_currentY = M(i,15);
    BackNearFoot_currentP = M(i,16);

    BackFarFoot_currentX = M(i,17);
    BackFarFoot_currentY = M(i,18);
    BackFarFoot_currentP = M(i,19);
    
    tipofabdomen_currentX = M(i,20);
    tipofabdomen_currentY = M(i,21);
    tipofabdomen_currentP =  M(i,22);
    
    tipofProboscis_currentX = M(i,23);
    tipofProboscis_currentY = M(i,24);
    tipofProboscis_currentP = M(i,25);
    
    CenterofFecalMatter_currentX = M(i,26);
    CenterofFecalMatter_currentY = M(i,27);
    CenterofFecalMatter_currentP = M(i,28);

    CenterofEgg_currentX = M(i,29);
    CenterofEgg_currentY = M(i,30);
    CenterofEgg_currentP = M(i,31);

    FrontNearFoot_prevX = M(i-1,2);
    FrontNearFoot_prevY = M(i-1,3);
    FrontNearFoot_prevP = M(i-1,4);
    
    FrontFarFoot_prevX = M(i-1,5);
    FrontFarFoot_prevY = M(i-1,6);
    FrontFarFoot_prevP = M(i-1,7);

    MiddleNearFoot_prevX = M(i-1,8);
    MiddleNearFoot_prevY = M(i-1,9);
    MiddleNearFoot_prevP = M(i-1,10);
    
    MiddleFarFoot_prevX = M(i-1,11);
    MiddleFarFoot_prevY = M(i-1,12);
    MiddleFarFoot_prevP = M(i-1,13);
    
    BackNearFoot_prevX = M(i-1,14);
    BackNearFoot_prevY = M(i-1,15);
    BackNearFoot_prevP = M(i-1,16);

    BackFarFoot_prevX = M(i-1,17);
    BackFarFoot_prevY = M(i-1,18);
    BackFarFoot_prevP = M(i-1,19);
    
    tipofabdomen_prevX = M(i-1,20);
    tipofabdomen_prevY = M(i-1,21);
    tipofabdomen_prevP =  M(i-1,22);
    
    tipofProboscis_prevX = M(i-1,23);
    tipofProboscis_prevY = M(i-1,24);
    tipofProboscis_prevP = M(i-1,25);
    
    CenterofFecalMatter_prevX = M(i-1,26);
    CenterofFecalMatter_prevY = M(i-1,27);
    CenterofFecalMatter_prevP = M(i-1,28);

    CenterofEgg_prevX = M(i-1,29);
    CenterofEgg_prevY = M(i-1,30);
    CenterofEgg_prevP = M(i-1,31);
    
    FrontNearFootCoord_1=[FrontNearFoot_prevX,FrontNearFoot_prevY];  % The first coordinate
    FrontNearFootCoord_2=[FrontNearFoot_currentX,FrontNearFoot_currentY];  % The second coordinate
    FrontNearFootpair=[FrontNearFootCoord_1; FrontNearFootCoord_2];
    FrontNearFootdistance=pdist(FrontNearFootpair,'euclidean');

    FrontFarFootCoord_1=[FrontFarFoot_prevX,FrontFarFoot_prevY];  % The first coordinate
    FrontFarFootCoord_2=[FrontFarFoot_currentX,FrontFarFoot_currentY];  % The second coordinate
    FrontFarFootpair=[FrontFarFootCoord_1; FrontFarFootCoord_2];
    FrontFarFootdistance=pdist(FrontFarFootpair,'euclidean');

    MiddleNearFootCoord_1=[MiddleNearFoot_prevX,MiddleNearFoot_prevY];  % The first coordinate
    MiddleNearFootCoord_2=[MiddleNearFoot_currentX,MiddleNearFoot_currentY];  % The second coordinate
    MiddleNearFootpair=[MiddleNearFootCoord_1; MiddleNearFootCoord_2];
    MiddleNearFootdistance=pdist(MiddleNearFootpair,'euclidean');

    MiddleFarFootCoord_1=[MiddleFarFoot_prevX,MiddleFarFoot_prevY];  % The first coordinate
    MiddleFarFootCoord_2=[MiddleFarFoot_currentX,MiddleFarFoot_currentY];  % The second coordinate
    MiddleFarFootpair=[MiddleFarFootCoord_1; MiddleFarFootCoord_2];
    MiddleFarFootdistance=pdist(MiddleFarFootpair,'euclidean');

    BackNearFootCoord_1=[BackNearFoot_prevX,BackNearFoot_prevY];  % The first coordinate
    BackNearFootCoord_2=[BackNearFoot_currentX,BackNearFoot_currentY];  % The second coordinate
    BackNearFootpair=[BackNearFootCoord_1; BackNearFootCoord_2];
    BackNearFootdistance=pdist(BackNearFootpair,'euclidean');

    BackFarFootCoord_1=[BackFarFoot_prevX,BackFarFoot_prevY];  % The first coordinate
    BackFarFootCoord_2=[BackFarFoot_currentX,BackFarFoot_currentY];  % The second coordinate
    BackFarFootpair=[BackFarFootCoord_1; BackFarFootCoord_2];
    BackFarFootdistance=pdist(BackFarFootpair,'euclidean');

    tipofabdomenCoord_1=[tipofabdomen_prevX,tipofabdomen_prevY];  % The first coordinate
    tipofabdomenCoord_2=[tipofabdomen_currentX,tipofabdomen_currentY];  % The second coordinate
    tipofabdomenpair=[tipofabdomenCoord_1; tipofabdomenCoord_2];
    tipofabdomendistance=pdist(tipofabdomenpair,'euclidean');

    tipofProboscisCoord_1=[tipofProboscis_prevX,tipofProboscis_prevY];  % The first coordinate
    tipofProboscisCoord_2=[tipofProboscis_currentX,tipofProboscis_currentY];  % The second coordinate
    tipofProboscispair=[tipofProboscisCoord_1; tipofProboscisCoord_2];
    tipofProboscisdistance=pdist(tipofProboscispair,'euclidean');

    CenterofFecalMatterCoord_1=[CenterofFecalMatter_prevX,CenterofFecalMatter_prevY];  % The first coordinate
    CenterofFecalMatterCoord_2=[CenterofFecalMatter_currentX,CenterofFecalMatter_currentY];  % The second coordinate
    CenterofFecalMatterpair=[CenterofFecalMatterCoord_1; CenterofFecalMatterCoord_2];
    CenterofFecalMatterdistance=pdist(CenterofFecalMatterpair,'euclidean');

    CenterofEggCoord_1=[CenterofEgg_prevX,CenterofEgg_prevY];  % The first coordinate
    CenterofEggCoord_2=[CenterofEgg_currentX,CenterofEgg_currentY];  % The second coordinate
    CenterofEggpair=[CenterofEggCoord_1; CenterofEggCoord_2];
    CenterofEggdistance=pdist(CenterofEggpair,'euclidean');
    
   
    movement = tipofabdomendistance + BackFarFootdistance + BackNearFootdistance + MiddleFarFootdistance + MiddleNearFootdistance...
        + FrontFarFootdistance + FrontNearFootdistance;

    FrameNumber = [FrameNumber; M(i,1)];
    FrontNearFoot = [FrontNearFoot; FrontNearFootdistance];
    FrontFarFoot = [FrontFarFoot; FrontFarFootdistance];
    MiddleNearFoot = [MiddleNearFoot; MiddleNearFootdistance];
    BackNearFoot = [BackNearFoot; BackNearFootdistance];
    BackFarFoot = [BackFarFoot; BackFarFootdistance];
    tipofabdomen = [tipofabdomen; tipofabdomendistance];
    tipofProboscis = [tipofProboscis; tipofProboscisdistance];
    CenterofFecalMatter = [CenterofFecalMatter; CenterofFecalMatterdistance];
    CenterofEgg = [CenterofEgg; CenterofEggdistance];
    mvmtscore = [mvmtscore; movement];

end

varNames2 = {'FrameNumber', 'FrontNearFoot', 'FrontFarFoot', 'MiddleNearFoot', ...
     'BackNearFoot', 'BackFarFoot', 'tipofabdomen', 'tipofProboscis'...
    'CenterofFecalMatter', 'CenterofEgg', 'Movement'};

DistanceOutput = table(FrameNumber,FrontNearFoot,FrontFarFoot,MiddleNearFoot,...
    BackNearFoot,BackFarFoot,tipofabdomen,tipofProboscis,CenterofFecalMatter,CenterofEgg, mvmtscore);

DistanceOutput.Properties.VariableNames = varNames2; % names of columns
expdate1 = extractAfter(DLC_csv,"fc2_save_");
expdate2 = extractBefore(DLC_csv,"-_");

writetable(DistanceOutput, outname)
plot(mvmtscore)
