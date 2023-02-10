
%%%%adapted from Cynthia's code to allow for the quantification of signal from
%%%%individual slices of data in addition to the average of the whole z stack.
%%%%also aligns the DLC behavioral data rather than data resulting from the 
%%%%frame subtraction method. Finally, it includes a step to subtract the
%%%%background of Channel 2 in addition to normalizing to the Channel 1 
%%%%(i.e, the tdtomato signal). This code will always subtract the Ch2
%%%%background first, and if the variable readChannel1 = 1 below, it will 
%%%%also normalize the background subtracted data to the Ch1 data

%before running this code, you must have: 
%1. run writewholevidmp4.m
%2. analyze video using DLC 
%3. convert DLC h5 file to csv
%4. run code called calcDLCbodypartMvmt.m to create the mvmtfile 
%5. run getpixeldat_singlevid.m to create diffArray file
%5. run both write2PhotonProJmask_slices.m and write2PhotonProJmask.m, 
% write2ProtonProJmask code performance, you may have also run uiDrawMasks.m 
%6. put the output of all of these codes into one centralized folder

close all;
%clear variables;

%%%%%%%%%%%%%START parameters to modify:
rootdir = 'F:\011722_processed121422';
Tseriesdir = 'F:\';
savedir = 'F:\011722_processed121422';
Tseriesroot = 'TSeries-01172022-1611-688';
Mvmtfilename = 'fc2_save_2022-01-17-163110-DLCmvmt.csv';
vidRootName = 'fc2_save_2022-01-17-163110-';
diffArray = 'fc2_save_2022-01-17-163110-_2Xspeed_frame.mat';
TProjMask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__userDrawnMask';
laserThresh = 20000; %This value may have to be adjusted as appropriate. 
%In general, it has been working well
cycleStart = 1;
cycleList = 1 ; %If Piezo is not used, each stack is a separate cycle. %In 
% high speed T-series (Piezo is used, typical for high frequency imaging), 
% just one cycle for the entire series.
slicesPerStack = 7;
writeVid = 0;
upperStretchLim = 0.1; %this variable is related to writing the video
readChannel1 = 0; %%% if set to 0 choosing to only subtract Ch2 background ; 
% if set to 1, choosing to subtract channel2 background and normalize to Ch1 data
IntOutName = Tseriesroot +  "_fullMovementAndBrainSignal.mat";
medfiltSize = 1; %% this is related to the image filtering
recordedFPS = 30; %frames per second of behavioral video
FalseStart = 0 ; %if the laser was turned on and the laser flashed but you 
% had to restart the tseries - set false start to 1 and count the number of
% flashes 
numlaserflashes = 1;% if FalseStart = 1, input the # of flashes before the 
% real experiment began
Slice1Mask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__slice1_userDrawnMask.mat';
Slice2Mask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__slice2_userDrawnMask.mat';
Slice3Mask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__slice3_userDrawnMask.mat';
Slice4Mask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__slice4_userDrawnMask.mat';
Slice5Mask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__slice5_userDrawnMask.mat';
Slice6Mask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__slice6_userDrawnMask.mat';
Slice7Mask = 'TSeries-01172022-1611-688_Cycle00001_Ch2__slice7_userDrawnMask.mat';
imagedimensions = 256*256; %this could change if the experiment is run at a 
% resolution other than 256 x 256
backgrounddimension = 25; %this is the size of one edge of the square that will be pulled for background %subtraction
backgroundarea = backgrounddimension*backgrounddimension;
% the top left corner that will be used to subtract the background
%%%%%%%%%%%%%END parameters to modify:


% plot the movement data
cd(rootdir);
mvmtdata = readmatrix(Mvmtfilename);
mvmtscorearray = mvmtdata(:,11);
avgmvmt = mean(mvmtscorearray);
baseline = 2*avgmvmt;
BLadjustedmvmtscore = mvmtscorearray-baseline;

BLadjustedfullMovement = NaN(length(mvmtdata),2);
BLadjustedfullMovement(:,1) = mvmtdata(:,1);
BLadjustedfullMovement(:,2) = BLadjustedmvmtscore;

lastFrameNum = length(mvmtdata);
framesPerMin = 60*recordedFPS;

%plot(fullMovement(:,2));


cycleStartString = ['Cycle' num2str(cycleStart,'%05.0f')];

TImgRoot = [Tseriesroot '_' cycleStartString '_Ch2_'];%000001.ome


if(writeVid)
    vidObj = VideoWriter([TImgRoot '_upperStretchLim ' num2str(upperStretchLim) '.avi']);
end

%load in the masks and calculate the size of the mask (i.e. # pixels
%included in mask)


A = load(TProjMask);
tProjMaskMat = A.bwFrame;
OneDmask = reshape(tProjMaskMat.',1,[]);
pixelsinmask = sum(OneDmask);

MaskbySliceArray = NaN(length(OneDmask), slicesPerStack);  %later we will 
% use this to apply the mask to all 7 slices of data at once

B = load(Slice1Mask);
Slice1MaskMat = B.bwFrame1;
S1OneDmask = reshape(Slice1MaskMat.',1,[]);
pixelsinmask1 = sum(S1OneDmask);
MaskbySliceArray(:,1) = S1OneDmask;

C = load(Slice2Mask);
Slice2MaskMat = C.bwFrame2;
S2OneDmask = reshape(Slice2MaskMat.',1,[]);
pixelsinmask2 = sum(S2OneDmask);
MaskbySliceArray(:,2) = S2OneDmask;

D = load(Slice3Mask);
Slice3MaskMat = D.bwFrame3;
S3OneDmask = reshape(Slice3MaskMat.',1,[]);
pixelsinmask3 = sum(S3OneDmask);
MaskbySliceArray(:,3) = S3OneDmask;

E = load(Slice4Mask);
Slice4MaskMat = E.bwFrame4;
S4OneDmask = reshape(Slice4MaskMat.',1,[]);
pixelsinmask4 = sum(S4OneDmask);
MaskbySliceArray(:,4) = S4OneDmask;

F = load(Slice5Mask);
Slice5MaskMat = F.bwFrame5;
S5OneDmask = reshape(Slice5MaskMat.',1,[]);
pixelsinmask5 = sum(S5OneDmask);
MaskbySliceArray(:,5) = S5OneDmask;

G = load(Slice6Mask);
Slice6MaskMat = G.bwFrame6;
S6OneDmask = reshape(Slice6MaskMat.',1,[]);
pixelsinmask6 = sum(S6OneDmask);
MaskbySliceArray(:,6) = S6OneDmask;

H = load(Slice7Mask);
Slice7MaskMat = H.bwFrame7;
S7OneDmask = reshape(Slice7MaskMat.',1,[]);
pixelsinmask7 = sum(S7OneDmask);
MaskbySliceArray(:,7) = S7OneDmask;

% I = load(Slice8Mask);
% Slice8MaskMat = I.bwFrame8;
% S8OneDmask = reshape(Slice8MaskMat.',1,[]);
% pixelsinmask8 = sum(S8OneDmask);
% MaskbySliceArray(:,8) = S8OneDmask;

sliceMaskPxCount = [pixelsinmask1, pixelsinmask2, pixelsinmask3, pixelsinmask4...
pixelsinmask5, pixelsinmask6, pixelsinmask7];%, pixelsinmask8]; %later we will use this array to 
%normalize the slice data by the number of pixels in their respective masks
%all at once

if(writeVid)
    open(vidObj);
end

disp('Extracting Raw Neural Data and Normalizing to size of the mask')

for ci = 1:numel(cycleList)
    cd([Tseriesdir Tseriesroot]);
    cycleText = ['Cycle' num2str(cycleList(ci),'%05.0f')];
    imgRootForCycle = strrep(TImgRoot,cycleStartString,cycleText);
    tiffList = dir([imgRootForCycle '*.ome.tif']);
    TotalStacks = length(tiffList)/slicesPerStack;
    IntTotalStacks = floor(TotalStacks);
    IntOutName = strrep(TProjMask, '.mat',['_' cycleText '_' 'Intensities.mat']);
    maxFrame2read = numel(tiffList);
    stackNum = 0;
    %initialize arrays to hold intensity data averaged across all slices
    AVGIntData = NaN(IntTotalStacks, 1);
    Ch1AVGIntData = NaN(IntTotalStacks, 1);
    Ch2bkgndAVGIntData = NaN(IntTotalStacks, 1);
    %initialize array to hold the intensity data for each slice across all stacks
    pixelsbyslicearray = NaN(IntTotalStacks,slicesPerStack);
    Ch1pixelsbyslicearray = NaN(IntTotalStacks,slicesPerStack);
    Ch2bkgndpixelsbyslicearray = NaN(IntTotalStacks,slicesPerStack);
    %before comparing the raw pixels to the background pixels, both should
    %be divided by the size of the mask/background image... 1

    for ti = 1:maxFrame2read
        %if mod(ti,1000) == 0
            %display(ti)
        %end
   %for ti = 1:slicesPerStack
        slicenum = mod(ti,slicesPerStack);
        if slicenum == 0
            %slicenum = 8;
            slicenum = 7;
        end
        
%         if(ti<=200 || mod(ti,100)==0)
%             disp(ti);
%             if(ti==200)
%                 disp('Have demonstrated first 200 frames are read. Now outputting once every 100 frames.');
%             end
%         end
        
        imgName = [imgRootForCycle num2str(ti,'%06.0f') '.ome.tif'];

        if(exist(imgName,'file'))
      
            stackIndex = mod(ti,slicesPerStack);
            try
                rawImg = imread(imgName); %This rawImg contains the data that we probably want to extract data from (GCamp6m)
            catch
                rawImg(:) = NaN; %(size(rawImg,1),size(rawImg,2));
            end

            if readChannel1 
          
                ctrlImgName = strrep(imgName,'Ch2','Ch1');
                if(exist(ctrlImgName,'file'))
                    try
                        ctrlImg = imread(ctrlImgName);
                    catch
                        ctrlImg(:) = NaN;
                        display(['Could not read ' ctrlImgName]);
                    end
                else
                    ctrlImg(:) = NaN;
                    disp('Could not find Channel 1 data')
                    return
                end

                Ch2bkgndImg = imcrop(rawImg,[0 0 backgrounddimension backgrounddimension]);

            else
                Ch2bkgndImg = imcrop(rawImg,[0 0 backgrounddimension backgrounddimension]);
                %create bkgnd image with pixels from the top left corner of
                %the raw img... this should be a 25 x 25 cropped area from
                %the top right corner of the raw img

            end

            if slicenum==1
                stackNum = stackNum+1;

                %initalize Gcamp data arrays
                FiltImageDataInStack = NaN(size(rawImg,1)*size(rawImg,2),slicesPerStack);
                noFiltimageDataInStack = NaN(size(rawImg,1)*size(rawImg,2),slicesPerStack);

                if(ti~=1 && writeVid)
                    
                    I = getframe(h);
                    writeVideo(vidObj,I); 
                    close(figure(1));
               
                end
                 
                if readChannel1 %initialize ctrl arrays 
                 
                    Ch1ctrlFiltDataInStack = NaN(size(ctrlImg,1)*size(ctrlImg,2),slicesPerStack);
                    Ch1ctrlNoFiltDataInStack = NaN(size(ctrlImg,1)*size(ctrlImg,2),slicesPerStack);

                    Ch2bkgndFiltDataInStack = NaN(size(Ch2bkgndImg,1)*size(Ch2bkgndImg,2),slicesPerStack); %initialize Ch2 background arrays
                    Ch2bkgndNoFiltDataInStack = NaN(size(Ch2bkgndImg,1)*size(Ch2bkgndImg,2),slicesPerStack);
               
                else

                    Ch2bkgndFiltDataInStack = NaN(size(Ch2bkgndImg,1)*size(Ch2bkgndImg,2),slicesPerStack); %initialize Ch2 background arrays
                    Ch2bkgndNoFiltDataInStack = NaN(size(Ch2bkgndImg,1)*size(Ch2bkgndImg,2),slicesPerStack);
               
                end


            elseif slicenum==0 && writeVid
                sumProjectionVector = nansum(FiltImageDataInStack,2);
                sumProjection = reshape(sumProjectionVector,size(rawImg,1),size(rawImg,2));
                h = figure(1);
                imagesc(sumProjection,[0 4000]); 
            end

            if(slicenum~=slicesPerStack) %Just save the data into the stack.
  
                noFiltimageDataInStack(:,slicenum) = double(rawImg(:)); 
                FiltFrame = double(medfilt2(uint16(reshape(rawImg,size(rawImg,1),size(rawImg,2))),[medfiltSize medfiltSize]));
                FiltImageDataInStack(:,slicenum) = double(FiltFrame(:));  

                if readChannel1 %add Ch1 data to control (Ch1) stack array
            
                    Ch1ctrlNoFiltDataInStack(:,slicenum) = double(ctrlImg(:)); 
                    Ch1ctrlFiltFrame = double(medfilt2(uint16(reshape(ctrlImg,size(ctrlImg,1),size(ctrlImg,2))),[medfiltSize medfiltSize]));
                    Ch1ctrlFiltDataInStack(:,slicenum) = double(Ch1ctrlFiltFrame(:));
                
                    Ch2bkgndNoFiltDataInStack(:,slicenum) = double(Ch2bkgndImg(:));
                    Ch2bkgndFiltFrame = double(medfilt2(uint16(reshape(Ch2bkgndImg,size(Ch2bkgndImg,1),size(Ch2bkgndImg,2))),[medfiltSize medfiltSize]));
                    Ch2bkgndFiltDataInStack(:,slicenum) = double(Ch2bkgndFiltFrame(:));
                    
                else %add Ch2 background to control array
                    Ch2bkgndNoFiltDataInStack(:,slicenum) = double(Ch2bkgndImg(:));
                    Ch2bkgndFiltFrame = double(medfilt2(uint16(reshape(Ch2bkgndImg,size(Ch2bkgndImg,1),size(Ch2bkgndImg,2))),[medfiltSize medfiltSize]));
                    Ch2bkgndFiltDataInStack(:,slicenum) = double(Ch2bkgndFiltFrame(:));
                    
                end

            else %we are at the end of a stack and we want to add the last row of data from this stack to the array
           
                noFiltimageDataInStack(:,slicenum) = double(rawImg(:));
%                 FiltFrame = double(medfilt2(uint16(reshape(rawImg,size(rawImg,1),size(rawImg,2))),[medfiltSize medfiltSize]));
%                 FiltImageDataInStack(:,slicenum) = double(FiltFrame(:));
%                 multiply the data by the binary mask. everything outside 
%                 mask becomes 0 
%                 maskedDataInStack = FiltImageDataInStack.*MaskbySliceArray; 


                % for the whole dFB average data....Find the mean image for the whole stack.
                meanImg = mean(noFiltimageDataInStack,2);
                %filter that image
                meanFiltFrame = double(medfilt2(meanImg, [medfiltSize medfiltSize]));
                %multiply the data by the binary mask. everything outside 
                % mask becomes 0 
                meanSignal = meanFiltFrame.*OneDmask.'; 
                %get the sum of pixels in each slice once the mask is applied
                SumofpixelAVGincurrentstack = sum(meanSignal);
                %normalize the total intensity by the number of pixels
                %in the mask
                CurrAVGIntNormtoPxinMask = SumofpixelAVGincurrentstack/pixelsinmask;
                AVGIntData(stackNum)= CurrAVGIntNormtoPxinMask;

                  %for the slice data
                %each column of the noFiltimg data by stack is the data for
                %that slice..we do not want to avg it together
                %instead immediately multiply by the slice mask array so that 
                % for each slice, everything outside its mask is 0...
                maskedSliceData = noFiltimageDataInStack.*MaskbySliceArray;
                SumofPixelAvgBySliceInCurrStack = sum(maskedSliceData,1);
                 %normalize the total intensity by the number of pixels
                %in each mask
                CurrSliceIntNormtoPxinMask = SumofPixelAvgBySliceInCurrStack./sliceMaskPxCount;
                pixelsbyslicearray(stackNum,:) = CurrSliceIntNormtoPxinMask;
            end

            if readChannel1
        
                %add Ch1 control data to stack and save that data to
                %control arrays
                Ch1ctrlNoFiltDataInStack(:,slicenum) = double(ctrlImg(:)); 
                Ch1ctrlFiltFrame = double(medfilt2(uint16(reshape(ctrlImg,size(ctrlImg,1),size(ctrlImg,2))),[medfiltSize medfiltSize]));
                Ch1ctrlFiltDataInStack(:,slicenum) = double(Ch1ctrlFiltFrame(:));
                
               
                %for AVG Control data
                %Find the mean.
                CtrlmeanImg = mean(Ch1ctrlNoFiltDataInStack,2); 
                %filter that image
                Ch1CtrlmeanFiltFrame = double(medfilt2(CtrlmeanImg, [medfiltSize medfiltSize]));
                %multiply the data by the binary mask. everything outside 
                % mask becomes 0 
                Ch1meanSignal = Ch1CtrlmeanFiltFrame.*OneDmask.'; 
                %get the sum of pixels in each slice once the mask is applied
                Ch1SumofpixelAVGincurrentstack = sum(Ch1meanSignal);
                %normalize the total intensity by the number of pixels
                %in the mask
                Ch1CurrStackIntNormtoPxinMask = Ch1SumofpixelAVGincurrentstack/pixelsinmask;
                %append the normalized data to an array
                Ch1AVGIntData(stackNum)= Ch1CurrStackIntNormtoPxinMask;
                
                
                %for slice data
                %each column of the noFiltimg data by stack is the data for
                %that slice..we do not want to avg it together
                %instead immediately multiply by the slice mask array so that 
                % for each slice, everything outside its mask is 0...
                maskedCh1CtrlSliceData = Ch1CtrlmeanFiltFrame.*MaskbySliceArray;
                SumofPixelAvgCh1CtrlBySliceInCurrStack = sum(maskedCh1CtrlSliceData,1);
                %normalize the total intensity by the number of pixels
                %in each mask
                CurrStackIntCh1CtrlNormtoPxinMask = SumofPixelAvgCh1CtrlBySliceInCurrStack./sliceMaskPxCount;
                Ch1pixelsbyslicearray(stackNum,:) = CurrStackIntCh1CtrlNormtoPxinMask;
                
                 %add last slice in data to the stack array
                Ch2bkgndNoFiltDataInStack(:,slicenum) = double(Ch2bkgndImg(:));
%                 Ch2bkgndFiltFrame = double(medfilt2(uint16(reshape(Ch2bkgndImg,size(Ch2bkgndImg,1),size(Ch2bkgndImg,2))),[medfiltSize medfiltSize]));
%                 Ch2bkgndFiltDataInStack(:,slicenum) = double(Ch2bkgndFiltFrame(:));
               
                %for AVG Control data
                %Find the mean image of the Ch2 data.
                CtrlmeanImg = mean(Ch2bkgndNoFiltDataInStack,2); 
                %filter that image
                Ch2bkgndmeanFiltFrame = double(medfilt2(CtrlmeanImg, [medfiltSize medfiltSize]));
    
                %dont need to use the binary mask here... using the
                %background signal per pixel to subtract the background

                %get the sum of pixels 
                Ch2bkgndSumofpixelAVGincurrentstack = sum(Ch2bkgndmeanFiltFrame);
                %norm the amount of background to the size of the background image
                Ch2bkgndCurrSliceIntNormtoPxinBkgnd = Ch2bkgndSumofpixelAVGincurrentstack/backgroundarea;
                % backgroundarea = to the size of the cropped background image
                %append the normalized data to an array
                Ch2bkgndAVGIntData(stackNum)= Ch2bkgndCurrSliceIntNormtoPxinBkgnd;

                %for slice data
                %each column of the noFiltimg data by stack is the data for
                %that slice..we do not want to avg it together
                
                %also dont need to use the binary mask here... using the
                %background signal per pixel to subtract the background
                
                Ch2bkgndSumofpixelbySliceincurrentstack = sum(Ch2bkgndmeanFiltFrame,1);
                %normalize the total intensity by the number of pixels
                %in each mask
                CurrStackIntCh2CtrlNormtoPxinMask = Ch2bkgndSumofpixelbySliceincurrentstack./sliceMaskPxCount;
                Ch2bkgndpixelsbyslicearray(stackNum,:) = CurrStackIntCh2CtrlNormtoPxinMask;

            else 
                %add last slice in data to the stack array
                Ch2bkgndNoFiltDataInStack(:,slicenum) = double(Ch2bkgndImg(:));
%                 Ch2bkgndFiltFrame = double(medfilt2(uint16(reshape(Ch2bkgndImg,size(Ch2bkgndImg,1),size(Ch2bkgndImg,2))),[medfiltSize medfiltSize]));
%                 Ch2bkgndFiltDataInStack(:,slicenum) = double(Ch2bkgndFiltFrame(:));
               
                %for AVG Control data
                %Find the mean image of the Ch2 data.
                CtrlmeanImg = mean(Ch2bkgndNoFiltDataInStack,2); 
                %filter that image
                Ch2bkgndmeanFiltFrame = double(medfilt2(CtrlmeanImg, [medfiltSize medfiltSize]));
    
                %dont need to use the binary mask here... using the
                %background signal per pixel to subtract the background

                %get the sum of pixels 
                Ch2bkgndSumofpixelAVGincurrentstack = sum(Ch2bkgndmeanFiltFrame);
                %norm the amount of background to the size of the background image
                Ch2bkgndCurrSliceIntNormtoPxinBkgnd = Ch2bkgndSumofpixelAVGincurrentstack/backgroundarea;
                % backgroundarea = to the size of the cropped background image
                %append the normalized data to an array
                Ch2bkgndAVGIntData(stackNum)= Ch2bkgndCurrSliceIntNormtoPxinBkgnd;

                %for slice data
                %each column of the noFiltimg data by stack is the data for
                %that slice..we do not want to avg it together
                
                %also dont need to use the binary mask here... using the
                %background signal per pixel to subtract the background
                
                Ch2bkgndSumofpixelbySliceincurrentstack = sum(Ch2bkgndmeanFiltFrame,1);
                %normalize the total intensity by the number of pixels
                %in each mask
                CurrStackIntCh2CtrlNormtoPxinMask = Ch2bkgndSumofpixelbySliceincurrentstack./sliceMaskPxCount;
                Ch2bkgndpixelsbyslicearray(stackNum,:) = CurrStackIntCh2CtrlNormtoPxinMask;
                
            end
      
        else
            disp('missing image:')
            disp(imgName)
            
        end
  
    if(writeVid)
        close(vidObj)
    end

   end
end
%% 


%initiate the array that will hold and align the behavioral and neural data
fullMovementAndBrainSignal = NaN(numel(BLadjustedfullMovement(:,2)),3);
%creating an array of NaNs that is the size of the movement file by 3
fullMovementAndBrainSignal(:,1) = BLadjustedfullMovement(:,1);%fill the first column of
% fullMovement with frame number
fullMovementAndBrainSignal(:,2) = BLadjustedfullMovement(:,2);%fill the 2nd column of
% fullMovement with movement data

fullMovementAndBrainSignalSlicedata = NaN(numel(BLadjustedfullMovement(:,2)),slicesPerStack+2);
fullMovementAndBrainSignalSlicedata(:,1) = BLadjustedfullMovement(:,1);
fullMovementAndBrainSignalSlicedata(:,2) = BLadjustedfullMovement(:,2);

%get the pixel intensity data for each slice or stack and normalize it:

if readChannel1
    disp('Subtracting background from Gcamp signal AND Normalizing Gcamp signal to Channel 1 data')
    
    slicedata = pixelsbyslicearray-Ch2bkgndpixelsbyslicearray;
    Ch2norm_AVG_Gcampdata = AVGIntData-Ch2bkgndAVGIntData;
    
    slicedata = pixelsbyslicearray./Ch1pixelsbyslicearray;
    Ch1norm_AVG_Gcampdata = AVGIntData./Ch1AVGIntData;

else
    disp('Subtracting background from Gcamp signal')
    slicedata = pixelsbyslicearray-Ch2bkgndpixelsbyslicearray;
    Ch2norm_AVG_Gcampdata = AVGIntData-Ch2bkgndAVGIntData;

end

disp('Aligning Neural Data')

%This final section of code aligns the neural data with the behavioral data
%(i.e. frame number) and interpolates data between the 2P measurements. 
cd(rootdir)
LaserData = load(diffArray); %Need to use this later for aligning 
% the timestamps.

LaserIsOn = NaN(length(LaserData.diffArray),1);

%figure out during what frames the laser is on vs off
indeces=1:length(LaserData.diffArray);
DetectOn=1;
for i = indeces
   if LaserData.diffArray(i,1) > laserThresh
       if DetectOn==1
           LaserIsOn(i,1) = 1;
           DetectOn=0;  
           if LaserData.diffArray(i+1,1) > laserThresh %we've already determined 
                    % that the laser is on. 
                DetectOn=1;
           end       
       else
           LaserIsOn(i,1)=0;
           DetectOn=1;
           if LaserData.diffArray(i+1,1) > laserThresh 
                LaserIsOn(i,1)=1;
                DetectOn=0; %if the laser takes a few frames to turn off, 
                % resort to the latest frame for the off time. 
           end
       end       
   else
       if DetectOn==1
           LaserIsOn(i,1)=0;
       else
           LaserIsOn(i,1)=1;
       end
   end
end

%get a list of the indeces where the laser is on. align the neural data to
%the onset of each laser flash

laserStartIndices = find(diff(LaserIsOn)==1); 
% laserStartIndices is finding all of the image #s where the laser goes
% from off to on. note: the laser will stay on for a few frames
% (2-5ish)and will vary

laserOffIndices = find(diff(LaserIsOn)==-1); 
%at the last laser off index, we want to stop the interpolation, because
%there is no more measured neural signal after that point


if FalseStart %if the T-series was started and then stopped again while the behav video stayed on
    for i = numlaserflashes+1:length(laserStartIndices)
        currentindx = laserStartIndices(i);

        for j = 3:(slicesPerStack+2)
            fullMovementAndBrainSignalSlicedata(currentindx, j) = slicedata(i,j-2); 
            %array with cols = frame number, behavioral output,
            %neural signal for each slice in order

        end

        if readChannel1
            fullMovementAndBrainSignal(currentindx,3) = Ch1norm_AVG_Gcampdata(i);
            %array with cols = frame number, behavioral output, avg neural
            %signal 
        else
            fullMovementAndBrainSignal(currentindx,3) = Ch2norm_AVG_Gcampdata(i);
        end
    end

else
     for i = 1:3363
        currentindx = laserStartIndices(i);
        for j = 3:(slicesPerStack+2)
            
            fullMovementAndBrainSignalSlicedata(currentindx, j) = slicedata(i,j-2);
        end

        if readChannel1
            fullMovementAndBrainSignal(currentindx,3) = Ch1norm_AVG_Gcampdata(i);
   
        else
            fullMovementAndBrainSignal(currentindx,3) = Ch2norm_AVG_Gcampdata(i);
        end
    end

 for i = 3364:length(laserStartIndices)-1
        currentindx = laserStartIndices(i+1);
        for j = 3:(slicesPerStack+2)
            
            fullMovementAndBrainSignalSlicedata(currentindx, j) = slicedata(i,j-2);
        end

        if readChannel1
            fullMovementAndBrainSignal(currentindx,3) = Ch1norm_AVG_Gcampdata(i);
   
        else
            fullMovementAndBrainSignal(currentindx,3) = Ch2norm_AVG_Gcampdata(i);
        end
    end

end


%the code above aligns as much real neural data as we have to the onset of
%the laser. But now, we want to interpolate the data between readings
%so find the values that are still a NaN
isNumIndices = find(~isnan(fullMovementAndBrainSignal(:,3)));

%Interpolate the brain signal at these values so there is data for every
%behavior frame (some are measured, some are interpolated)
nanx = isnan(fullMovementAndBrainSignal(:,3));
t    = 1:length(fullMovementAndBrainSignal);
%t    = 1:162333;

interpBrainSignal = interp1(t(~nanx), fullMovementAndBrainSignal(~nanx, 3), t);

%insert the interpolated data into the array
fullMovementAndBrainSignal(1:length(interpBrainSignal),3) = interpBrainSignal(:);


%do the same for the slice data
for i =1:7
    disp(i)
    nanx = isnan(fullMovementAndBrainSignalSlicedata(:,i+2));
    t    = 1:length(fullMovementAndBrainSignalSlicedata (:, i+2));
    sliceinterpBrainSignal = interp1(t(~nanx), fullMovementAndBrainSignalSlicedata(~nanx, i+2), t);
    fullMovementAndBrainSignalSlicedata(1:numel(sliceinterpBrainSignal),i+2) = sliceinterpBrainSignal(:);
end


if readChannel1
    exportname = [Tseriesroot, '_Ch1Norm_', 'mvmtandbrain.xlsx'];
    sliceexportname = [Tseriesroot, '_Ch1Norm_', 'mvmtandbrainSLICES.xlsx'];
else
    exportname = [Tseriesroot, '_Ch2Norm_', 'mvmtandbrain.xlsx'];
    sliceexportname = [Tseriesroot, '_Ch2Norm_', 'mvmtandbrainSLICES.xlsx'];
end


disp('Applying Filters and Exporting Data')

cd(savedir)

%apply Savitzky-Golay filter
order = 3; %can play with these parameters if necessary
framelen = 201;
RawfiltwithSgo = sgolayfilt(fullMovementAndBrainSignal(:,3),order,framelen);
SliceRawfiltwithSgo = sgolayfilt(fullMovementAndBrainSignalSlicedata(:,3:end),order,framelen);

time_sec = 0:1/30:length(fullMovementAndBrainSignal)/30;
time_min = 0:1/(framesPerMin):length(fullMovementAndBrainSignal)/(framesPerMin);

time_min(1) = [];
time_sec(1) = [];


AVGSignalwithFilters = [time_sec' time_min' fullMovementAndBrainSignal RawfiltwithSgo]; 

SLICESignalwithFilters = [time_sec' time_min' fullMovementAndBrainSignalSlicedata SliceRawfiltwithSgo];

%create and export figures of data
figure(1)
yyaxis left
plot(AVGSignalwithFilters(:,2), AVGSignalwithFilters(:,5))%plotting raw pixel int on left y axis
%change to :,5 if there are 
ylabel('Average Raw Pixel Intensity')
xlabel('Time (min)')

yyaxis right
plot(AVGSignalwithFilters(:,2), AVGSignalwithFilters(:,4)) %plotting movement data on right y axis
ylabel('Movement Score')
title('Avg Raw Neural and Movement Data')

if readChannel1
    savefig([Tseriesroot,  '_Ch1Norm_', 'Avg_Raw_Neural_and_Movement_Data'])
else
    savefig([Tseriesroot,  '_Ch2Norm_', 'Avg_Raw_Neural_and_Movement_Data'])
end


figure(2)
yyaxis left
plot(AVGSignalwithFilters(:,2), AVGSignalwithFilters(:,6)) %plotting filtered pixel int on left y axis
ylabel('Filtered Avg Pixel Intensity')
xlabel('Time (min)')

yyaxis right
plot(AVGSignalwithFilters(:,2), AVGSignalwithFilters(:,4)) %plotting movement data on right y axis
ylabel('Movement Score')
title('Avg Filtered Neural and Movement Data')
if readChannel1
    savefig([Tseriesroot,  '_Ch1Norm_', 'Avg_Filtered_Neural_and_Movement_Data'])
else
    savefig([Tseriesroot,  '_Ch2Norm_', 'Avg_Filtered_Neural_and_Movement_Data'])
end

figure(3)
plot(SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,5),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,6),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,7),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,8),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,9),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,10),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,11),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,12));
   
hold on
ylabel('Raw Slice Pixel Intensity')
xlabel('Time (min)')

yyaxis right
plot(AVGSignalwithFilters(:,2), AVGSignalwithFilters(:,4)) %plotting movement data on right y axis
ylabel('Movement Score')
legend('Slice 1', 'Slice 2', 'Slice 3', 'Slice 4','Slice 5','Slice 6','Slice 7','Slice 8')

title('Raw Slice Neural and Movement Data')
savefig([Tseriesroot 'Raw_Slice_Neural_and_Movement_Data'])
if readChannel1
    savefig([Tseriesroot,  '_Ch1Norm_', 'Raw_Slice_Neural_and_Movement_Data'])
else
    savefig([Tseriesroot,  '_Ch2Norm_', 'Raw_Slice_Neural_and_Movement_Data'])
end


figure(4)
plot(SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,12),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,13),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,14),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,15),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,16),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,17),...
    SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,18)...
    );
% 
% use if there are 8 slices per stack
% plot(SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,13),...
%     SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,14),...
%     SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,15),...
%     SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,16),...
%     SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,17),...
%     SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,18),...
%     SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,19);
    %SLICESignalwithFilters(:,2), SLICESignalwithFilters(:,20))
hold on
ylabel('Filtered Slice Pixel Intensity')
xlabel('Time (min)')

yyaxis right
plot(AVGSignalwithFilters(:,2), AVGSignalwithFilters(:,4)) %plotting movement data on right y axis
ylabel('Movement Score')
legend('Slice 1', 'Slice 2', 'Slice 3', 'Slice 4','Slice 5','Slice 6','Slice 7','Slice 8')

title('Filtered Slice Neural and Movement Data')
if readChannel1
    savefig([Tseriesroot,  '_Ch1Norm_', 'Filtered_Slice_Neural_and_Movement_Data'])
else
    savefig([Tseriesroot,  '_Ch2Norm_', 'Filtered_Slice_Neural_and_Movement_Data'])
end

%export csv files of data
mvmtandbrain = array2table(AVGSignalwithFilters, ...
    'VariableNames',{'time_sec','time_min', 'FrameNum','MovementData','RawPixelIntensity'...
    ,'sgolayfiltdata'});

VarNames = {'time_sec','time_min', 'FrameNum','MovementData','Slice1_RawPixelIntensity',...
    'Slice2_RawPixelIntensity', 'Slice3_RawPixelIntensity', 'Slice4_RawPixelIntensity', ...
    'Slice5_RawPixelIntensity', 'Slice6_RawPixelIntensity', 'Slice7_RawPixelIntensity', ...
    'Slice1_sgolayfiltdata',...
    'Slice2_sgolayfiltdata', 'Slice3_sgolayfiltdata', 'Slice4_sgolayfiltdata', ...
    'Slice5_sgolayfiltdata', 'Slice6_sgolayfiltdata', 'Slice7_sgolayfiltdata'};

 slicemvmtandbrain = array2table(SLICESignalwithFilters, 'VariableNames', VarNames);

writetable(mvmtandbrain,exportname);
writetable(slicemvmtandbrain,sliceexportname);
writematrix(slicedata,'RawSliceDataBeforeNorm');
writematrix(AVGIntData,'RawAVGDataBeforeNorm');


