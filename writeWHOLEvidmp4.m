
close all; clear all;

%%%%user notes: if you want to use this code to a shorter clip, change the variable frameStart in variable list %%%%below to the first frame you want in your video, then scroll down into the code to line 24, where the %%%%allVidsEnd variable is defined and change it from allVidsEnd = numFrames-1; to allVidsEnd = last frame you %%%%want in your vid

%%%%%%%%%%%%%START parameters to modify:
rootdir = 'F:\'; 
flyTiffFolder = '2022.02.01 behavior'; %Location of the frames. Assuming this is a subfolder of rootdir (above).
behaviorImgRoot = 'fc2_save_2022-02-01-112946-'; %everything but the zero padded 4 to 6 digit frame rank order.
recordedFPS = 30;
fps2write = 60; %at 60 this is set to write a video at 2x speed, change back to 30 for real time
frameStart = 1;
%%%%%%%%%%%%%END parameters to modify:

xSpeedVal = round(fps2write/recordedFPS*100)/100;
framesPerVid = 5*60*recordedFPS;
cd(rootdir);
% First, want to take a look at how many frames we have:
cd(flyTiffFolder);
tiffList = dir([behaviorImgRoot '*.pgm']);
%display(numel(tiffList));
numFrames = numel(tiffList);

allVidsEnd = numFrames-1; 
%allVidsEnd = 136800;

numZeroPad = 4; %Number of zeros in first (zeroth) frame - note that Flycap will automatically increase the number of digits once this is exceeded.
currentUpperBound = 10^(numZeroPad+1);

vidObj = VideoWriter([behaviorImgRoot '_' num2str(xSpeedVal) 'Xspeed_frame' '_frame' num2str(frameStart) '_toframe' num2str(allVidsEnd) '.mp4'],'MPEG-4');
vidObj.FrameRate = fps2write;
open(vidObj);
    
cd(rootdir)
cd(flyTiffFolder);
 
for(fi = frameStart: allVidsEnd) %numel(tiffList)),
%for (fi= frameStart: frameEnd)
    if(fi>currentUpperBound)
        numZeroPad = numZeroPad + 1;
        currentUpperBound = 10^(numZeroPad+1);
    end

    imgName = [behaviorImgRoot sprintf(['%0' num2str(numZeroPad) '.0f'],fi) '.pgm'];
    jpegName = [behaviorImgRoot sprintf(['%0' num2str(numZeroPad) '.0f'],fi) '.jpg'];

    if(exist(imgName,'file'))
        thisFrame = imread(imgName);
        graythisFrame = im2gray(thisFrame);
        writeVideo(vidObj,graythisFrame);
%         cd(rootdir)
%         cd(jpegsavedir)
%         imwrite(graythisFrame, jpegName, 'jpg');
%         cd(rootdir)
%         cd(flyTiffFolder)
    end
%         if(mod(fi,100)==0),
        %display(fi);
%             toc;
%         end;
end
   
close(vidObj);
display(['frameStart = ' num2str(frameStart) ', last Frame Written: ' imgName]);

