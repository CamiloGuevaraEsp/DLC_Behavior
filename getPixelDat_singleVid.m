
%getPixelDat_singleVid.m

%%%%NOTE: in its current form, this code cannot be used to measure movement 
%%%%using a frame subtraction method because the background image was set to
%%%%be the very first frame and never updates. This allows it to measure the 
%%%%contrast that occurs when the laser turns on. however to measure movement 
%%%%of the fly, an average frame wouldneed to be continuously re-measured 

%%%%%%%%%%%%%START parameters to modify:
beh_rootdir = 'F:\020222_processed121522'; %the location of the video
vidRoot = 'fc2_save_2022-02-02-110737-_2Xspeed_frame'; 
vidType = '.mp4';
frameStart = 1;
recordedFPS = 30;
%%%%%%%%%%%%%END parameters to modify:

cd(beh_rootdir);

vid2read = [vidRoot vidType];
bgImgName = strrep(vid2read, vidType,'.png');

disp('reading in video')
vid = VideoReader(vid2read);
frameEnd = vid.NumFrames;
thisFrame = read(vid,1); 

if ndims(thisFrame) ==3 %if the image is being read as a RGB, convert to grayscale
    disp('converting to grayscale')
    thisFrame = rgb2gray(thisFrame);
end

%just using the first frame since this code 
% is only being used to detect if the laser is on and not movement
BGFrame = double(thisFrame);
display(['Writing ' bgImgName]);
imwrite(uint8(BGFrame),bgImgName);
bgImg = imread(bgImgName);

vid2write = strrep(bgImgName,'.png','.mp4'); %Need this for also creating the diffArray output.
diffArray = NaN(frameEnd,1);

for(ti = frameStart:(frameEnd-1))

    %framePair = read(vid,ti:ti+1);
    frame1 = read(vid, ti);
    gframe1 = rgb2gray(frame1);
    frame2 = read(vid, ti+1);
    gframe2 = rgb2gray(frame2);

    %if it is still not working, maybe read each image individually convert
    %to grayscale and then define the frame pair..
    framePair = [gframe1, gframe2];

    bgSubImgA = imsubtract(gframe1, bgImg); 
    %take first of the pics in the pair and subtract the bgImg
    bgSubImgB = imsubtract(gframe2, bgImg); 


    if(ti>frameStart)
        %Previously the new ref image as declared bgSubImgB.
        bgImgDiffA =  double(bgSubImgA)-double(bgSubImgB);
        IA = uint8(abs(bgImgDiffA));
        bwIA = double(im2bw(IA,0.1));
        %disp(sum(bwIA(:)))
        diffArray(ti) = sum(bwIA(:));

    else
        %take 2nd of the pics in the pair and subtract the bgImg
        bgImgDiffB = double(bgSubImgB)-double(bgSubImgA);
        %this is the diff that we want to append to the diffArray data, but
        %still in shape of image
        IB = uint8(abs(bgImgDiffB)); %abs val of diffs
        bwIB = double(im2bw(IB,0.1)); %converted to a binary.. 
        % not sure why that is necessary...but fine. 
        diffArray(ti) = sum(bwIB(:));

    end
end

save(strrep(vid2write,'.mp4','.mat'),'diffArray');
