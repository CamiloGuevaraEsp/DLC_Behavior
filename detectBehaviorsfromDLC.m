close all; clear all;


%%%%%%%%%%%%%START parameters to modify:
movement_csv = 'fc2_save_2022-01-07-162840-DLCmvmt.csv';
file_path = '/Users/theresapatten/Downloads/';
%%%%%%%%%%%%%END parameters to modify:


pat = digitsPattern(2,4);
dateOfExp = extract(movement_csv,pat);
dateOfExpString = strjoin(dateOfExp((1:3),1));
cd(file_path);
M = readtable(movement_csv);
mvmtscore = M(:,11);
mvmtscorearray = table2array(mvmtscore);
avgmvmt = mean(mvmtscorearray,1);
baseline = 2*avgmvmt;
BLadjustedmvmtscore = mvmtscorearray-baseline;

for ii = 1:length(BLadjustedmvmtscore)
    if BLadjustedmvmtscore(ii)<0
       BLadjustedmvmtscore(ii) = 0;
    end
end

RestArray = nan(length(BLadjustedmvmtscore), 4);
framenums = 1:length(BLadjustedmvmtscore);
secs = framenums*(1/30);
mins = secs/60;
RestArray(:,1) = framenums;
RestArray(:,2) = secs;
RestArray(:,3) = mins;

plot(mins,BLadjustedmvmtscore, 'color', 'k')
xlabel('Time (mins)')
ylabel('Movement Score (AU)')
xlim([0 90])
hold on

resting = 0;
framesatrest = 0;
EventArray = [];

for ii = 1:length(BLadjustedmvmtscore)
    if resting == 0 %fly is awake
        if BLadjustedmvmtscore(ii)<20 %restperiodstarts
            %disp('fly is at rest')
            resting = 1;
            startRestFrame = RestArray(ii,1);
            startRestTime = RestArray(ii,3);
            RestArray(ii,4) = 0;
            framesatrest = framesatrest +1;
        else %fly is still awake
           RestArray(ii,4) = 1;
           %disp('fly is still awake')
        end
    end

    if resting == 1
        if BLadjustedmvmtscore(ii)<20 %fly is still asleep
            %disp('fly is still resting')
            RestArray(ii,4) = 0;
            framesatrest = framesatrest +1;
        else %fly has moved
            %disp('fly has moved')
            RestArray(ii,4) = 1;
            %how long has the fly been at rest...only record rest periods >2
            %mins long
            restduration = RestArray(ii,3) - startRestTime;
            activecount = 0;
            if ii+60<length(BLadjustedmvmtscore)
               for j = ii+30:ii+60
                   if BLadjustedmvmtscore(j)>20 %is there any movement 
                       % in the 2nd second after the initial movement
                       activecount = activecount+1;
                   end
               end  
           else 
               for j = ii+30:length(BLadjustedmvmtscore)
                   if BLadjustedmvmtscore(j)>20 %is there any movement 
                       % in the 2nd second after the initial movement
                       activecount = activecount+1; 
                   end
               end  
           end

           if restduration > 3 && activecount > 0 %it's a rest period worth noting
               %disp('rest duration is long enough and movement persisted for > 1 min...recording rest to wake transition')
               endRestFrame = RestArray(ii,1);
               endRestTime = RestArray(ii,3);
               %add this rest event to the array of events
               EventArray = [EventArray; startRestFrame, endRestFrame, startRestTime, endRestTime, restduration];
               resting = 0; %the fly is no longer resting
               %disp('fly is now awake')
               rectangle('Position',[startRestTime, max(BLadjustedmvmtscore)+0.05, restduration, 15], 'FaceColor','k')

               framesatrest = 0;

           elseif restduration < 3 && activecount > 0 
               %disp('fly movement persisted past 1 second - fly is awake')
               resting = 0;
               framesatrest = 0;

           elseif restduration < 3 && activecount == 0 
               %disp('fly twitched, but movement did not last more than one second - fly is still at rest')
               resting = 1;
           %then the fly twitched, 
               % but it is still considered to be resting (resting still
               % =1)
           end
        end
    end
end

savefigname = [dateOfExpString, '_Movmement.jpg'];
savefig(savefigname)
savename = [dateOfExpString, '_ListofRestPeriods.csv'];
varnames = {'Rest Start Frame', 'Rest End Frame', 'Rest Start Time', 'Rest End Time', 'Rest Duration'};
EventList = array2table(EventArray, "VariableNames", varnames);
writetable(EventList, savename);
