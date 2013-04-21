function [ons_secs, outliersHigh, outliersLow] = physio_correct_cardiac_pulses_manually(ons_secs,percentile,upperThresh,lowerThresh)

% this function takes the onsets from ECG measure and controls for
% outliers (more or less than a threshold given by a percentile increased
% or decreased by upperTresh or lowerThresh percent respectively.

[outliersHigh,outliersLow,fh] = selectOutliers(ons_secs.t, ons_secs.cpulse, percentile, upperThresh, lowerThresh);
if any(outliersHigh)
    additionalPulse=[];
    fh2=figure;
    for outk=1:length(outliersHigh)
        s=0;
        while ~(s==1)
            indStart = outliersHigh(outk)-1; indEnd = outliersHigh(outk)+2;
            ind=find(ons_secs.t==ons_secs.cpulse(indStart))-100:find(ons_secs.t==ons_secs.cpulse(indEnd))+100;
            figure(fh2); clf;
            plot(ons_secs.t(ind),ons_secs.c(ind),'r')
            hold on;
            plot(ons_secs.cpulse(indStart:indEnd),ones(4,1)*max(ons_secs.c(ind)),'ok')
            inpNum=input('How many triggers do you want to set? Enter a number between 0 and 10 : ');
            I1=[];
            for ii=1:inpNum
                figure(fh2);
                [I1(ii), J1] = ginput(1);
                plot(I1(ii),J1, 'b*', 'MarkerSize',10);
                
            end
            s=input('If you agree with the selected triggers, press 1 (then enter) : ');
            if isempty(s)
                s=0;
            end
        end
        additionalPulse=[additionalPulse;I1'];
    end
    ons_secs.cpulse = sort([ons_secs.cpulse;additionalPulse]);
    close(fh2);
end
close(fh);


[outliersHigh,outliersLow,fh] = selectOutliers(ons_secs.t, ons_secs.cpulse, percentile, upperThresh, lowerThresh);
finalIndex=1:length(ons_secs.cpulse);
if any(outliersLow)
    fh3=figure;
    for outk=1:length(outliersLow)
        s=0;
        while ~(s==1)
            indStart = outliersLow(outk)-2; indEnd = outliersLow(outk)+2;
            ind=find(ons_secs.t==ons_secs.cpulse(indStart))-100:find(ons_secs.t==ons_secs.cpulse(indEnd))+100;
            figure(fh3); clf;
            plot(ons_secs.t(ind),ons_secs.c(ind),'r')
            hold on;
            plot(ons_secs.cpulse(indStart:indEnd),ones(5,1)*max(ons_secs.c(ind)),'ok','MarkerFaceColor','r');
            alreadyDeleted=intersect(indStart:indEnd,setdiff(1:length(ons_secs.cpulse),finalIndex));
            plot(ons_secs.cpulse(alreadyDeleted),ones(size(alreadyDeleted))*max(ons_secs.c(ind)),'or','filled');
            for kk=indStart:indEnd
                text(ons_secs.cpulse(kk),max(ons_secs.c(ind))*1.05,int2str(kk));
            end
            
            delInd= [];
            
            delInd=input('Enter the number of a pulse you want to delete (0 if none): ');
            plot(ons_secs.cpulse(delInd),max(ons_secs.c(ind)), 'rx', 'MarkerSize',20);
            
            s=input('If you agree with the deleted triggers, press 1 (then enter) : ');
            if isempty(s)
                s=0;
            end
            finalIndex=setdiff(finalIndex,delInd');
        end
        
        close(fh3);
    end
    ons_secs.cpulse = sort(ons_secs.cpulse(finalIndex));
end
close(fh);


end


function [outliersHigh,outliersLow,fh] = selectOutliers(t,tCardiac,percentile,deviationPercentUp,deviationPercentDown)

fh = get_default_fig_params();
set(fh, 'Name','Diagnostics raw phys time series');
dt = diff(tCardiac);

plot(tCardiac(2:end), dt);
xlabel('t (seconds)');
ylabel('lag between heartbeats (seconds)');
title('temporal lag between heartbeats');

nBins = length(dt)/10;
[dtSort,dtInd]=sort(dt);
percentile=percentile/100;
upperThresh=(1+deviationPercentUp/100)*dtSort(ceil(percentile*length(dtSort)));
lowerThresh=(1-deviationPercentDown/100)*dtSort(ceil((1-percentile)*length(dtSort)));
outliersHigh=dtInd(find(dtSort>upperThresh));
outliersLow=dtInd(find(dtSort<lowerThresh));

if nnz(dt > upperThresh)
    text(t( find(dt==max(dt))+1 ),max(dt),...
        {'Warning: There seem to be skipped heartbeats R-waves in the scanner-log', ...
        sprintf('first at timepoint %01.1f s',tCardiac(min(outliersHigh+1))), ...
        'rerun read\_physlog with ECG\_min set to 1'}, ...
        'Color', [1 0 0])
end

end