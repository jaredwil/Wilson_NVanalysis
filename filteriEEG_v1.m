function [filtData] = filteriEEG_v1(data, )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    %%%%%%%%%%%%%%%%%%%%%%%Removed: use external filtering function%%%%%%%%%%%%
    %         if filtFlag
    %             for c = 1:numel(channels)
    %                 [b a] = butter(3,[1/(fs/2)],'high');
    %                 d1 = filtfilt(b,a,blockData(:,c));
    %                 try
    %                 [b a] = butter(3,[70/(fs/2)],'low');
    %                 d1 = filtfilt(b,a,d1);
    %                 catch 
    %                 end
    %                 [b a] = butter(3,[58/(fs/2) 62/(fs/2)],'stop');
    %                 d1 = filtfilt(b,a,d1);
    % 
    %                 if filtCheck 
    %                     figure;
    %                     subplot(2,1,1); plot(blockData(:,c)); title('original');
    %                     subplot(2,1,2); plot(d1); title('d1');
    %                     linkaxes;
    %                     chk = input('Check plot for correct filter results, enter y for yes to stop check: ','s');
    %                     if strcmp(chk,'y'); filtCheck = 0; end
    %                 end
    %                 blockData(:,c) = d1;
    %             end
    %         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

