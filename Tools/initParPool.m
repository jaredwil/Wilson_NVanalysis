function [p, poolsize ] = initParPool( N )
%initParPool( N ) starts a parallel pool with a max of N workers (depending
%on if N workers is available) This function is used to limit workers when
%working on a shared resource.
%   Function returns both p - info about parallel pool, and the poolsize
%   for reference in par inplementation
%%%%%%%%%%%%%%%%%%%%%% start par pool%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Establish parpool
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)        %if not p will contain all info about current pool;
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
%initialize paralell pool if available
if(poolsize == 0)

    myCluster = parcluster('local');
    numWork = myCluster.NumWorkers;

    if(numWork <= 2)  %processing is done on a laptop so don't do it in parallel
        parpool('local',1)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    elseif(numWork > N) %limit the number of workers to 6
        parpool('local',N)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    else  %set up a parallel pool with max number of workers available between 2 and 6
        parpool(myCluster)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

