% Script chkpt.m
% Purpose:
% This script m-file is part of the checkpoint restart package for
% users to periodically save data to a file during a batch job
% In the event that the batch job got terminated, a user can rerun the job
% from the point of the last save instead of starting from the beginning.
% Input:
% matfile -- all data need saved go to this mat-file (e.g., myfile.mat)
% s is a struct array that stores the various data that the application
% code need save for restarting

for k=1:nNames
   s.(chkNames{k}) = eval(chkNames{k});  % update all variables first
end
save(matfile, '-struct', 's');