function [] = figLog(filelocation,nameoffile) 
% creates a simple log file with date stamp and direcory of m-file

% to be used in conjunction with saving figures and data, so I can trace
% where they came from.


outdir = '/Users/kanestone/Dropbox/Work_Share/MITwork/figLogFiles/';
outfile = [outdir,nameoffile,'_log.txt'];

fid = fopen(outfile,'a');
if fid == -1
  error('Cannot open log file.');
end
fprintf(fid, '%s: %s\n', datestr(now, 0), filelocation);
fclose(fid);

end
