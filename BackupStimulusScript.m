function BackupStimulusScript( scriptPath, scriptName, save_dir )
% Create a backup copy of the running script in the folder:
% save_dir\usedCode\
% Entry this line at the beginning of the caller script:
% [scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if ~exist('save_dir', 'var')
    save_dir = '.';
end

folder = [save_dir filesep 'usedCode' filesep];
datesuffix = datestr(now,'yyyymmdd_HHMMSS');
backupname = [folder scriptName '-backup-' datesuffix '.m'];
if ~exist(folder, 'dir')
    mkdir( folder )
end
copyfile( [scriptPath filesep scriptName '.m'] , backupname );
disp(['Backup copy of the code : ' scriptName]);