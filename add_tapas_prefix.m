function add_tapas_prefix()

% directory where the renaming of functions is to be done
% pathToolbox = 'C:\Users\kasperla\Documents\Promotion\experiments\monitoring\smoothing\trunk\tSNR_fMRI_SPM\CheckPhysRETROICOR\PhysIOToolbox';
pathToolbox = 'C:\Users\kasperla\Documents\Promotion\experiments\monitoring\smoothing\trunk\tSNR_fMRI_SPM\CheckPhysRETROICOR\test3';


% folder with all the function names to be replaced in all
% functions/scripts
dirCode = fullfile(pathToolbox, 'code');


% folders where scripts/functions with renamed function names occur
pathScripts = pathToolbox;
dirScripts  = {
    dirCode %comment this line, if functions itself are already renamed
    fullfile(pathScripts, 'examples', 'GE', 'Pulseoxy_15m4')
    fullfile(pathScripts, 'examples', 'Philips', 'ECG3T')
    fullfile(pathScripts, 'examples', 'Philips', 'ECG7T')
    fullfile(pathScripts, 'examples', 'Philips', 'PPU')
    };

% prefix to be added
prefix = 'tapas_';

% collect all filenames from code and script directories
filelistCode= dir(fullfile(dirCode,'*.m'));

filelistScripts = [];
for d = 1:length(dirScripts)
    tmp = dir(fullfile(dirScripts{d},'*.m'));
    for i = 1:length(tmp)
        tmp(i).name = fullfile(dirScripts{d}, tmp(i).name);
    end
    filelistScripts = [filelistScripts;tmp];    
end


% go through each file and replace filename
for i=1:length(filelistScripts)
    fnin = filelistScripts(i).name;
    fnout = [fnin '.tmp'];
    fin = fopen(fnin);
    fout = fopen(fnout,'w');
    while ~feof(fin) % cycle through all lines of the script
        s = fgetl(fin);
        for j=1:length(filelistCode) %cycle through all code-function names to be replaced
            [fpathstr,fname,fext] = fileparts(fullfile(dirCode,filelistCode(j).name));
            s = strrep(s, fname, strcat(prefix,fname));
            s = strrep(s, [prefix prefix], prefix); %undo double replacements
        end
        fprintf(fout,'%s\n',s);
    end
    fclose(fin);
    fclose(fout);
    % rename file
    movefile(fnout, fnin);
end

% 
% 
% % rename all code files
isVersionControlled = false;
for j=1:length(filelistCode)
    movefile(fullfile(dirCode,filelistCode(j).name),fullfile(dirCode,[prefix,filelistCode(j).name]))
    if isVersionControlled
        unix(['svn mv ' fullfile(dirCode,filelistCode(j).name), fullfile(dirCode,[prefix,filelistCode(j).name])]);
    end
end
