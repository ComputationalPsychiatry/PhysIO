function physio_update_tapas_fork_git(pathTargetLocalRepository, urlSource)
% Updates PhysIO-subfolder of private TAPAS fork with changes from public
% folder of this repository
%
%   output = physio_update_tapas_fork_git(input)
%
% Uses git subtree and a split, as specified in:
%   http://stackoverflow.com/questions/23937436/add-subdirectory-of-remote-repo-with-git-subtree
% updates of subtree as in
%   https://www.kernel.org/pub/software/scm/git/docs/howto/using-merge-subtree.html
%
% IN
%   pathTargetLocalRepository   local path of repository, whose remote is a
%                               fork of tapas
%   urlSource                   remote location of repository, whose
%                               subfolder should be merged into TAPAS
% OUT
%
% EXAMPLE
%   physio_update_tapas_fork_git
%
%   See also
%
% Author: Lars Kasper
% Created: 2016-10-11
% Copyright (C) 2016 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id: teditRETRO.m 775 2015-07-17 10:52:58Z kasperla $


%% default input parameters
if nargin < 2
    urlSource = 'git@tnurepository.ethz.ch:lkasper/physio.git';
end

if nargin < 1
    pathTargetLocalRepository = '/Users/kasperla/Documents/code/matlab/tapas';
end

%% deployment options
branchNameSource = 'publicPhysio';
relativePathSourcePhysIO = 'public'; % name of subdir of urlSource-repository which should be transferred
relativePathTargetPhysIO = 'PhysIO'; % name of subdir of Tapas where Physio resides

% different execution depending on whether this is the first time we try
% it...
isFirstDeployment = true; % true only, if PhysIO Sub-directory does not exist
doRemoveCommitHistory = false;
useSubTree = true;

%% Execution

if doRemoveCommitHistory
    flagCommitHistory = '--squash';
else
    flagCommitHistory = '';
end

% Go to local repository of Tapas
pathTmp = pwd;
cd(pathTargetLocalRepository);

%% Do this the first time:

if useSubTree % does somehow not properly include history of all commits...
    
    if isFirstDeployment
        unix(sprintf('git remote add -f -t master --no-tags %s %s', branchNameSource, urlSource));
        unix(sprintf('git checkout %s/master', branchNameSource));
        unix(sprintf('git subtree split -P %s -b temporary-split-branch', relativePathSourcePhysIO));
        unix(sprintf('git checkout master'));
        
        if exist(relativePathTargetPhysIO, 'dir')
            unix(sprintf('git rm -rf %s', relativePathTargetPhysIO));
            unix(sprintf('git commit -m ''Deleted PhysIO folder for reinsertion with commit history'''));
        end
        
        unix(sprintf(...
            'git subtree add %s -P %s temporary-split-branch', ...
            flagCommitHistory, relativePathTargetPhysIO));
        
        unix(sprintf('git branch -D temporary-split-branch'))
        
    else
        %% In future, you can merge in additional changes as follows:
        unix(sprintf('git checkout %s/master', branchNameSource));
        unix(sprintf('git subtree split -P %s -b temporary-split-branch', relativePathSourcePhysIO));
        unix(sprintf('git checkout master'));
        
        unix(sprintf('git subtree merge %s -P %s temporary-split-branch', ...
            flagCommitHistory, relativePathTargetPhysIO));
        % Now fix any conflicts if you'd modified third_party/git-completion.
        % maybe needed for conflict resolution
        % unix(sprintf('git checkout --theirs *'));
        % unix(sprintf('git add -u'));
        
        unix(sprintf('git branch -D temporary-split-branch'));
    end
    
else % don't use subtree
    %% TODO: third strategy: merge into a new branch physio
    % git checkout -b physio master
    
    if isFirstDeployment
        unix(sprintf('git remote add -f -t master --no-tags %s %s', branchNameSource, urlSource));
        unix(sprintf('git merge -s ours --no-commit %s/master', branchNameSource));
        
        if exist(relativePathTargetPhysIO, 'dir')
            unix(sprintf('git rm -rf %s', relativePathTargetPhysIO));
        end
        
        unix(sprintf('git read-tree --prefix=%s -u %s/master:%s', ...
            relativePathTargetPhysIO, branchNameSource, relativePathSourcePhysIO));
        unix(sprintf('git commit -m ''Merged PhysIO changes into TAPAS'''));
    else
        % In future, you can *overwrite* with the latest changes as follows:
                unix(sprintf('git merge -s ours --no-commit %s/master', branchNameSource));
                unix(sprintf('git rm -rf %s', relativePathTargetPhysIO));
                unix(sprintf('git read-tree --prefix=%s -u %s/master:%s', ...
                    relativePathTargetPhysIO, branchNameSource, relativePathSourcePhysIO));
                unix(sprintf('git commit -m ''Merged PhysIO changes into TAPAS'''));
    %or, untested...does not seem to work with subdirectories:
        %unix(sprintf('git pull -s subtree %s master', branchNameSource));
    end
    
    
end
cd(pathTmp);