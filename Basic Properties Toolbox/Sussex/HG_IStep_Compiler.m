%% HG_IStep_Compiler

% Description: 

% This is a short script that reads through all of the output files of the
% IStep.m function in a given directory and merges the data into a more
% easily accessible format. 

% Throughout, sections are in bold (%%) with notes and suggestions in each
% section (%). Bits to be deleted and will be listed as blanks ([]). 


% You will need to search the functions and work out their syntax, MATLAB
% provides a lot of help for this. You can type help and then the function
% name into the command window, or typing the function name into the search
% bar in the top right, or even into google. 



%% Limit the scope of the search, then list all the files of interest

% Why? 
%   It is important to limit the scope of this script, so that it doesn't
%   search your entire computer (potentially taking an age and finding
%   files you're not interested in merging into your spreadsheet).

% How? 
%   The best way to do this would be create a variable of the full 
%   filepath of the folder you've done your analysis in. Then we can use
%   the 'dir' function to search that directory for a file of interest. The
%   'pwd' funcition may be useful to get the full filepath in the way
%   matlab likes it written.

% Create the pathstr variable that will have your filepath in it
pathstr = []; % eg. pathstr = 'C:\OGSteele\Somefolder\Someotherfolder'

% Change directory to the search directory listed above with the 'cd'
% function
cd([]);

% List all of the files of interest in your directory using the 'dir'
% function. Do this by creating the list variable that will be a list of
% all of the contents of the directory that contain files ending in a
% certain pattern (ie, the file type such as '.xlsx', '.pptx' etc etc), so 
% be sure to define the correct pattern. 
list = dir([]); 

% The file type you've chosen can actually created for lots of other
% reasons, so it's important to apply a filter to this. We can now go down
% our list created above and search for whether the file that the 'dir'
% function found is actually one of ours, ie, does it start with a
% convention? Is there a pattern to the names of the files that IStep
% creates? This patten will need to be consistent across all of the
% files created, but also selective to only the files created by IStep and
% not created by another function, or just in your directory ... Once
% you've figured this out, create the pat variable with this in it. 
pat = []; % eg. 'filenamepattern' 

% Using the above pattern, we now actually need to go down the list one by
% one and get MATLAB to remove the files that are not part of this naming
% convention. To do this we'll create a for loop that is the same size as
% the 'list variable, that asks on each loog whether or not the file is one
% of ours, and then delete the ones we're not interested in from our list.
% The 'size' function may help to tell you the size of the list variable.
% Then inside the loop we'll use the 'startsWith' function to check whether
% a part of the list variable starts with the pattern you created in the
% 'pat' variable.
% Note, when using the 'size' function, it will often spit out two numbers -
% and you're only interested in one ... is there a way to get the 'size'
% function to only spit out one dimension? 
for i = 1:[]
    list([]).isfile = startsWith([],[]); 
end

% We will then create a variable that is the list of whether or not the
% file is one we want to include, and then delete the ones we're not
% interested in.
% Note, i've done this bit for you, as it's a little counter intuitive.
fileFlags = [list.isfile]; % list of whether or not the file is one we want to include
list = list(fileFlags); % save only the ones we're interested in

% Congratulations, you SHOULD now have a list that is only the files of
% interest! Check this visually to make sure there are none missing or
% there are some naughty files that have snuck through your filter. You may
% need to change your filter if that's the case ... 

%% Create a master structure by loading the files in the list. 

% Why? 
%   This section will then go down your list one by one and load the saved
%   output, store it, load the next and then squash them together into a
%   new variable and carry on.

% How? 
%   This bit will use for loops like we've used above and the 'load' 
%   function - nothing more complex than that. It's just a little odd
%   sometimes thinking about for loops that stitch things onto the end each
%   time, as usually you want a for loop to stay the same size (as we did
%   in the above example, whereas in this one we want it to grow.

% Start by loading the first file of interest and then renaming the
% variable to allow the for loop to work properly. We will need to use the
% 'fullfile' function here to stitch parts of the list variable together
% for the 'load' function to work reliably.
load(fullfile([],[])); % load the first file
compiled = []; % rename the first file to facilitate the for loop (note,
                % what is the file called in your workspace once loaded?)

% Create another for loop now that will run down the list variable and load
% the files in order and stitch them together. 
for i = []:[]
    load(fullfile(list([]).[],list([]).[]));
    compiled = [compiled [] ]; % stitch next file onto the existing file
end

% Congratulations, you should now have a variable called 'compiled' that is
% one giant structure of all the data in one section. 

%% Correct the input mismatch
% This is a bit of an odd one, but i'm not entirely sure why some of the
% waves are wrong - this corrects for them just by pushing the array back
% one. I've done this for you, as it would be a bit unfair to get you to do
% this! 

% correct for input mismatch
x = zeros(31,1); % empty array of cells
for i = 1:size(compiled,2)
    if size(compiled(i).numSpikes,1) == 30
        spikesHolder = x; % assign holder variable the empty array
        spikesHolder(2:end) = compiled(i).numSpikes;
        compiled(i).numSpikes = spikesHolder;
    else
end

%% Organise action potential spike numbers into a similar format to excel

% Why? 
%   Looking at your newly created master structure, you'll notice that the
%   numSpikes array is not organised very conveniently for Excel. Let's fix
%   that! 

% How? 
%   Yet another for loop, and get the data transposed ready for excel. This
%   section will teach you how to accces data within structures. Accessing
%   a variable in a structure is done with the '.' seperator between the
%   structure name and the field within the structure. If you then want to
%   only look at a particular row, add the number of the row before the
%   seperator in brackets ie, '(x).' 

% Let's create a new blank variable (numSpikes) that will be filled with 
% the number of spikes for each recordings that is is roughly the right 
% size (ie, the number of waves, and then the number of recordings).
% Remember MATLAB does collumns by row by default ... think about the way
% you want the data organised to fit with your excel template sheet
numWaves = size([]); % number of waves in our current step protocol
numRecordings = size([]); % number of recordings analysed
numSpikes = zeros([],[]); % our empty variable that will eventually include
                            % the number of spikes in each wave of each
                            % recording

% Now, let's create that for loop that runs will reorganise our data.
% Essentially, we want this to go down our master structure and copy the
% data from our field that has the number of action potentials in it and
% paste it into our new variable (numSpikes, created above)
% Note, i've left this one completely blank for you to put together all of
% your knowledge gained from the above excercises to do this last step. 

