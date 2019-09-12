function filename = uigetfile_with_preview(filterSpec, prompt, folder, callback, multiSelect)
% uigetfile_with_preview - GUI dialog window with a preview panel
%
% Syntax:
%    filename = uigetfile_with_preview(filterSpec, prompt, folder, callback, multiSelect)
%
% Description:
%    uigetfile_with_preview displays a file-selection dialog with an
%    integrated preview panel.
%
% Inputs:
%    filterSpec - similar to uigetfile's filterSpec. Examples:
%                    'defaultFile.mat'
%                    '*.mat'
%                    {'MAT files', '*.mat'}
%                    {'Data or m files (*.m,*.mat)', {'m';'mat'}}
%    prompt     - optional dialog prompt (default: 'Select file')
%    folder     - optional default folder (default: current folder)
%    callback   - optional handle to a callback function with the following
%                 interface:
%                    callbackFunction(hPanel, filename)
%                 The callback function should check the specified file and
%                 update the specified panel accordingly.
%                 (default: [])
%    multiSelect - optional flag specifying whether multiple files can be
%                  selected (default: false)
%
% Output:
%    filename - full-path filename of the selected file(s) (or empty array)
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Change log:
%    2016-11-02: First version posted on <a href="https://www.mathworks.com/matlabcentral/fileexchange/?term=authorid%3A27420">MathWorks File Exchange</a>
%
% See also:
%    uigetfile

% Programmed by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.0 $  $Date: 2016/11/02 23:43:54 $

  if nargin < 1,  filterSpec = {'All files','*.*'};  end
  if nargin < 2,  prompt = '';          end
  if nargin < 3,  folder = pwd;         end
  if nargin < 4,  callback = '';        end
  if nargin < 5,  multiSelect = false;  end

  % Prepare the figure dialog window
  if isempty(prompt),  prompt = 'Select file';  end
  hFig = dialog('Name',prompt, 'units','pixel', 'pos',[200,200,500,600]);

  % Display the file-selection component
  javaComponentName = 'javax.swing.JFileChooser'; %#ok<NASGU>
  javaComponentName = 'com.mathworks.hg.util.dFileChooser';  % this is a MathWorks extension of JFileChooser - no good on Win7
  %javaComponentName = 'com.mathworks.mwswing.MJFileChooserPerPlatform';
  [hjFileChooser, hContainer] = javacomponent(javaComponentName, [0,200,500,400], hFig); %#ok<ASGLU>
  hjFileChooser.setCurrentDirectory(java.io.File(folder));
  hjFileChooser.setMultiSelectionEnabled(multiSelect);
  drawnow;

  % Java default background colors are different from Matlab's - fix this
  %bgcolor = hjFileChooser.Background.getComponents([]);
  %set(hFig,'color',bgcolor(1:3));
  bgcolor = get(hFig,'color');
  hjFileChooser.Background = java.awt.Color(bgcolor(1),bgcolor(2),bgcolor(3));

  % Prepare the allowable file types filter (similar to the uigetfile function)
  if ~isempty(filterSpec)
      if ischar(filterSpec)
          origFilterSpec = filterSpec;
          [dummy,dummy,ext] = fileparts(filterSpec); %#ok<ASGLU>
          filterSpec = {[ext(2:end) ' files'], ['*' ext]};
          if ~any(origFilterSpec=='*')
              try hjFileChooser.setSelectedFile(java.io.File(origFilterSpec)); catch, end
          end
      end
      hjFileChooser.setAcceptAllFileFilterUsed(false);
      fileFilter = {};
      %fileFilter{1} = AddFileFilter(hjFileChooser, 'Data files (*.mat)', {'mat'});
      %fileFilter{2} = AddFileFilter(hjFileChooser, 'Data or m files (*.m,*.mat)', {'m';'mat'});
      %fileFilter{3} = AddFileFilter(hjFileChooser, 'All Files', '*.*');
      for filterIdx = 1 : size(filterSpec,1)
          fileFilter{end+1} = AddFileFilter(hjFileChooser, filterSpec{filterIdx,:}); %#ok<AGROW>
      end
      try
          hjFileChooser.setFileFilter(fileFilter{1});  % use the first filter by default
      catch
          % never mind - ignore...
      end
  end

  % Prepare the preview panel
  hPreviewPanel = uipanel('parent',hFig, 'title','Preview', 'tag','PreviewPanel', ...
      'units','pixel', 'pos',[10,10,480,180], 'Background',bgcolor(1:3));

  % Prepare the figure callbacks
  hjFileChooser.PropertyChangeCallback  = {@PreviewCallback,hPreviewPanel,callback};
  hjFileChooser.ActionPerformedCallback = {@ActionPerformedCallback,hFig};

  % Key-typed callback
  try
      hFn = handle(hjFileChooser.getComponent(2).getComponent(2).getComponent(2).getComponent(1),'CallbackProperties');
      % hFn should be a javahandle_withcallbacks.com.sun.java.swing.plaf.windows.WindowsFileChooserUI$7
      hFn.KeyTypedCallback = {@KeyTypedCallback,hjFileChooser};
  catch
      % maybe the file-chooser editbox changed location
  end

  % Wait for user input
  uiwait(hFig);
  
  % We get here if the figure is either deleted or Cancel/Open were pressed
  if ishghandle(hFig)
      % Open were clicked
      filename = getappdata(hFig,'selectedFile');
      close(hFig);
  else  % figure was deleted/closed
      filename = '';
  end
end  % uigetfile_with_preview

%% Add a file filter type
function fileFilter = AddFileFilter(hjFileChooser, description, extension)
  try
      if ~iscell(extension)
          extension = {extension};
      end
      if strcmp(extension{1},'*.*') %any(extension{1} == '*')
          jBasicFileChooserUI = javax.swing.plaf.basic.BasicFileChooserUI(hjFileChooser.java);
          fileFilter = javaObjectEDT('javax.swing.plaf.basic.BasicFileChooserUI$AcceptAllFileFilter',jBasicFileChooserUI);
      else
          extension = regexprep(extension,'^.*\*?\.','');
          fileFilter = com.mathworks.mwswing.FileExtensionFilter(description, extension, false, true);
      end
      javaMethodEDT('addChoosableFileFilter',hjFileChooser,fileFilter);
  catch
      % ignore...
      fileFilter = [];
  end
end  % AddFileFilter

%% Preview callback function
function PreviewCallback(hjFileChooser, eventData, hPreviewPanel, callback) %#ok<INUSL>
  try
      % Get the selected file
      filename = char(hjFileChooser.getSelectedFile);
      if isempty(filename) || ~exist(filename,'file')
          return;  % bail out
      end
  
      % Analyze the file based on its type
      ClearPanel(hPreviewPanel)
      if isempty(callback)  % default callback
          [fpath,fname,fext] = fileparts(filename); %#ok<ASGLU>
          switch fext
              case '.mat'
                   hax = axes('Parent',hPreviewPanel, 'units','norm', 'LooseInset',[0,0,0,0]);
                  % Load the first data element in the file
%                   data = load(filename);
%                   fnames = fieldnames(data);
%                   data = data.(fnames{1});

                  SSS=open(filename);
    %           uiopen('load')
                  TempPart=SSS.Psave.GUISavePart ;
                  
                  % Update the preview panel
                 DrawPartH( TempPart,'',1,gca ,'')
%                   plot(hax, data.x, data.y);
                  box(hax, 'off');

              otherwise
                  % Do something useful here...
          end
      else  % user-specified callback that accepts hPreviewPanel,filename
          callback(hPreviewPanel,filename);
      end
  catch
      % Never mind - bail out...
  end
end  % PreviewCallback

%% Clear the specified panel
function ClearPanel(hPanel)
  try
      hObjs = findall(hPanel);
      hObjs = setdiff(hObjs, findall(hPanel,'string',get(hPanel,'title')));  % keep the title
      %delete(hObjs(2:end));
      hObjs = setdiff(hObjs, hPanel);  % this is safer...
      delete(hObjs);
  catch
      % ignore...
  end
end  % ClearPanel

%% Figure actions (Cancel & Open)
function ActionPerformedCallback(hjFileChooser, eventData, hFig)
  switch char(eventData.getActionCommand)
      case 'CancelSelection'
          close(hFig);

      case 'ApproveSelection'
          files = cellfun(@char, cell(hjFileChooser.getSelectedFiles), 'uniform',0);
          %msgbox(['Selected file: ' files{:}], 'Selected!');
          if numel(files)==1
              %files = files{1};
          end
          if isempty(files)
              files = char(hjFileChooser.getSelectedFile);
          end
          setappdata(hFig,'selectedFile', files);
          uiresume(hFig);

      otherwise
          % should never happen
  end
end  % ActionPerformedCallback

%% Key-types callback in the file-name editbox
function KeyTypedCallback(hEditbox, eventData, hjFileChooser) %#ok<INUSL>
  text = char(get(hEditbox,'Text'));
  [wasFound,idx,unused,folder] = regexp(text,'(.*[:\\/])'); %#ok<ASGLU>
  if wasFound
      % This will silently fail if folder does not exist
      hjFileChooser.setCurrentDirectory(java.io.File(folder));
  end
end
