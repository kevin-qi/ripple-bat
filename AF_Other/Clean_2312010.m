
% % Open the file to count total number of lines
% fileID = fopen('231210_cdp_1.txt', 'r');
% totalLines = 0;
% 
% % Count the number of lines
% while ~feof(fileID)
%     fgetl(fileID);
%     totalLines = totalLines + 1;
% end
% fclose(fileID);

% Parameters
filename = '231210_cdp_1.txt';
numLinesToRead = 1597757;%16708706; % Change this to the number of lines you want to read

% Open the file
fileID = fopen(filename, 'r');

% Preallocate cell array to hold the data
data = cell(numLinesToRead, 1);

% Read the first numLinesToRead lines
for i = 1:numLinesToRead
    line = fgetl(fileID);
    if ~ischar(line)
        data = data(1:i-1); % Trim cell array if end of file is reached early
        break; % End of file reached
    end
    data{i} = line;
    
     if ~mod(i,10000)
        disp((i/numLinesToRead)*100)
    end
end

% Close the file
fclose(fileID);

numericData_size = cellfun(@(x) size(str2num(x),2), data, 'UniformOutput', false); 
data = data([numericData_size{:}]==8,1);


% Convert the cell array to a matrix
% Assuming the data is numeric and space-separated
numericData = cellfun(@(x) str2num(x), data, 'UniformOutput', false); %#ok<ST2NM>
numericData = vertcat(numericData{:});

% Assuming numericData is your array
filename = '231210_cdp_1_restored.txt';

% Save the numeric data to a text file with comma delimiter
writematrix(numericData, filename, 'Delimiter', ',');
disp(['Data saved to ', filename]);

