%% Option Data

clear
clc
close all

R = [];

%% Create connection to WRDS:

driver = eval('org.postgresql.Driver');
dbURL = 'jdbc:postgresql://wrds-pgdata.wharton.upenn.edu:9737/wrds?ssl=require&sslfactory=org.postgresql.ssl.NonValidatingFactory';

try
    % Attempt to load credentials from config.m
    run('config.m');
catch
    % If config.m is missing or causes an error, prompt the user for input
    disp('Could not load config.m. Please enter your credentials:');
    username = input('Enter your WRDS username: ', 's');
    password = input('Enter your WRDS password: ', 's');
end

WRDS = java.sql.DriverManager.getConnection(dbURL, username, password);

%% ZC Data

q = strcat('SELECT DISTINCT ON (date) ',...
                      ' date, days, rate FROM wrds.optionm_all.zerocd',...
                      ' WHERE date > ''2008-01-01''',...
                      ' AND days < 30',...
                      ' ORDER BY date, days DESC');

q = WRDS.prepareStatement(q);
rs = q.executeQuery();
MetaData = rs.getMetaData;
numCols = MetaData.getColumnCount;
dataL = cell(0,numCols);
for colIdx = numCols : -1 : 1
    ColumnNames{colIdx} = char(MetaData.getColumnLabel(colIdx));
    ColumnType{colIdx}  = char(MetaData.getColumnClassName(colIdx));
end
ColumnType = regexprep(ColumnType,'.*\.','');
 
% Loop through result set and save data into a MATLAB cell array:
rowIdx = 1;
while rs.next
    for colIdx = 1 : numCols
        switch ColumnType{colIdx}
            case {'Float','Double'}
                dataL{rowIdx,colIdx} = rs.getDouble(colIdx);
            case {'Long','Integer','Short','BigDecimal'}
                dataL{rowIdx,colIdx} = double(rs.getDouble(colIdx));
            case 'Boolean'
                dataL{rowIdx,colIdx} = logical(rs.getBoolean(colIdx));
            otherwise
                dataL{rowIdx,colIdx} = char(rs.getString(colIdx));
        end
    end
    rowIdx = rowIdx + 1;
end

q = strcat('SELECT DISTINCT ON (date) ',...
                      ' date, days, rate FROM wrds.optionm_all.zerocd',...
                      ' WHERE date > ''2008-01-01''',...
                      ' AND days > 29',...
                      ' ORDER BY date, days');

q = WRDS.prepareStatement(q);
rs = q.executeQuery();
MetaData = rs.getMetaData;
numCols = MetaData.getColumnCount;
dataU = cell(0,numCols);
for colIdx = numCols : -1 : 1
    ColumnNames{colIdx} = char(MetaData.getColumnLabel(colIdx));
    ColumnType{colIdx}  = char(MetaData.getColumnClassName(colIdx));
end
ColumnType = regexprep(ColumnType,'.*\.','');
 
% Loop through result set and save data into a MATLAB cell array:
rowIdx = 1;
while rs.next
    for colIdx = 1 : numCols
        switch ColumnType{colIdx}
            case {'Float','Double'}
                dataU{rowIdx,colIdx} = rs.getDouble(colIdx);
            case {'Long','Integer','Short','BigDecimal'}
                dataU{rowIdx,colIdx} = double(rs.getDouble(colIdx));
            case 'Boolean'
                dataU{rowIdx,colIdx} = logical(rs.getBoolean(colIdx));
            otherwise
                dataU{rowIdx,colIdx} = char(rs.getString(colIdx));
        end
    end
    rowIdx = rowIdx + 1;
end

R = [str2double(erase(string(dataU(:,1)), '-')),cell2mat(dataL(:,2)),cell2mat(dataU(:,2)),cell2mat(dataL(:,3)),cell2mat(dataU(:,3))];
         
%% Close Connection and Save

rs.close();
WRDS.close();

dataFolder = getPath('Data'); % Get the path to the Data folder

% Save the file in the Data folder
save(fullfile(dataFolder, 'ZC_2008_30d.mat'), 'R');