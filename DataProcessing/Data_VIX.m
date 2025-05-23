%% Option Data

clear
clc
close all

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

%% Tickers

ticker = {'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY'};

% NB: to get secid for ticker VIX:
SECID = {};
for i=1:9
    q = strcat('SELECT secid FROM wrds.optionm_all.securd WHERE ticker = ''',ticker{i},''' AND issue_type = ''%''');
    q = WRDS.prepareStatement(q);
    rs = q.executeQuery();
    MetaData = rs.getMetaData;
    numCols = MetaData.getColumnCount;
    datai = cell(0,numCols);
    for colIdx = numCols : -1 : 1
        ColumnNames{colIdx} = char(MetaData.getColumnLabel(colIdx));
        ColumnType{colIdx}  = char(MetaData.getColumnClassName(colIdx));
    end
    ColumnType = regexprep(ColumnType,'.*\.','');
    rowIdx = 1;
    while rs.next
        for colIdx = 1 : numCols
            switch ColumnType{colIdx}
                case {'Float','Double'}
                    datai{rowIdx,colIdx} = rs.getDouble(colIdx);
                case {'Long','Integer','Short','BigDecimal'}
                    datai{rowIdx,colIdx} = double(rs.getDouble(colIdx));
                case 'Boolean'
                    datai{rowIdx,colIdx} = logical(rs.getBoolean(colIdx));
                otherwise
                    datai{rowIdx,colIdx} = char(rs.getString(colIdx));
            end
        end
        rowIdx = rowIdx + 1;
    end
    SECID = [SECID;datai];
end

SECID = cell2mat(SECID);

%% GetData

dataFolder = getPath('Data'); % Get the path to the Data folder

for j=1:9
    P = [];
    data = {};
    imin = 2008;
    imax = 2020;
    for i = imin:imax
        q = strcat('SELECT secid, date, low, high, close',...
                          ' FROM wrds.optionm_all.secprd',num2str(i),...
                          ' WHERE secid = ', num2str(SECID(j)));
        q = WRDS.prepareStatement(q);
        rs = q.executeQuery();
        MetaData = rs.getMetaData;
        numCols = MetaData.getColumnCount;
        datai = cell(0,numCols);
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
                        datai{rowIdx,colIdx} = rs.getDouble(colIdx);
                    case {'Long','Integer','Short','BigDecimal'}
                        datai{rowIdx,colIdx} = double(rs.getDouble(colIdx));
                    case 'Boolean'
                        datai{rowIdx,colIdx} = logical(rs.getBoolean(colIdx));
                    otherwise
                        datai{rowIdx,colIdx} = char(rs.getString(colIdx));
                end
            end
            rowIdx = rowIdx + 1;
        end
        Pi= [cell2mat(datai(:,1)),str2double(erase( string(datai(:,2)), '-')), cell2mat(datai(:,3:end))];
        P = [P;Pi];
    end
    
    save(fullfile(dataFolder, strcat(ticker{j},'_',num2str(imin),num2str(imax))),'P');
end
       
rs.close();
WRDS.close();