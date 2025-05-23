%% Option Data

clear
clc
close all

O = [];

cpflag = 'C';
m = 30; %% percent moneyness

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

%% Dates

% NB: to get secid for ticker SPY:
% q = 'SELECT ticker, secid FROM wrds.optionm_all.securd s WHERE s.ticker = ''SPY''';
% rs = q.executeQuery();
% SPY = rs.getMetaData;

q = strcat('SELECT sp.date',...
                      ' FROM wrds.optionm_all.secprd2020 sp',...
                      ' WHERE sp.secid = ''109820''',...
                      ' AND sp.date BETWEEN ''2020-01-01'' AND ''2020-12-31''');
q = WRDS.prepareStatement(q);
rs = q.executeQuery();
MetaData = rs.getMetaData;
numCols = MetaData.getColumnCount;
dates = cell(0,numCols);
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
                dates{rowIdx,colIdx} = rs.getDouble(colIdx);
            case {'Long','Integer','Short','BigDecimal'}
                dates{rowIdx,colIdx} = double(rs.getDouble(colIdx));
            case 'Boolean'
                dates{rowIdx,colIdx} = logical(rs.getBoolean(colIdx));
            otherwise
                dates{rowIdx,colIdx} = char(rs.getString(colIdx));
        end
    end
    rowIdx = rowIdx + 1;
end

%% Option Price Data

n = length(dates);

for dd = 1:n
    d = string(dates(dd));
    q = strcat(' with otab(secid, date, exdate, cp_flag, strike_price, best_bid, best_offer) as',...
                        '(SELECT o.secid, o.date, o.exdate, o.cp_flag, o.strike_price,',...
                          ' o.best_bid,	o.best_offer',...
                          ' FROM wrds.optionm_all.opprcd2020 o',...
                          ' WHERE o.secid = ''109820''',...
                          ' AND o.date = ''',d,'''),',...
                    ' stab(secid, date, close) as',...
                        '(SELECT sp.secid, sp.date, sp.close',...
                          ' FROM wrds.optionm_all.secprd2020 sp',...
                          ' WHERE sp.secid = ''109820''',...
                          ' AND sp.date = ''',d,''')',...
              ' SELECT otab.secid, otab.date, otab.exdate, stab.close, otab.strike_price/1000 strike,',...
                         ' otab.best_bid,	otab.best_offer FROM otab',...
               ' INNER JOIN stab ON (stab.secid = otab.secid AND stab.date = otab.date)',...
               ' WHERE otab.cp_flag = ''',cpflag,'''');

    q = WRDS.prepareStatement(q);
    rs = q.executeQuery();
    
    % Get the column names and data types from the ResultSet's metadata
    MetaData = rs.getMetaData;
    numCols = MetaData.getColumnCount;
    data = cell(0,numCols);  % initialize
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
                    data{rowIdx,colIdx} = rs.getDouble(colIdx);
                case {'Long','Integer','Short','BigDecimal'}
                    data{rowIdx,colIdx} = double(rs.getDouble(colIdx));
                case 'Boolean'
                    data{rowIdx,colIdx} = logical(rs.getBoolean(colIdx));
                otherwise
                    data{rowIdx,colIdx} = char(rs.getString(colIdx));
            end
        end
        rowIdx = rowIdx + 1;
    end
    
    Od = [cell2mat(data(:,1)),str2double(erase( string(data(:,2)), '-')),...
                str2double(erase( string(data(:,3)), '-')), cell2mat(data(:,4:end))];
    
    indd = logical( (Od(:,4)-Od(:,5) >= -(m/100)*Od(:,4)) .* (Od(:,5)-Od(:,4) >= -(m/100)*Od(:,4)) );
    Od = Od(indd,:);
    
    Td = datenum(datetime(Od(:,3),'ConvertFrom','yyyymmdd'))-datenum(datetime(Od(:,2),'ConvertFrom','yyyymmdd'));
    Od = [Od,Td];
    
    TTd = min(( Td - 30 ).^2);
    indd = ( (Td-30).^2 == TTd );
    
    Od = Od(indd,:);
    O = [O;Od];
end

% Add time to maturity to ColumnNames
ColumnNames = [ColumnNames,{'Time to Maturity'}];

%% Close Connection and Save

rs.close();
WRDS.close();

dataFolder = getPath('Data'); % Get the path to the Data folder

% Save the file in the Data folder
save(fullfile(dataFolder, strcat('SPY_',cpflag,'_T1M_MONEY',num2str(m),'_2020')),'O');
