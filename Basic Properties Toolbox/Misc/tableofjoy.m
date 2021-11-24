function [TestData,Table] = tableofjoy(data)
% Author: O.G. Steele on 28.05.18
% input variable 'data' should be in the form of .csv with 7 collumns.
   % in the .csv the first must be time, 2nd must be subtracted from the
   % remaining 5 collumns. Produces an X by 6 table.
   
   
TestData = importdata(data);
    Time_data = TestData.data(:,1);                      % first collumn unchanged
    ROI2 = TestData.data(:,3) - TestData.data(:,2);      % 2nd ROI collumn - 2 collumn
    ROI3 = TestData.data(:,4) - TestData.data(:,2);      % 3rd ROI collumn - 2 collumn
    ROI4 = TestData.data(:,5) - TestData.data(:,2);      % 4th ROI collumn - 2 collumn
    ROI5 = TestData.data(:,6) - TestData.data(:,2);      % 5th ROI collumn - 2 collumn
    ROI6 = TestData.data(:,7) - TestData.data(:,2);      % 6th ROI collumn - 2 collumn
TestData.data = [Time_data,ROI2,ROI3,ROI4,ROI5,ROI6];    % replaces the TestData.data file with above

    Time_text = TestData.textdata(:,1);                  % first collumn unchanged
    rest_text = TestData.textdata(:,3:7);                % collumns 3:7 unchanged
TestData.textdata = [Time_text,rest_text];               % replaces the TestData.textdata file with above

    Time_col = TestData.colheaders(:,1);                 % first collumn unchanged
    rest_col = TestData.colheaders(:,3:7);               % collumns 3:7 unchanged
TestData.colheaders = [Time_col,rest_col];

Table = table((TestData.data(:,1)), ...                  % Plots data into table as desired.
    (TestData.data(:,2)), ...
    (TestData.data(:,3)), ...
    (TestData.data(:,4)), ...
    (TestData.data(:,5)), ...
    (TestData.data(:,6)));
    Table.Properties.VariableNames = ({'Time','ROI1','ROI2','ROI3','ROI4','ROI5'});
end