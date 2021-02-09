clear all
clc

global e workbook file eSheets eSheet1 GT_type

% Establishes COM server. Excel is the server and Matlab the client.
% Opens the workbook

e = actxserver('Excel.Application');  %Create a COM server running Microsoft Excel

workbook = e.Workbooks;   %Create an interface workbook
%e.Visible = 1;

% %GE LM2500+G4
% GT_type = 1;
% % % file = workbook.Open('C:\Users\riboldi\Desktop\HES-OFF\Modelling\Design\hybrid system\GT map LM2500+G4.xlsm'); %Open a file not on the MATLAB path by including the complete file specification
% %with H2
% file = workbook.Open('C:\Users\riboldi\Desktop\HES-OFF\Modelling\Design\hybrid system\GT map LM2500+G4 H2.xlsm'); %Open a file not on the MATLAB path by including the complete file specification

%GE LM6000 PF
GT_type = 2;
% file = workbook.Open('C:\Users\riboldi\Desktop\HES-OFF\Modelling\Design\hybrid system\GT map LM6000 PF.xlsx'); %Open a file not on the MATLAB path by including the complete file specification
%with H2
file = workbook.Open('C:\Users\riboldi\Desktop\HES-OFF\Modelling\Design\hybrid system\GT map LM6000 PF H2.xlsx'); %Open a file not on the MATLAB path by including the complete file specification


%Gets all the sheets of excel into eSheets
eSheets = e.ActiveWorkbook.Sheets;
