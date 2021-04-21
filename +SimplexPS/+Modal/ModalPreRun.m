% This .m file is used to prepare an excel file for user to configure
% participation/sensitivity analysis. If the file existed, this program
% will clear the original contents. The excel file is built based on the
% running results of the toolbox.

% Author: Yue Zhu

% Warning added by Yitong:
%
% This function may give an error if the current path of Matlab is not the
% toolbox root path, this should be improved.
%
% Running this function needs to close all excel forms. Then, the
% UserData.xlsx will be re-opened. This might influence the unsaved work of
% users.


FileModal=[cd '\Examples\ParticipationAnalysis\ModalConfig_' UserData];
SimplexPS.Modal.ExcelPrep(FileModal); %create a new excel file, or clear old contents.
% write contents in the excel file.

% *For State PF analysis, Auto-Select will select all states at each device. It will
% also select two modes lower than 100Hz, with most small damping coefficient, but not
% around 0Hz (0.1Hz tollerence) or Fbase(1Hz tollerence).
% *For Impedance PF analysis, Auto-select will select all devices for
% Layer-1 analysis, d-d axis for bode plot, and device-1 for Layer-3
% analysis. As well, it will select two modes with lowest damping.
% Write 1 to enable auto select.
AutoSel = 1;

[AutoSelResult] = SimplexPS.Modal.ExcelWrite(N_Bus,N_Device,DeviceType,...
    DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr, GminSS, GsysDSS, AutoSel, Fbase, FileModal);

fprintf('%s is now ready.\nPlease open the file and select the states and devices you are interested.\n',FileModal);
fprintf('After selection, save the excel file and run Modal Analysis.m.\n');
winopen(UserData);
if AutoSel == 0
    winopen(FileModal);
end
if AutoSel == 1 && AutoSelResult == 1
    SimplexPS.Modal.ModalAnalysis
elseif AutoSelResult == 0
elseif AutoSel ==1 && AutoSelResult == 0
    error(['Error: Mode Auto-Selection failed. Please open ModalConfig.xlsx file to select the mode manually.'])
else
end