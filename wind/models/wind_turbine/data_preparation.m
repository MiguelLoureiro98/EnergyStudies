% This script loads the data for the NREL 5MW wind turbine.
%
% The data comprises power coefficients, tip-speed ratios, and blade
% angles of attack.
%
% The raw data is rearranged and the columns are renamed.
% All missing values are removed.
%
% Outputs: 3 arrays ready for surface fitting (TSR, beta, Cp)

%% Load the data.

pkg_path = "C:\Users\User\Desktop\Project_repos\";
%pkg_path = "J:\Projectos_e_relatorios\Project_repos\";

if(strcmp(fit_var_name, "Cp") == 1)

    filename = pkg_path + ...
               "Controlab\Controlab\data\NREL5MW_data\Cp" + ...
               "\CpVersusTSR&Pitch_WT_Perf_updated.xlsx";

    sheetname = "CpVersusTSR&Pitch-SortedByTSR";

elseif(strcmp(fit_var_name, "Ct") == 1)

    filename = pkg_path + ...
               "Controlab\Controlab\data\NREL5MW_data\Ct" + ...
               "\CtVersusTSR&Pitch_WT_Perf_updated.xlsx";

    sheetname = "CtVersusTSR&Pitch-SortedByTSR";

else

    fprintf(2, "A valid variable should be selected.\n");

end

data = readtable(filename, "Sheet", sheetname);
data = removevars(data, "Var1");


%% Rename columns.

col_list = string(zeros(1, 16));
name_list = string(zeros(1, 16));
name_list(1) = "TSR";

for i=2:17

    col_list(i-1) = strcat("Var", num2str(i));
    name_list(i) = strcat("Beta", num2str(i-7));

end

name_list = name_list(1:16);
data = renamevars(data, col_list, name_list);

%% Remove missing values.

data = data(~any(ismissing(data), 2), :);

%% Rearrange rows and columns for surface fitting.

if(strcmp(fit_var_name, "Cp") == 1)

    data_points = 182;

else

    data_points = 211;

end

%data_points = 182;
beta_values = 15;
TSR = zeros(data_points * beta_values, 1);
beta = zeros(data_points * beta_values, 1);
Cp = zeros(data_points * beta_values, 1);

for i=1:15

    TSR(data_points*(i-1)+1:data_points*i) = table2array(data(:, "TSR"));
    beta(data_points*(i-1)+1:data_points*i) = ones(data_points, 1) * (i-6);
    Cp(data_points*(i-1)+1:data_points*i) = table2array(data(:, strcat("Beta", num2str(i-6))));

end