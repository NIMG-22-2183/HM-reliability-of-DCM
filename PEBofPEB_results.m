% Reliability of dynamic causal modelling of resting state magnetoencephalography 
% This use PEB approach to assess relaibltiy of posterior proablity distributions.
%  Here PEB establishes the model of each group thereby capture the model of distribution
% from which data were drawn. Then PEB of PEB assess whether any differences 
% between each distribution. 
%-------------PEB of PEB to model population--------------------------------
clc
clear all
load('TRdata.mat')
a              = [1:14]; % subjects
Test           = D(a,1);
Re_Test        = D(a,2);
M              = struct();
N              = length(a);
M.X            = ones(N,1);
field          = {'T','A','AN','H','L','J','D','CV','a','d'} ;
[PEB1]         = spm_dcm_peb(Test,M,field);
[PEB2]         = spm_dcm_peb(Re_Test,M,field);
X              = [1 1;1 0];
PEBs           = {PEB1; PEB2};
PEB3           = spm_dcm_peb(PEBs,X);
PEB            = PEB3; 

%==================plot_PEB_results========================================
Ep   = [] ; 
Cp   = [];
np   = length(PEB.Pnames); % Parameters
nc   = size(PEB.M.X,2);    % Covariates

if size(PEB.Ep,2) ~= nc
    PEB.Ep = reshape(PEB.Ep,np,nc);
end

if isvector(PEB.Cp)
    PEB.Cp = diag(PEB.Cp);
end

corr            = spm_cov2corr(PEB.Cp);
np              = length(PEB.Pnames); % Parameters
ns              = length(PEB.Snames); % Subjects
nc              = size(PEB.Ep,2);     % Covariates
effect          = 2;
effect_idx      = 1:np:(np*nc);
peb_param_idx   = effect_idx(effect) : (effect_idx(effect) + np - 1);

% Posterior means / covariance
Ep              = PEB.Ep(:,effect);
Cp              = diag(PEB.Cp);
Cp              = Cp(peb_param_idx);
T               = 0;
Pp              = 1 - spm_Ncdf(T,abs(Ep),Cp);
Ep              = Ep .* (Pp(:) > 0.95);
Cp              = Cp .* (Pp(:) > 0.95);
names           = PEB1.Pnames;

FS_labels=10; FS_ticks=10; fs_ticks=10;
figure('color','white','units','centimeters','position',[10 10 50 7],'papersize',[10 7],'filename','E.pdf')
set(gca,'fontsize',fs_ticks)
spm_plot_ci(Ep,Cp);
set(gca,'fontsize',FS_labels)
xticks([1:np]);
set(gca,'XTickLabel',names, 'fontsize',8); % 'FontWeight','bold'
xtickangle(45+45)
xlabel('Parameter','fontsize',15)
ylabel('Effect size','fontsize',15)
box off
axis tight;

%=================Create Table of Parameters===============================
% % Define your 156 names
%  names = cell(156, 1);
% % Populate names 
% for i = 1:156
%     names{i} = sprintf('Name %d', i);
% end

% Define the number of rows and columns
num_rows = 12;
num_cols = 13;

% Open the LaTeX file for writing
fid = fopen('table.tex', 'w');

% Write the beginning of the LaTeX document
fprintf(fid, '\\documentclass{article}\n');
fprintf(fid, '\\begin{document}\n');
fprintf(fid, '\\begin{center}\n');

% Write the beginning of the table environment
fprintf(fid, '\\begin{tabular}{|');
for j = 1:num_cols
    fprintf(fid, 'c|');
end
fprintf(fid, '}\n');
fprintf(fid, '\\hline\n');

% Populate the table with names
for i = 1:num_rows
    for j = 1:num_cols
        index = (i - 1) * num_cols + j;
        if index <= numel(names)
            fprintf(fid, '%s', names{index});
        end
        if j < num_cols
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' \\\\\n');
        end
    end
    fprintf(fid, '\\hline\n');
end

% Write the end of the table and the document
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\\end{document}\n');

% Close the file
fclose(fid);

% Display message
disp('LaTeX table generation complete.');

%%
% Define the number of rows and columns
num_rows = 12;
num_cols = 13;

% Create a cell array to hold the table data
table_data = cell(num_rows, num_cols);

% Populate the table data with names
for i = 1:num_rows
    for j = 1:num_cols
        index = (i - 1) * num_cols + j;
        if index <= numel(names)
            table_data{i, j} = names{index};
        else
            table_data{i, j} = ''; % Empty string for cells without names
        end
    end
end

% Create a figure to display the table
figure;
uitable('Data', table_data, 'ColumnName', {}, 'RowName', {}, ...
    'Units', 'Normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

% Set figure title
title('Names Table Visualization');
%%
% Populate the table data with names
for i = 1:num_rows
    for j = 1:num_cols
        index = (i - 1) * num_cols + j;
        if index <= numel(names)
            table_data{i, j} = names{index};
        else
            table_data{i, j} = ''; % Empty string for cells without names
        end
    end
end

% Create a figure to display the table
fig = figure('Units', 'inches', 'Position', [1.5, 0.5, 9.5, 3.15]); % Adjust size 

% Create a uitable to display the data
uitable('Parent', fig, 'Data', table_data, 'ColumnName', {}, 'RowName', {}, ...
    'Units', 'Normalized', 'Position', [0.065, 0.065, 0.9, 0.9]);

% % Set figure title
% title('Names Table Visualization');

table_data = cell(num_rows, num_cols);

% Populate the table data with names
for i = 1:num_rows
    for j = 1:num_cols
        index = (i - 1) * num_cols + j;
        if index <= numel(names)
            table_data{i, j} = names{index};
        else
            table_data{i, j} = ''; % Empty string for cells without names
        end
    end
end

% Create a figure to display the table
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 2.68]); % Adjust size as needed
box off
axis off;

% Create a uitable to display the data
uitable('Parent', fig, 'Data', table_data, 'ColumnName', {}, 'RowName', {}, ...
    'Units', 'Normalized', 'Position', [0.15, 0.065, 0.89, 0.89], ...
    'FontSize', 9); % Set FontSize property to adjust font size
% Set figure title
% Rotate the x-axis labels clockwise by 90 degrees
% xticklabel_rotate([],90,[],'FontSize',12,'FontWeight','bold');
print('table_image_without_axes.png', '-dpng', '-r300'); % Specify the resolution (300 dpi)

num_rows = 12;
num_cols = 13;
num_names = num_rows * num_cols;

% Check if the number of names matches the table size
if numel(names) ~= num_names
    error('Number of names does not match the table size.');
end

% Reshape the names into a 12x13 table
table_data = reshape(names, num_rows, num_cols);

% Display the table
disp(table_data);
%%

fid = fopen('table_data.tex', 'w');

% Write the beginning of the LaTeX document
fprintf(fid, '\\documentclass{article}\n');
fprintf(fid, '\\begin{document}\n');
fprintf(fid, '\\begin{center}\n');

% Write the beginning of the table environment
fprintf(fid, '\\begin{tabular}{|');
for j = 1:num_cols
    fprintf(fid, 'c|');
end
fprintf(fid, '}\n');
fprintf(fid, '\\hline\n');

% Populate the table with names
for i = 1:num_rows
    for j = 1:num_cols
        fprintf(fid, '%s', table_data{i, j});
        if j < num_cols
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' \\\\\n');
        end
    end
    fprintf(fid, '\\hline\n');
end

% Write the end of the table and the document
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\\end{document}\n');

% Close the file
fclose(fid);

% Publish the LaTeX file to PDF
publish('table_data.tex', 'pdf');