%modified from Chang Liu **NOT PLOTTING PROPERLY ** MRI views aren't
%working for certain brain views
%requires ft_plot_dipole_CL() to be in fieldtrip plotting folder
%Inputs
% EEG       EEGlab data structure
% compi     component index to plot, empty = all
function plot_dipoles(EEG, compi)
%% Figure 2 - representative participant brain component
%%%%%%%%%% Standard mri
load(EEG.dipfit.mrifile);

c = (othercolor('BuOr_8')); %matlab add-on
dipole_color_FEM = c(1,:);%From red to blue

figure('color','black','position',[0 0 800 200]);
subplot(1,3,1);hold on;
pos = [0 0 0];%so that the MRI slices interset at the mid of the brain

ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [1 0 0], 'resolution', 1,'unit', 'mm');
view(90, 0)
axis equal
%ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [0 1 0], 'resolution', 1,'unit', 'mm');
%ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [0 0 1], 'resolution', 1,'unit', 'mm');

if isempty(compi)
    for i = 1:size(EEG.dipfit.model,2)
        if ~isempty(EEG.dipfit.model(i).posxyz)
           compi = [compi, i];
        end
    end
end

for i = 1:length(compi)
    IC = compi(i);
    ft_plot_dipole_CL(EEG.dipfit.model(IC).posxyz(1, :), [0 0 0],'unit', 'mm','diameter',6,'color', [252,146,114]/255);
    hold on;
    
    %saveas(gcf,fullfile(strcat(outputFolder,'Dipole_loc_'EEG.subject,'IC',num2str(IC)]));
end


subplot(1,3,2);hold on;
pos = [0 0 0];%so that the MRI slices interset at the mid of the brain

%ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [1 0 0], 'resolution', 1,'unit', 'mm');
%ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [0 1 0], 'resolution', 1,'unit', 'mm');
ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [0 0 1], 'resolution', 1,'unit', 'mm');
view(0, 90) 
axis equal
for i = 1:length(compi)
    IC = compi(i);  
    ft_plot_dipole_CL(EEG.dipfit.model(IC).posxyz(1, :), [0 0 0],'unit', 'mm','diameter',6,'color', [252,146,114]/255);
    hold on;
end
%axis equal


subplot(1,3,3);hold on;
pos = [0 0 0];%so that the MRI slices interset at the mid of the brain

%ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [1 0 0], 'resolution', 1,'unit', 'mm');
ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [0 1 0], 'resolution', 1,'unit', 'mm');
%ft_plot_slice(mri.anatomy,'transform', mri.transform,'location', pos, 'orientation', [0 0 1], 'resolution', 1,'unit', 'mm');
axis equal
view([0 -90 0])

for i = 1:length(compi)
    IC = compi(i);  
    ft_plot_dipole_CL(EEG.dipfit.model(IC).posxyz(1, :), [0 0 0],'unit', 'mm','diameter',6,'color', [252,146,114]/255);
    hold on;
end

set(gcf,'color','black')