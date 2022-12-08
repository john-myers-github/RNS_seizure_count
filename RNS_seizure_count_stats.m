%% RNS Data - Predict Seizure Count
% Import Data
load('RNS_time_stimparams_data.mat'); 

% Initialize Parameters and Variable Names
var_names = {'Subjects', 'time', 'charge_density', 'gamma_LH_ant', 'gamma_LH_post', 'gamma_RH_ant', 'gamma_RH_post', ...
             'theta_LH_ant', 'theta_LH_post', 'theta_RH_ant', 'theta_RH_post', 'seizure_count'};
nsubj = 12;
seizure_count  = rns_data(:,12);
time           = rns_data(:,2); 
charge_density = rns_data(:,3);
gamma_LH_ant   = rns_data(:,4);
gamma_LH_post  = rns_data(:,5);
gamma_RH_ant   = rns_data(:,6);
gamma_RH_post  = rns_data(:,7);
theta_LH_ant   = rns_data(:,8);
theta_LH_post  = rns_data(:,9);
theta_RH_ant   = rns_data(:,10);
theta_RH_post  = rns_data(:,11);

% Set Up Data Table
T = table(seizure_count, rns_data(:,2), rns_data(:,3),  rns_data(:,4), rns_data(:,5),  rns_data(:,6),  rns_data(:,7), ...
          rns_data(:,8), rns_data(:,9), rns_data(:,10), rns_data(:,11), rns_data(:,12), 'VariableNames', var_names);                            

% Run Stats 
seizure_count_time_mdl    = fitlm(T, 'seizure_count ~ time + Subjects'); 
seizure_count_cd_mdl      = fitlm(T, 'seizure_count ~ charge_density + Subjects'); 
seizure_count_LH_ant_mdl  = fitlm(T, 'seizure_count ~ time * charge_density * theta_LH_ant  * gamma_LH_ant  + Subjects');   
seizure_count_LH_post_mdl = fitlm(T, 'seizure_count ~ time * charge_density * theta_LH_post * gamma_LH_post + Subjects');                     
seizure_count_RH_ant_mdl  = fitlm(T, 'seizure_count ~ time * charge_density * theta_RH_ant  * gamma_RH_ant  + Subjects');                       
seizure_count_RH_post_mdl = fitlm(T, 'seizure_count ~ time * charge_density * theta_RH_post * gamma_RH_post + Subjects');                          

% Compute Pvals
[~, seizure_count_RH_ant_pval]  = corrcoef(seizure_count, seizure_count_RH_ant_mdl.Fitted,  'rows', 'complete'); 
[~, seizure_count_RH_post_pval] = corrcoef(seizure_count, seizure_count_RH_post_mdl.Fitted, 'rows', 'complete');  
[~, seizure_count_LH_ant_pval]  = corrcoef(seizure_count, seizure_count_LH_ant_mdl.Fitted,  'rows', 'complete'); 
[~, seizure_count_LH_post_pval] = corrcoef(seizure_count, seizure_count_LH_post_mdl.Fitted, 'rows', 'complete'); 

% Store R-squared Data
full_model_r2    = [seizure_count_LH_ant_mdl.Rsquared.Adjusted; seizure_count_LH_post_mdl.Rsquared.Adjusted; ... 
                    seizure_count_RH_ant_mdl.Rsquared.Adjusted; seizure_count_RH_post_mdl.Rsquared.Adjusted];
full_model_pvals = [seizure_count_LH_ant_pval(1,2); seizure_count_LH_post_pval(1,2);...
                    seizure_count_RH_ant_pval(1,2); seizure_count_RH_post_pval(1,2)]; 
         
%% 3D Visualization - PLANES
% SEIZURE_COUNT VISUALIZATIONS - CHARGE DENSITY
figure; 
dot_size = 300;
spectral_data = rns_data(:,11);
scatter3(spectral_data, charge_density, seizure_count, dot_size, 'filled'); 
% [X,Y,Z] = meshgrid(peak_frequencies, y_coord, z_coord);

% Fit a Plane and Plot
hold on; 
x = spectral_data; 
y = charge_density; 
z = seizure_count; 

B = [x(:) y(:) ones(size(x(:)))] \ z(:); 

xv = linspace(min(x), max(x), size(rns_data, 1))';
yv = linspace(min(y), max(y), size(rns_data, 1))';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);
s = mesh(X, Y, Z, 'FaceAlpha', 0.5);
s.FaceColor = 'Flat';
% set(gca, 'xdir', 'reverse');

xlabel('Spectral Power'); 
ylabel('Charge Density'); 
zlabel('Seizure Count'); 
colormap(jet);
view([-20.12 34.91]);
zlim([min(z), max(z)]);      
caxis([0 10]);
            
%% 3D Visualization  - PLANES
% SEIZURE_COUNT VISUALIZATIONS - TIME INTERACTIONS
figure; 
dot_size = 300;
theta_data = theta_LH_post;
gamma_data = gamma_LH_post;

scatter3(time, gamma_data, seizure_count, dot_size, 'filled'); 
% [X,Y,Z] = meshgrid(peak_frequencies, y_coord, z_coord);

% Fit a Plane and Plot
hold on; 
x = time; 
y = gamma_data; 
z = seizure_count; 

B = [x(:) y(:) ones(size(x(:)))] \ z(:); 

xv = linspace(min(x), max(x), size(rns_data, 1))';
yv = linspace(min(y), max(y), size(rns_data, 1))';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);
s = mesh(X, Y, Z, 'FaceAlpha', 0.5);
s.FaceColor = 'Flat';
% set(gca, 'xdir', 'reverse');

xlabel('Time (months)'); 
ylabel('Gamma Power'); 
zlabel('Seizure Count'); 
colormap(jet);
view([-20.12 34.91]);
zlim([min(z), max(z)]);      
% caxis([0 5]);

%% 3D Visualization  - PLANES
% SEIZURE_COUNT VISUALIZATIONS - POWER INTERACTIONS
figure; 
dot_size = 300;
theta_data = theta_LH_post;
gamma_data = gamma_LH_post;

scatter3(gamma_data, theta_data, seizure_count, dot_size, 'filled'); 
% [X,Y,Z] = meshgrid(peak_frequencies, y_coord, z_coord);

% Fit a Plane and Plot
hold on; 
x = theta_data; 
y = gamma_data; 
z = seizure_count; 

B = [x(:) y(:) ones(size(x(:)))] \ z(:); 

xv = linspace(min(x), max(x), size(rns_data, 1))';
yv = linspace(min(y), max(y), size(rns_data, 1))';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);
s = mesh(X, Y, Z, 'FaceAlpha', 0.5);
s.FaceColor = 'Flat';
% set(gca, 'xdir', 'reverse');

xlabel('Theta Power'); 
ylabel('Gamma Power'); 
zlabel('Seizure Count'); 
colormap(jet);
view([-20.12 34.91]);
zlim([min(z), max(z)]);      
caxis([0 10]);

%% 3D Visualization - SMOOTH SURFACES via LOWESS Regression
% Comparing Plots with Time vs. Charge Density

% Time x Gamma_LH_posterior
figure; 
plot(TIME_mdl_gamma_LH_post_fittedmodel, [time, gamma_LH_post], seizure_count);

% Charge Density x Gamma_LH_posterior
figure; 
plot(CD_mdl_gamma_LH_post_fittedmodel, [charge_density, gamma_LH_post], seizure_count);
       
%% Run Stats (Alternate Models) 
seizure_count_LH_ant_mdl2  = fitlm(T, 'seizure_count ~ charge_density * theta_LH_ant  * gamma_LH_ant  + Subjects');   
seizure_count_LH_post_mdl2 = fitlm(T, 'seizure_count ~ charge_density * theta_LH_post * gamma_LH_post + Subjects');                     
seizure_count_RH_ant_mdl2  = fitlm(T, 'seizure_count ~ charge_density * theta_RH_ant  * gamma_RH_ant  + Subjects');                       
seizure_count_RH_post_mdl2 = fitlm(T, 'seizure_count ~ charge_density * theta_RH_post * gamma_RH_post + Subjects');  

% Compute Pvals
[~, seizure_count_RH_ant_pval]  = corrcoef(seizure_count, seizure_count_RH_ant_mdl2.Fitted,  'rows', 'complete'); 
[~, seizure_count_RH_post_pval] = corrcoef(seizure_count, seizure_count_RH_post_mdl2.Fitted, 'rows', 'complete');  
[~, seizure_count_LH_ant_pval]  = corrcoef(seizure_count, seizure_count_LH_ant_mdl2.Fitted,  'rows', 'complete'); 
[~, seizure_count_LH_post_pval] = corrcoef(seizure_count, seizure_count_LH_post_mdl2.Fitted, 'rows', 'complete'); 

full_model_r2_2    = [seizure_count_LH_ant_mdl2.Rsquared.Adjusted; seizure_count_LH_post_mdl2.Rsquared.Adjusted; ... 
                      seizure_count_RH_ant_mdl2.Rsquared.Adjusted; seizure_count_RH_post_mdl2.Rsquared.Adjusted];
full_model_pvals_2 = [seizure_count_LH_ant_pval(1,2); seizure_count_LH_post_pval(1,2);...
                      seizure_count_RH_ant_pval(1,2); seizure_count_RH_post_pval(1,2)]; 
