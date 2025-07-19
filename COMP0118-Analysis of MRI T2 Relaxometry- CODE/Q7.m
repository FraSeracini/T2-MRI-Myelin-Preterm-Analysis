%% Q7
% Investigate how parameters change between subjects over all the brain and
% in specific parts
% I used in advance the NLLS algorithm for a 3 compartment model to estimate all the
% parameters and I saved the results

%% Load the data
load('all_data.mat');
load("info2.mat");
TEs=load('TEs.txt');
load("all_T2_data.mat");

%% Compute the mean for the brain, the WM and all the parcellations

num_subject = 143;

% Initialize
mean_brain = zeros(num_subject, 4);
mean_wm = zeros(num_subject, 4);
mean_frontal = zeros(num_subject, 4);
mean_parietal = zeros(num_subject, 4);
mean_temporal = zeros(num_subject, 4);
mean_occipital = zeros(num_subject, 4);
mean_insula = zeros(num_subject, 4);
mean_cingulate = zeros(num_subject, 4);
mean_corpus_callosum = zeros(num_subject, 4);
mean_internal_capsule = zeros(num_subject, 4);


% Calculate the average values for the different parameters and different
% masks
for i = 1 : 143
    if i ~= 90 && i ~= 107 && i ~= 115 && i ~= 129

        data = all_subjects_data{i};

        % Get masks for this subject
        mask = all_mask(:,:,:,i) > 0;
        wm_mask = all_par1(:,:,:,i) > 0;
        frontal_mask = all_par1(:,:,:,i) == 1 | all_par1(:,:,:,i) == 2;
        temporal_mask = all_par1(:,:,:,i) == 3 | all_par1(:,:,:,i) == 4;
        parietal_mask = all_par1(:,:,:,i) == 5 | all_par1(:,:,:,i) == 6;
        occipital_mask = all_par1(:,:,:,i) == 7 | all_par1(:,:,:,i) == 8;
        cingulate_mask = all_par1(:,:,:,i) == 9 | all_par1(:,:,:,i) == 10;
        insula_mask = all_par1(:,:,:,i) == 11 | all_par1(:,:,:,i) == 12;
        corpus_callosum_mask = all_par2(:,:,:,i) == 2 | all_par2(:,:,:,i) == 4;
        internal_capsule_mask = all_par2(:,:,:,i) == 1 | all_par2(:,:,:,i) == 3;


        % Compute means and append to results
        mean_brain(i,:) = [mean(data.S0_map(mask), 'omitnan'), ...
                   mean(data.T2_1_map(mask), 'omitnan'), ...
                   mean(data.T2_2_map(mask), 'omitnan'), ...
                   mean(data.T2_3_map(mask), 'omitnan')];
        
        mean_wm(i,:) = [mean(data.S0_map(wm_mask), 'omitnan'), ...
            mean(data.T2_1_map(wm_mask), 'omitnan'), ...
            mean(data.T2_2_map(wm_mask), 'omitnan'), ...
            mean(data.T2_3_map(wm_mask), 'omitnan')];

        mean_frontal(i,:) = [mean(data.S0_map(frontal_mask), 'omitnan'), ...
            mean(data.T2_1_map(frontal_mask), 'omitnan'), ...
            mean(data.T2_2_map(frontal_mask), 'omitnan'), ...
            mean(data.T2_3_map(frontal_mask), 'omitnan')];

        mean_temporal(i,:) = [mean(data.S0_map(temporal_mask), 'omitnan'), ...
            mean(data.T2_1_map(temporal_mask), 'omitnan'), ...
            mean(data.T2_2_map(temporal_mask), 'omitnan'), ...
            mean(data.T2_3_map(temporal_mask), 'omitnan')];

        mean_parietal(i,:) = [mean(data.S0_map(parietal_mask), 'omitnan'), ...
            mean(data.T2_1_map(parietal_mask), 'omitnan'), ...
            mean(data.T2_2_map(parietal_mask), 'omitnan'), ...
            mean(data.T2_3_map(parietal_mask), 'omitnan')];

        mean_occipital(i,:) = [mean(data.S0_map(occipital_mask), 'omitnan'), ...
            mean(data.T2_1_map(occipital_mask), 'omitnan'), ...
            mean(data.T2_2_map(occipital_mask), 'omitnan'), ...
            mean(data.T2_3_map(occipital_mask), 'omitnan')];

        mean_cingulate(i,:) = [mean(data.S0_map(cingulate_mask), 'omitnan'), ...
            mean(data.T2_1_map(cingulate_mask), 'omitnan'), ...
            mean(data.T2_2_map(cingulate_mask), 'omitnan'), ...
            mean(data.T2_3_map(cingulate_mask), 'omitnan')];

        mean_insula(i,:) = [mean(data.S0_map(insula_mask), 'omitnan'), ...
            mean(data.T2_1_map(insula_mask), 'omitnan'), ...
            mean(data.T2_2_map(insula_mask), 'omitnan'), ...
            mean(data.T2_3_map(insula_mask), 'omitnan')];

        mean_corpus_callosum(i,:) = [mean(data.S0_map(corpus_callosum_mask), 'omitnan'), ...
            mean(data.T2_1_map(corpus_callosum_mask), 'omitnan'), ...
            mean(data.T2_2_map(corpus_callosum_mask), 'omitnan'), ...
            mean(data.T2_3_map(corpus_callosum_mask), 'omitnan')];

        mean_internal_capsule(i,:) = [mean(data.S0_map(internal_capsule_mask), 'omitnan'), ...
            mean(data.T2_1_map(internal_capsule_mask), 'omitnan'), ...
            mean(data.T2_2_map(internal_capsule_mask), 'omitnan'), ...
            mean(data.T2_3_map(internal_capsule_mask), 'omitnan')];

    end
end

% Save results
save('regional_T2_means.mat', ...
     'mean_brain', 'mean_wm', 'mean_frontal', 'mean_temporal', 'mean_parietal', ...
     'mean_occipital', 'mean_cingulate', 'mean_insula', ...
     'mean_corpus_callosum', 'mean_internal_capsule');

%% Overall Brain

% Load data
load('info2.mat');                 
load('regional_T2_means.mat');     

% Define group and gender vectors 
% EPT if gestational age < 36
group = nan(size(division,1),1);
group(division(:,5) < 36) = 1;   % EPT
group(division(:,5) >= 36) = 2;  % FT

% Gender
gender = nan(size(division,1),1);
gender(division(:,6) == 77) = 1;  % Male
gender(division(:,6) == 70) = 2;  % Female

% Initialize mean/std results for 2 groups (EPT/FT) 
param_names = {'S0', 'T2_1', 'T2_2', 'T2_3'};
mean_ept = zeros(1,4);
std_ept = zeros(1,4);
mean_term = zeros(1,4);
std_term = zeros(1,4);

% Compute mean/std for EPT and FT for each parameter 
for i = 1:4
    param = mean_brain(:,i);
    mean_ept(i)  = mean(param(group == 1 & param ~= 0), 'omitnan');
    std_ept(i)   = std(param(group == 1 & param ~= 0), 'omitnan');
    mean_term(i) = mean(param(group == 2 & param ~= 0), 'omitnan');
    std_term(i)  = std(param(group == 2 & param ~= 0), 'omitnan');
end

% Format output
mean_std_ept = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_ept, std_ept, 'UniformOutput', false);
mean_std_term = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_term, std_term, 'UniformOutput', false);

% Create table for EPT vs FT
summary_table_simple = table(mean_std_ept', mean_std_term', 'VariableNames', {'EPT', 'Term'}, 'RowNames', param_names);
disp('Mean (std) for mean\_brain parameters grouped by EPT and Term:');
disp(summary_table_simple);

% Group also by gender
labels = {'EPT-M', 'EPT-F', 'FT-M', 'FT-F'};
group_code = group + (gender-1)*2;  % [1=EPT-M, 2=EPT-F, 3=FT-M, 4=FT-F]

% Initialize 
mean_4groups = strings(4,4);
for p = 1:4
    param = mean_brain(:,p);
    for g = 1:4
        values = param(group_code == g & param ~= 0);
        m = mean(values, 'omitnan');
        s = std(values, 'omitnan');
        mean_4groups(g,p) = sprintf('%.2f (%.2f)', m, s);
    end
end

% Create summary table 
summary_table_extended = array2table(mean_4groups, 'VariableNames', param_names, 'RowNames', labels);
disp('Mean (std) for mean\_brain parameters grouped by EPT/FT and Gender:');
disp(summary_table_extended);

% Boxplot for T2_1 - EPT vs FT
T2_1 = mean_brain(:,2);
valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group);
figure;
boxplot(T2_1(valid), group(valid), 'Labels', {'EPT', 'FT'}, 'Symbol', '');
title('T2\_1 in mean\_brain: EPT vs FT');
ylabel('T2\_1 (ms)');
grid on;

% Boxplot for T2_1 - EPT/FT 
valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group_code);
figure;
boxplot(T2_1(valid), group_code(valid), 'Labels', labels, 'Symbol', '');
title('T2\_1 in mean\_brain: EPT/FT by Gender');
ylabel('T2\_1 (ms)');
grid on;

% Correlations
% Get predictors
height = division(:,2);
weight = division(:,3);
iq     = division(:,4);
gest   = division(:,5);
sex    = gender;  

predictors = {height, weight, iq, gest, sex};
predictor_names = {'Height', 'Weight', 'IQ', 'GestWeeks', 'Sex'};

% Correlate each predictor with each parameter in mean_brain
correlation_results = table();
idx = 1;
for p = 1:length(predictors)
    pred = predictors{p};
    for j = 1:4
        param = mean_brain(:,j);
        valid = ~isnan(pred) & ~isnan(param) & param ~= 0;
        [r, pval] = corr(pred(valid), param(valid), 'type', 'Pearson');
        correlation_results(idx,:) = {predictor_names{p}, param_names{j}, r, pval, abs(r)};
        idx = idx + 1;
    end
end

correlation_results.Properties.VariableNames = {'Predictor', 'Parameter', 'r_value', 'p_value', 'abs_r'};
correlation_results = sortrows(correlation_results, 'abs_r', 'descend');

disp('Correlation between mean\_brain parameters and demographic predictors:');
disp(correlation_results);

%% White Matter

%% White Matter

% Load data
load('info2.mat');                 
load('regional_T2_means.mat');     

% Define group and gender vectors 
% EPT if gestational age < 36
group = nan(size(division,1),1);
group(division(:,5) < 36) = 1;   % EPT
group(division(:,5) >= 36) = 2;  % FT

% Gender
gender = nan(size(division,1),1);
gender(division(:,6) == 77) = 1;  % Male
gender(division(:,6) == 70) = 2;  % Female

% Initialize mean/std results for 2 groups (EPT/FT) 
param_names = {'S0', 'T2_1', 'T2_2', 'T2_3'};
mean_ept = zeros(1,4);
std_ept = zeros(1,4);
mean_term = zeros(1,4);
std_term = zeros(1,4);

% Compute mean/std for EPT and FT for each parameter 
for i = 1:4
    param = mean_wm(:,i);
    mean_ept(i)  = mean(param(group == 1 & param ~= 0), 'omitnan');
    std_ept(i)   = std(param(group == 1 & param ~= 0), 'omitnan');
    mean_term(i) = mean(param(group == 2 & param ~= 0), 'omitnan');
    std_term(i)  = std(param(group == 2 & param ~= 0), 'omitnan');
end

% Format output
mean_std_ept = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_ept, std_ept, 'UniformOutput', false);
mean_std_term = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_term, std_term, 'UniformOutput', false);

% Create table for EPT vs FT
summary_table_simple = table(mean_std_ept', mean_std_term', 'VariableNames', {'EPT', 'Term'}, 'RowNames', param_names);
disp('Mean (std) for mean\_wm parameters grouped by EPT and Term:');
disp(summary_table_simple);

% Group also by gender
labels = {'EPT-M', 'EPT-F', 'FT-M', 'FT-F'};
group_code = group + (gender-1)*2;  % [1=EPT-M, 2=EPT-F, 3=FT-M, 4=FT-F]

% Initialize 
mean_4groups = strings(4,4);
for p = 1:4
    param = mean_wm(:,p);
    for g = 1:4
        values = param(group_code == g & param ~= 0);
        m = mean(values, 'omitnan');
        s = std(values, 'omitnan');
        mean_4groups(g,p) = sprintf('%.2f (%.2f)', m, s);
    end
end

% Create summary table 
summary_table_extended = array2table(mean_4groups, 'VariableNames', param_names, 'RowNames', labels);
disp('Mean (std) for mean\_wm parameters grouped by EPT/FT and Gender:');
disp(summary_table_extended);

% Boxplot for T2_1 - EPT vs FT
T2_1 = mean_wm(:,2);
valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group);
figure;
boxplot(T2_1(valid), group(valid), 'Labels', {'EPT', 'FT'}, 'Symbol', '');
title('T2\_1 in mean\_wm: EPT vs FT');
ylabel('T2\_1 (ms)');
grid on;

% Boxplot for T2_1 - EPT/FT by Gender
valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group_code);
figure;
boxplot(T2_1(valid), group_code(valid), 'Labels', labels, 'Symbol', '');
title('T2\_1 in mean\_wm: EPT/FT by Gender');
ylabel('T2\_1 (ms)');
grid on;

% Correlations
% Get predictors
height = division(:,2);
weight = division(:,3);
iq     = division(:,4);
gest   = division(:,5);
sex    = gender;  

predictors = {height, weight, iq, gest, sex};
predictor_names = {'Height', 'Weight', 'IQ', 'GestWeeks', 'Sex'};

% Correlate each predictor with each parameter in mean_wm
correlation_results = table();
idx = 1;
for p = 1:length(predictors)
    pred = predictors{p};
    for j = 1:4
        param = mean_wm(:,j);
        valid = ~isnan(pred) & ~isnan(param) & param ~= 0;
        [r, pval] = corr(pred(valid), param(valid), 'type', 'Pearson');
        correlation_results(idx,:) = {predictor_names{p}, param_names{j}, r, pval, abs(r)};
        idx = idx + 1;
    end
end

correlation_results.Properties.VariableNames = {'Predictor', 'Parameter', 'r_value', 'p_value', 'abs_r'};
correlation_results = sortrows(correlation_results, 'abs_r', 'descend');

disp('Correlation between mean\_wm parameters and demographic predictors:');
disp(correlation_results);

%% WM Lobs

% Load data
load('info2.mat');                 
load('regional_T2_means.mat');     

% Define group and gender vectors 
group = nan(size(division,1),1);
group(division(:,5) < 36) = 1;   % EPT
group(division(:,5) >= 36) = 2;  % FT

gender = nan(size(division,1),1);
gender(division(:,6) == 77) = 1;  % Male
gender(division(:,6) == 70) = 2;  % Female

% Labels and setup
lobes = {
    'mean_frontal',     'Frontal';
    'mean_temporal',    'Temporal';
    'mean_occipital',   'Occipital';
    'mean_parietal',    'Parietal';
    'mean_insula',      'Insula';
    'mean_cingulate',   'Cingulate'
};
param_names = {'S0', 'T2_1', 'T2_2', 'T2_3'};
labels = {'EPT-M', 'EPT-F', 'FT-M', 'FT-F'};
group_code = group + (gender-1)*2;  % [1=EPT-M, 2=EPT-F, 3=FT-M, 4=FT-F]

% Loop over each lobe
for l = 1:size(lobes,1)
    var_name = lobes{l,1};
    region_label = lobes{l,2};
    data = eval(var_name);  % access mean_<region>

    % --- Mean and std by group (EPT vs FT) ---
    mean_ept = zeros(1,4);
    std_ept = zeros(1,4);
    mean_term = zeros(1,4);
    std_term = zeros(1,4);

    for i = 1:4
        param = data(:,i);
        mean_ept(i)  = mean(param(group == 1 & param ~= 0), 'omitnan');
        std_ept(i)   = std(param(group == 1 & param ~= 0), 'omitnan');
        mean_term(i) = mean(param(group == 2 & param ~= 0), 'omitnan');
        std_term(i)  = std(param(group == 2 & param ~= 0), 'omitnan');
    end

    mean_std_ept = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_ept, std_ept, 'UniformOutput', false);
    mean_std_term = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_term, std_term, 'UniformOutput', false);

    summary_table_simple = table(mean_std_ept', mean_std_term', 'VariableNames', {'EPT', 'Term'}, 'RowNames', param_names);
    disp(['Mean (std) for ', region_label, ' grouped by EPT and Term:']);
    disp(summary_table_simple);

    % --- Mean and std by group+gender ---
    mean_4groups = strings(4,4);
    for p = 1:4
        param = data(:,p);
        for g = 1:4
            values = param(group_code == g & param ~= 0);
            m = mean(values, 'omitnan');
            s = std(values, 'omitnan');
            mean_4groups(g,p) = sprintf('%.2f (%.2f)', m, s);
        end
    end

    summary_table_extended = array2table(mean_4groups, 'VariableNames', param_names, 'RowNames', labels);
    disp(['Mean (std) for ', region_label, ' grouped by EPT/FT and Gender:']);
    disp(summary_table_extended);

    % --- Boxplot T2_1: EPT vs FT ---
    T2_1 = data(:,2);
    valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group);
    figure;
    boxplot(T2_1(valid), group(valid), 'Labels', {'EPT', 'FT'}, 'Symbol', '');
    title(['T2\_1 in ', region_label, ': EPT vs FT']);
    ylabel('T2\_1 (ms)');
    grid on;

    % --- Boxplot T2_1: EPT/FT by gender ---
    valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group_code);
    figure;
    boxplot(T2_1(valid), group_code(valid), 'Labels', labels, 'Symbol', '');
    title(['T2\_1 in ', region_label, ': EPT/FT by Gender']);
    ylabel('T2\_1 (ms)');
    grid on;

    % --- Correlations with demographic predictors ---
    height = division(:,2);
    weight = division(:,3);
    iq     = division(:,4);
    gest   = division(:,5);
    sex    = gender;

    predictors = {height, weight, iq, gest, sex};
    predictor_names = {'Height', 'Weight', 'IQ', 'GestWeeks', 'Sex'};

    correlation_results = table();
    idx = 1;
    for p = 1:length(predictors)
        pred = predictors{p};
        for j = 1:4
            param = data(:,j);
            valid = ~isnan(pred) & ~isnan(param) & param ~= 0;
            [r, pval] = corr(pred(valid), param(valid), 'type', 'Pearson');
            correlation_results(idx,:) = {predictor_names{p}, param_names{j}, r, pval, abs(r)};
            idx = idx + 1;
        end
    end

    correlation_results.Properties.VariableNames = {'Predictor', 'Parameter', 'r_value', 'p_value', 'abs_r'};
    correlation_results = sortrows(correlation_results, 'abs_r', 'descend');
    disp(['Correlation between ', region_label, ' parameters and demographic predictors:']);
    disp(correlation_results);
end

%% Corpus Callosum and Internal Capsule

% Load data
load('info2.mat');                 
load('regional_T2_means.mat');     

% Define group and gender vectors 
group = nan(size(division,1),1);
group(division(:,5) < 36) = 1;   % EPT
group(division(:,5) >= 36) = 2;  % FT

gender = nan(size(division,1),1);
gender(division(:,6) == 77) = 1;  % Male
gender(division(:,6) == 70) = 2;  % Female

% Setup for both regions
regions = {
    'mean_corpus_callosum', 'Corpus Callosum';
    'mean_internal_capsule', 'Internal Capsule';
};

param_names = {'S0', 'T2_1', 'T2_2', 'T2_3'};
labels = {'EPT-M', 'EPT-F', 'FT-M', 'FT-F'};
group_code = group + (gender-1)*2;

% Loop through each region
for r = 1:size(regions,1)
    region_var = regions{r,1};
    region_label = regions{r,2};
    data = eval(region_var);

    % Mean/std by group (EPT vs FT)
    mean_ept = zeros(1,4);
    std_ept = zeros(1,4);
    mean_term = zeros(1,4);
    std_term = zeros(1,4);

    for i = 1:4
        param = data(:,i);
        mean_ept(i)  = mean(param(group == 1 & param ~= 0), 'omitnan');
        std_ept(i)   = std(param(group == 1 & param ~= 0), 'omitnan');
        mean_term(i) = mean(param(group == 2 & param ~= 0), 'omitnan');
        std_term(i)  = std(param(group == 2 & param ~= 0), 'omitnan');
    end

    mean_std_ept = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_ept, std_ept, 'UniformOutput', false);
    mean_std_term = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), mean_term, std_term, 'UniformOutput', false);

    summary_table_simple = table(mean_std_ept', mean_std_term', 'VariableNames', {'EPT', 'Term'}, 'RowNames', param_names);
    disp(['Mean (std) for ', region_label, ' grouped by EPT and Term:']);
    disp(summary_table_simple);

    % Mean/std by group × gender
    mean_4groups = strings(4,4);
    for p = 1:4
        param = data(:,p);
        for g = 1:4
            values = param(group_code == g & param ~= 0);
            m = mean(values, 'omitnan');
            s = std(values, 'omitnan');
            mean_4groups(g,p) = sprintf('%.2f (%.2f)', m, s);
        end
    end

    summary_table_extended = array2table(mean_4groups, 'VariableNames', param_names, 'RowNames', labels);
    disp(['Mean (std) for ', region_label, ' grouped by EPT/FT and Gender:']);
    disp(summary_table_extended);

    % Boxplot for T2_1 - EPT vs FT
    T2_1 = data(:,2);
    valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group);
    figure;
    boxplot(T2_1(valid), group(valid), 'Labels', {'EPT', 'FT'}, 'Symbol', '');
    title(['T2\_1 in ', region_label, ': EPT vs FT']);
    ylabel('T2\_1 (ms)');
    grid on;

    % Boxplot for T2_1 - EPT/FT by gender
    valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group_code);
    figure;
    boxplot(T2_1(valid), group_code(valid), 'Labels', labels, 'Symbol', '');
    title(['T2\_1 in ', region_label, ': EPT/FT by Gender']);
    ylabel('T2\_1 (ms)');
    grid on;

    % Correlations with predictors
    height = division(:,2);
    weight = division(:,3);
    iq     = division(:,4);
    gest   = division(:,5);
    sex    = gender;

    predictors = {height, weight, iq, gest, sex};
    predictor_names = {'Height', 'Weight', 'IQ', 'GestWeeks', 'Sex'};

    correlation_results = table();
    idx = 1;
    for p = 1:length(predictors)
        pred = predictors{p};
        for j = 1:4
            param = data(:,j);
            valid = ~isnan(pred) & ~isnan(param) & param ~= 0;
            [r, pval] = corr(pred(valid), param(valid), 'type', 'Pearson');
            correlation_results(idx,:) = {predictor_names{p}, param_names{j}, r, pval, abs(r)};
            idx = idx + 1;
        end
    end

    correlation_results.Properties.VariableNames = {'Predictor', 'Parameter', 'r_value', 'p_value', 'abs_r'};
    correlation_results = sortrows(correlation_results, 'abs_r', 'descend');
    disp(['Correlation between ', region_label, ' parameters and demographic predictors:']);
    disp(correlation_results);
end

%% Comparison of T2_1 values across brain regions for EPT and FT subjects. 

% Load required data
load('info2.mat');
load('regional_T2_means.mat');

% Define group vector
group = nan(size(division,1),1);
group(division(:,5) < 36) = 1;   % EPT
group(division(:,5) >= 36) = 2;  % FT

% Define region variable names and display names
region_vars = {
    'mean_brain',           'brain';
    'mean_wm',              'wm';
    'mean_frontal',         'frontal';
    'mean_temporal',        'temporal';
    'mean_parietal',        'parietal';
    'mean_occipital',       'occipital';
    'mean_cingulate',       'cingulate';
    'mean_insula',          'insula';
    'mean_corpus_callosum', 'corpus callosum';
    'mean_internal_capsule','internal capsule'
};

% Initialize output containers
n_regions = size(region_vars,1);
EPT_stats = strings(n_regions,1);
Term_stats = strings(n_regions,1);
p_values = nan(n_regions,1);

% Loop over each region
for i = 1:n_regions
    region_data = eval(region_vars{i,1});
    region_name = region_vars{i,2};

    T2_1 = region_data(:,2);  % T2_1 is column 2

    % Filter valid values
    valid = ~isnan(T2_1) & T2_1 ~= 0 & ~isnan(group);
    group_clean = group(valid);
    T2_1_clean = T2_1(valid);

    % Split by group
    ept_vals = T2_1_clean(group_clean == 1);
    term_vals = T2_1_clean(group_clean == 2);

    % Mean and std
    ept_m = mean(ept_vals, 'omitnan');
    ept_s = std(ept_vals, 'omitnan');
    term_m = mean(term_vals, 'omitnan');
    term_s = std(term_vals, 'omitnan');

    % Save formatted strings
    EPT_stats(i) = sprintf('%.2f (%.2f)', ept_m, ept_s);
    Term_stats(i) = sprintf('%.2f (%.2f)', term_m, term_s);

    % One-tailed t-test (EPT < FT)
    [~, p] = ttest2(ept_vals, term_vals, 'Vartype', 'unequal', 'Tail', 'left');
    p_values(i) = p;
end

% Create final table
T2_1_summary_table = table(EPT_stats, Term_stats, p_values, ...
    'VariableNames', {'EPT_mean_std', 'Term_mean_std', 'p_value'}, ...
    'RowNames', region_vars(:,2));

% Display the table
disp(T2_1_summary_table);