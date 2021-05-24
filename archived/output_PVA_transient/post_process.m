filenames = {
    '0.9_1_results.csv'    '0.9_2_results.csv'    'trans_0.9_results.csv'
    '0.1_1_results.csv'    '0.1_2_results.csv'    'trans_0.1_results.csv'
    '0.03_1_results.csv'   '0.03_2_results.csv'   'trans_0.03_results.csv'
    '0.01_1_results.csv'   '0.01_2_results.csv'   'trans_0.01_results.csv'
    '0.003_1_results.csv'  '0.003_2_results.csv'  'trans_0.003_results.csv'
    '0.0001_1_results.csv' '0.0001_2_results.csv' 'trans_0.0001_results.csv'
};

x_p = 5.85e-2;
bon = 37.78;

for i = 1:length(filenames)
    results_1 = readtable(filenames{i, 1});
    results_2 = readtable(filenames{i, 2});
    F_1 = results_1{:, 2};
    F_2 = results_2{:, 2};
    bson_1 = results_1{:, 6};
    bson_2 = results_2{:, 6};
    stress = bon*(x_p*interp1(F_2, bson_2, F_1) + (1 - x_p)*bson_1);
    results = [results_1{:, 1}, results_1{:, 2}, results_1{:, 3}, ...
        results_1{:, 4}, results_1{:, 5}, stress];
    csvwrite(filenames{i, 3}, results)
end