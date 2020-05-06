def genotypeCheck(standard_data):
    """
    function to perform genotype check
    :param standard_data: an instance of StandardDataFormat with attributes query_sheet_path, raw_count_path, log2_cpm,
                          output_dir, experiment_columns
    :returns: None. However, outputs in subdir of align_counts a series of IGV screenshots and histograms corresponding to the genes in the
    """

    # verify standard_data has correct attributes and the paths to the query_sheet and raw_counts exit TODO: update these checks
    if not (hasattr(standard_data, 'standardized_query_path') and os.path.exists(standard_data.standardized_query_path)):
        sys.exit('the StandardDataFormat instance does not have attribute query_sheet_path')
    elif not (hasattr(standard_data, 'raw_count_path') and os.path.exists(standard_data.raw_count_path)):
        sys.exit('the StandardDataFormat instance does not have attribute raw_count_path')
    elif not (hasattr(standard_data, 'log2_cpm_path') and os.path.exists(standard_data.log2_cpm_path)):
        sys.exit('the StandardDataFormat instance does not have attribute log2_cpm_path')
    elif not (hasattr(standard_data, 'output_dir') and os.path.exists(standard_data.output_dir)):
        sys.exit('the StandardDataFormat instance does not have attribute output_dir')
    # TODO: if no experiment columns, currently fails -- fix
    elif not standard_data.experiment_columns:
        print('no experiment columns provided. By default, genotype, timepoint, treatment will be used. Do you wish to continue (y/n)?\n')
        answer = input()
        if answer == 'n' | answer == 'No' | answer == 'N':
            sys.exit('relaunch script with appropriate experiment columns')
        setattr(standard_data, 'experiment_columns', ['genotype', 'timepoint', 'treatment'])

    else:
        # create script command for genotype_check_histogram.R
        script_cmd = 'genotype_check_histogram.R -q {} -c {} -o {}'.format(standard_data.standardized_query_path,
                                                                           standard_data.log2_cpm_path,
                                                                           standard_data.output_dir)
        exp_column_statement = ' -e '
        for i in standard_data.experiment_columns:
            exp_column_statement = exp_column_statement + '{},'.format(i)
        exp_column_statement = exp_column_statement[:-1]
        script_cmd = script_cmd + exp_column_statement
        # execute genotype_check_histogram.R
        utils.executeSubProcess(script_cmd)