def motif_search_wrapper(analysispath, fasta_dict, selected_bin_labels, bins):
    """
    Wrapper for searching fasta sequences from selected bins using Prosite and creating a summary of results: 

    Vars:
        analysispath: (str) path to directory containing the analysis output files
        fasta_dict: (dict) dictionary containing various useful information on the fasta inputs
        selected_bin_labels: (int, list[int]) the labels of the bins to continue analysis on
        bins: (dict) containing the bin and fasta_id
    Returns:
        fasta_dict: (dict) of key(fasta_id) with values(summary dict of analysis results) for further processing  
    """
    get_prosite_db()
    if not os.path.isdir(os.path.expanduser(analysispath+'/motifs')):
        os.mkdir(analysispath+'/motifs')
        if isinstance(selected_bin_labels, list):
            for selected_bin_label in selected_bin_labels:
                if not os.path.isdir(os.path.expanduser(analysispath+'/motifs/bin_'+str(selected_bin_label))):
                    os.mkdir(analysispath+'/motifs/bin_'+str(selected_bin_label))
                print('Creating individual fasta records for: bin_'+str(selected_bin_label))
                fasta_list_from_bin = get_selected_seq(selected_bin_label, bins[selected_bin_label], fasta_dict, analysispath+'/motifs/bin_'+str(selected_bin_label))
                search_motifs_from_list(fasta_list_from_bin, selected_bin_label, analysispath+'/motifs/bin_'+str(selected_bin_label))
                print('Prosite scan completed for bin:'+str(selected_bin_label))
                print(' ')
        else: 
            selected_bin_label = selected_bin_labels
            if not os.path.isdir(os.path.expanduser(analysispath+'/motifs/bin_'+str(selected_bin_label))):
                os.mkdir(analysispath+'/motifs/bin_'+str(selected_bin_label))
            print('Creating individual fasta records for: bin_'+str(selected_bin_label))
            fasta_list_from_bin = get_selected_seq(selected_bin_label, bins[selected_bin_label], fasta_dict, analysispath+'/motifs/bin_'+str(selected_bin_label))
            for fasta in fasta_list_from_bin:
                create_new_fasta_record(fasta, fasta_dict[fasta]['sequence'], analysispath+'/motifs/bin_'+str(selected_bin_label))
            search_motifs_from_list(fasta_list_from_bin, selected_bin_label, analysispath+'/motifs/bin_'+str(selected_bin_label))
            print('Prosite scan completed for bin:'+str(selected_bin_label))
            print(' ')
        for folder in os.listdir(os.path.expanduser(analysispath+'/motifs/')):
            report_motif_results(analysispath+'/motifs/'+folder)
        print('For more information on motifs per sequence, please see relevant files in folder:'+analysispath+'/motifs/')
    return fasta_dict 