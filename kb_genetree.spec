/*
A KBase module: kb_genetree
*/

module kb_genetree {

        
    /*
    ** Common types
    */
    typedef string workspace_name;
    typedef int    workspace_id;
    typedef string data_obj_ref;
    typedef string data_obj_name;
    typedef int    bool;


    /*
    ** Report Results
    */
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;


    /* run_genetree_genome_context()
    **
    ** create an interactive report image of a tree browser
    */
    typedef structure {
        workspace_name  workspace_name;
        workspace_id    workspace_id;
        data_obj_ref    input_genetree_ref;
        string          genome_disp_name_config;
	float           slice_width;
	int             max_rows;
	float           prevalence_color_threshold;
    } GeneTreeGenomeContext_Input;

    funcdef run_genetree_genome_context(GeneTreeGenomeContext_Input params)
        returns (ReportResults output)
        authentication required;
};
