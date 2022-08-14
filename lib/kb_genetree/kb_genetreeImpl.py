# -*- coding: utf-8 -*-
#BEGIN_HEADER
from __future__ import division
from __future__ import print_function

import logging
import os
#import copy
#import math
#import re
import sys
import uuid
#import random
from datetime import datetime
from pprint import pformat,pprint
from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class kb_genetree:
    '''
    Module Name:
    kb_genetree

    Module Description:
    A KBase module: kb_genetree
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_genetree"
    GIT_COMMIT_HASH = "fda774f8b63f9d53307a95a0ea0359039e8e5230"

    #BEGIN_CLASS_HEADER
    def now_ISO(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    def log(self, target, message):
        message = '['+self.now_ISO()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_genetree_genome_context(self, ctx, params):
        """
        :param params: instance of type "GeneTreeGenomeContext_Input"
           (run_genetree_genome_context() ** ** create an interactive report
           image of a tree browser) -> structure: parameter "workspace_name"
           of type "workspace_name" (** Common types), parameter
           "input_genetree_ref" of type "data_obj_ref", parameter
           "genome_disp_name_config" of String
        :returns: instance of type "ReportResults" (** Report Results) ->
           structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_genetree_genome_context
        console = []
        self.log(console, 'Running run_genetree_genome_context() with params=')
        self.log(console, "\n" + pformat(params))


        required_params = ['workspace_name',
                           'input_genetree_ref'
                           ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError("Must define required param: '" + arg + "'")
        

        # INIT VARS
        KBase_backend = True
        KBase_filesystem = True
        Display_plot = False
        Hide_controls = True
        genome_data_format = "KBase"
        #GenomeSet_ref = '16750/58/1'  # SAGs
        #GenomeSet_ref = '16750/39/3'  # Firmicutes
        GenomeSet_ref = None
        #FeatureSet_ref = '16750/62/1'
        #FeatureSet_ref = '16750/63/1'  # rpoB Firmicutes
        FeatureSet_ref = None
        GeneTree_ref = params['input_genetree_ref']
        Search_Terms = []
        domain_data_exists = True
        domain_data_format = "KBase_domains"
        domain_family_desc_base_path = None
        tree_data_format = 'newick'
        tree_data_base_path = None
        tree_data_file = None
        PrimaryAnchor_leafId = None
        PrimaryAnchor_locusTag = None

        self.log(console,"GOT TO A")  # DEBUG

        # RUN KGB
        #%pylab notebook
        #try:  
        #    from urllib2 import urlopen
        #except ImportError:  
        #    from urllib.request import urlopen 
        #KGB_url='https://raw.githubusercontent.com/dcchivian/KGB/master/KGB.py'
        #import_code = urlopen(KGB_url)
        #%pylab notebook
        #exec(import_code.read())


        #from __future__ import print_function
        #from __future__ import division


        ###############################################################################
        # KGB user input vars  (Preferably implement as separate upstream cell)
        ###############################################################################        

        """
        ################################################################################################
        # Config #1: local notebook with NCBI annotated genbank format genomes and gene tree           #
        #            Domain Annotations: NO (KBase domain annotations require KBase genome annotation) #
        #            Gene Tree: YES                                                                    #
        #            Search Terms: YES                                                                 #
        ################################################################################################

        %pylab notebook  # must occur prior to invoking KGB if running KGB as an exec()
        KBase_backend = False

        GenomeSet_names = ["Gsulf", "DvulH", "DdesulfG20", "EcoliK12", "Bsub", "DaudaxMP104C"]
        ContigSet_names = []
        PivotFeatures_IDs = ["GSU2863", "DVU2928", "Dde_2997", "b3987", "BSU01070", "Daud0216"]
        PrimaryAnchor_leafId = "Gsulf rpoB"
        PrimaryAnchor_locusTag = "GSU2863"

        genome_annotation_system = 'NCBI'
        genome_data_format = "Genbank"
        # if relative path to scaffolds is e.g.: ./data_example/NCBI_annot/genomes/Gsulf/scaffolds/scaf_1.gbk
        gbk_ext = "gbk"
        genome_data_base_path = "./data_example/NCBI_annot/genomes"
        genome_data_extra_subpath = "/scaffolds"

        domain_data_exists = False  # NCBI annotated genome
        domain_data_format = None
        domain_data_base_path = None
        domain_data_extra_subpath = None
        domain_family_desc_base_path = None

        tree_data_format = 'newick'
        # if relative path to tree is e.g.: ../data_example/NCBI_annot/trees/rpoB_tree-names.newickfeaturesetOA
        tree_data_base_path = './data_example/NCBI_annot/trees'
        tree_data_file = 'rpoB_tree-names.newick'

        Search_Terms = ['dna-directed polymerase',
                        '16S',
                        'DNA and methyltransferase',
                        '1.10.3.-',
                        'fucI',
                        'sulfate adenylyl transferase']

        # run KGB
        from urllib.request import urlopen
        import_code = urlopen('https://raw.githubusercontent.com/dcchivian/KGB/master/KGB.py')
        exec(import_code.read())

        """


        ###############################################################################
        # KGB
        ###############################################################################
        """
        ## KGB Genome Browser (KGB)                                                  
        ##                                                                              
        ##  An IPython/Jupyter Notebook genome browser that enables comparative         
        ##  browsing and searching of genome contig assemblies, with additional         
        ##  support for                                                                 
        ##                                                                              
        ##  * domain structure visualization                                            
        ##  * gene homology                                                             
        ##  * phylogenetic trees                                                        
        ##
        ## source code available at http://github.com/dcchivian/KGB
        ##                                                                              
        ## Copyright 2015-2017 Dylan Chivian  
        ##
        ##  Initial Author: Dylan Chivian (DCChivian@lbl.gov)
        ##  $Revision: 0.3 $
        ##  $Date: 2017/05/29 00:00:00 $
        ##  $Author: dylan $
        ##
        """
        ###############################################################################
        # INIT
        ###############################################################################

        self.log(console,"GOT TO B")  # DEBUG

        # Init independent of KBase_backend
        #
        #import numpy as np # comes with pylabfeatures
        import sys  # for io and exit
        import os
        from os import walk  # for dir reading
        from os import path  # for file and dir existence check
        import re
        import json
        import csv
        #from math import pi  # get with pylab
        #from math import sqrt  # get with pylab
        from math import acos
        #from Bio import Phylo
        #import re
        #matplotlib.use('nbagg')
        import matplotlib.pyplot as pyplot
        import matplotlib.transforms as mtransforms
        import matplotlib.patheffects as path_effects
        from matplotlib.patches import Rectangle
        from matplotlib.patches import Ellipse
        from matplotlib.patches import Arc
        from matplotlib.patches import FancyBboxPatch
        #from matplotlib.lines import Line2D

        # interact is cool but outside design.  Implement via event handling on plot elements.
        #from ipywidgets import interact, interactive, fixed
        #from ipywidgets import interact
        #import ipywidgets as widgets

        self.log(console,"GOT TO C")  # DEBUG


        # tree-based row order
        #
        def feature_id_order_from_genetree (newick_string):
            ordered_feature_ids = []
            
            import ete3
            self.log(console,"INSTANTIATE TREE")  # DEBUG
            treeObj = ete3.Tree(newick_string)
            treeObj.ladderize()  # read row order from leaves?

            #self.log(console,"TRAVERSING TREE")  # DEBUG
            ##for n in treeObj.traverse('inorder'):  # inorder not viable option
            #for n in treeObj.traverse('preorder'):  # THIS IS CORRECT
            ##for n in treeObj.traverse('postorder'):  # this is wrong
            #    if n.is_leaf():
            #        ordered_feature_ids.append(n.name)

            ordered_feature_ids.extend (treeObj.get_leaf_names())  # is this the same as traverse('preorder') to get top to bottom order

            return ordered_feature_ids
        
        # Extra Init for KBase
        #
        if KBase_backend:

            print ("USING KBASE_BACKEND")  # DEBUG

            # KBase-specific library
            from installed_clients.DataFileUtilClient import DataFileUtil as DFUClient
            
            # output directory
            if KBase_filesystem:
                timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds() * 1000)
                output_dir = os.path.join(self.shared_folder, 'output_' + str(timestamp))
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                html_dir = os.path.join(output_dir, 'html')
                if not os.path.exists(html_dir):
                    os.makedirs(html_dir)

            # silence whining
            import requests
            requests.packages.urllib3.disable_warnings()

            # special char
            genome_ref_feature_id_delim = '.f:'

            # object_info tuple   
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
            WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  

            # feature location tuple
            [KB_LOC_CTG_I, KB_LOC_BEG_I, KB_LOC_STRAND_I, KB_LOC_LEN_I] = range(4)      

            """ GenomeAnnotationAPI is gone.  Replaced with methods below
            # Standard setup for accessing Data API
            import doekbase.data_api
            #from doekbase.data_api.core import ObjectAPI
            from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI
            from doekbase.data_api.sequence.assembly.api import AssemblyAPI
            from doekbase.data_api.taxonomy.taxon.api import TaxonAPI

            # loss of doekbase.data_api means need to replace
            #. ga.get_assembly()
            #. ga.get_taxon()
            #. ga.get_feature_ids()
            #. ga.get_features()
            #. tax.get_scientific_name()
            #. ass.get_contig_ids()
            #. ass.get_contig_lengths()
            #. ass.get_contig_gc_content()
            """

            # Connect to Workspace
            from installed_clients.WorkspaceClient import Workspace as workspaceService
            import os

            services = {"workspace_service_url": "https://kbase.us/services/ws/",
                        "shock_service_url": "https://kbase.us/services/shock-api/"}
            session_token = os.environ["KB_AUTH_TOKEN"]

            try:
                ws = workspaceService(services['workspace_service_url'], token=session_token)
            except Exception as e:
                raise ValueError('Unable to establish workspace connection: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            #
            # Replacement API methods
            #

            # replacement method for GenomeAnnotationAPI ga
            def gaAPI_get_genome (genome_ref=None):
                genome_ref_noVER = '/'.join(genome_ref.split('/')[0:2])
                try:
                    genome_obj = ws.get_objects([{'ref':genome_ref_noVER}])[0]
                except Exception as e:
                    raise ValueError('Unable to fetch genome object '+genome_ref+' from workspace: ' + str(e))
                return genome_obj

            # replacement method for GenomeAnnotationAPI ga.get_assembly()
            def gaAPI_get_assembly (genome_obj=None):
                genome_obj_data = genome_obj['data']
                if not genome_obj_data.get('assembly_ref'):
                    raise ValueError("Genome object missing 'assembly_ref'")
                try:
                    assembly_obj_data = ws.get_objects([{'ref':genome_obj_data['assembly_ref']}])[0]['data']
                except Exception as e:
                    raise ValueError('Unable to fetch assembly object '+assembly_ref+' from workspace: ' + str(e))
                return assembly_obj_data

            # replacement method for GenomeAnnotationAPI tax.get_scientific_name()    
            def gaAPI_get_scientific_name (genome_obj=None):
                scientific_name = None
                genome_obj_data = genome_obj['data']
                genome_obj_info = genome_obj['info']
                if not genome_obj_data.get('scientific_name') \
                    or genome_obj_data['scientific_name'].lower() == 'unknown' \
                    or genome_obj_data['scientific_name'].lower() == 'unknown_taxon':
                    scientific_name = genome_obj_info[NAME_I]
                else:
                    scientific_name = genome_obj_data['scientific_name']    
                return scientific_name

            # replacement method for GenomeAnnotationAPI ga.get_feature_ids()
            #. not a huge fan of the structure returned, but keeping for now to avoid changing a lot
            def gaAPI_get_feature_ids (genome_obj=None):
                """ OLD CLIENT CODE
                #. slicing was deathly slow.  Maybe would be faster if I implement here?

        #                feature_slice_ids = ga.get_feature_ids(group_by='region', filters={ "region_list": [{'contig_id':scaffold_id, 'strand':'?', 'start':slice_beg, 'length':slice_end-slice_beg+1}]})
                        feature_slice_ids = ga.get_feature_ids(group_by='region')

                        #"by_region": dict<str contig_id, dict<str strand, dict<string range, list<string feature_id>>>>
                        for ctg_id in feature_slice_ids['by_region'].keys():
                            if ctg_id != scaffold_id:  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                continue
                            for strand in feature_slice_ids['by_region'][ctg_id].keys():
                                for f_range in feature_slice_ids['by_region'][ctg_id][strand].keys():                    
                                    #print ("%s\t%s\t%s"%(ctg_id, strand, f_range)) # B
                                    [f_range_beg,f_range_end] = f_range.split('-')  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                    if int(f_range_beg) > slice_end or int(f_range_end) < slice_beg:
                                        continue
                                    feature_slice_ids_list.extend(feature_slice_ids['by_region'][ctg_id][strand][f_range])

                """
                feature_ids = dict()  # key is ctg_id
                [KB_LOC_CTG_I, KB_LOC_BEG_I, KB_LOC_STRAND_I, KB_LOC_LEN_I] = range(4)      
                all_features = gaAPI_get_all_features(genome_obj)

                for f_i, feature in enumerate(all_features):
                    if not feature.get('id'):
                        raise ValueError ("missing id field for feature index="+str(f_i)+" in genome "+genome_obj['info'][NAME_I])
                    fid = feature['id']
                    if not feature.get('location'):
                        raise ValueError ("missing location field for feature "+fid+" in genome "+genome_obj['info'][NAME_I])

                    ctg_id = feature['location'][0][KB_LOC_CTG_I]
                    if not feature_ids.get(ctg_id):
                        feature_ids[ctg_id] = dict()
                    strand = feature['location'][0][KB_LOC_STRAND_I]
                    if not feature_ids[ctg_id].get(strand):
                        feature_ids[ctg_id][strand] = dict()
                    if strand == '+':
                        beg = feature['location'][0][KB_LOC_BEG_I]
                        end = beg + feature['location'][0][KB_LOC_LEN_I] - 1
                        f_range = str(beg)+'-'+str(end)
                    else:
                        beg = feature['location'][0][KB_LOC_BEG_I]
                        end = beg - feature['location'][0][KB_LOC_LEN_I] + 1
                        f_range = str(end)+'-'+str(beg)
                    if not feature_ids[ctg_id][strand].get(f_range):
                        feature_ids[ctg_id][strand][f_range] = list()

                    feature_ids[ctg_id][strand][f_range].append(fid)
                return feature_ids    


            # get all features for genome
            def gaAPI_get_all_genome_features (genome_obj_data=None):
                if not genome_obj_data.get('features'):
                    raise ValueError ("missing features field for genome "+genome_obj['info'][NAME_I])
                return genome_obj_data['features']


            # get all features for AMA
            def gaAPI_get_all_AMA_features (features_handle_ref=None):
                json_features_file_path = os.path.join (output_dir, features_handle_ref+".json")
                if (not os.path.exists(json_features_file_path) \
                    or os.path.getsize(json_features_file_path) == 0):
                    try:
                        dfu = DFUClient (self.callback_url)
                    except Exception as e:
                        raise ValueError('Unable to connect to DFU: ' + str(e))
                    try:
                        dfu.shock_to_file({'handle_id': features_handle_ref,
                                           'file_path': json_features_file_path+'.gz',
                                           'unpack': 'uncompress'
                           })
                    except Exception as e:
                        raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly features from SHOCK: ' + str(e))                    

                # read file into json structure
                with open(json_features_file_path, 'r') as f:
                    features_json = json.load(f)
                return features_json

            def gaAPI_get_all_features (genome_obj=None):
                # AMA
                if genome_obj['info'][TYPE_I].startswith("KBaseMetagenomes.AnnotatedMetagenomeAssembly"):
                    features_handle_ref = genome_obj['data']['features_handle_ref']
                    all_features = gaAPI_get_all_AMA_features(features_handle_ref)
                # Genome
                else: 
                    all_features = gaAPI_get_all_genome_features(genome_obj['data'])

                return all_features
                    
            # replacement method for GenomeAnnotationAPI ga.get_features()
            def gaAPI_get_features (genome_obj=None, feature_array_type='CDS', feature_id_list=None):
                features = dict()  # key is feature id
                all_features = gaAPI_get_all_features (genome_obj)

                # capture just those feature IDs requested
                for f_i, feature in enumerate(all_features):
                    if not feature.get('id'):
                        raise ValueError ("missing id field for feature index="+str(f_i)+" in genome "+genome_obj['info'][NAME_I])
                    fid = feature['id']
                    #print ("FID: '"+fid+"'")  # DEBUG
                    if feature_id_list and fid not in feature_id_list:
                        continue
                    features[fid] = feature
                    if not feature.get('type'):
                        features[fid]['type'] = feature_array_type
                return features    

            # get contig ids from features, for when not available in parent obj
            def gaAPI_get_contig_ids_from_features (genome_obj=None):
                contig_ids_dict = dict()
                all_features = gaAPI_get_all_features (genome_obj)
                for f_i, feature in enumerate(all_features):                
                    ctg_id = feature['location'][0][KB_LOC_CTG_I]
                    contig_ids_dict[ctg_id] = True
                return sorted(contig_ids_dict.keys())
            
            # get contig lens from features, for when not available in parent obj
            def gaAPI_get_contig_lens_from_features (genome_obj=None):
                contig_lens_dict = dict()
                all_features = gaAPI_get_all_features (genome_obj)
                for f_i, feature in enumerate(all_features):                
                    ctg_id = feature['location'][0][KB_LOC_CTG_I]
                    strand = feature['location'][0][KB_LOC_STRAND_I]
                    if strand == '+':
                        beg = feature['location'][0][KB_LOC_BEG_I]
                        end = beg + feature['location'][0][KB_LOC_LEN_I] - 1
                    else:
                        beg = feature['location'][0][KB_LOC_BEG_I]
                        end = beg - feature['location'][0][KB_LOC_LEN_I] + 1
                    bigger = end
                    if beg > end:
                        bigger = beg
                    if not contig_lens_dict.get(ctg_id) or bigger > contig_lens_dict[ctg_id]:
                        contig_lens_dict[ctg_id] = bigger
                return contig_lens_dict
            
            # replacement method for GenomeAnnotationAPI ass.get_contig_ids()
            def gaAPI_get_contig_ids (genome_obj=None):
                genome_obj_data = genome_obj['data']
                contig_ids = []
                if genome_obj_data.get('contig_ids') and len(genome_obj_data['contig_ids']) > 0:
                    contig_ids = genome_obj_data['contig_ids']
                else:
                    #raise ValueError("Genome object "+genome_obj['info'][NAME_I]+" missing 'contig_ids'")
                    contig_ids = gaAPI_get_contig_ids_from_features(genome_obj)
                    
                return contig_ids

            # replacement method for GenomeAnnotationAPI ass.get_contig_lengths()
            def gaAPI_get_contig_lengths (genome_obj=None):
                contig_lengths = dict()
                genome_obj_data = genome_obj['data']
                if genome_obj_data.get('contig_ids') \
                   and len(genome_obj_data['contig_ids']) > 0 \
                   and genome_obj_data.get('contig_lengths') \
                   and len(genome_obj_data['contig_lengths']) == len(genome_obj_data['contig_ids']):
                    for contig_i,contig_id in enumerate(genome_obj_data['contig_ids']):
                        contig_lengths[contig_id] = genome_obj_data['contig_lengths'][contig_i]
                else:
                    contig_lengths = gaAPI_get_contig_lens_from_features (genome_obj)

                return contig_lengths

            # replacement method for GenomeAnnotationAPI ass.get_contig_gc_content()
            def gaAPI_get_contig_gc_content (genome_obj=None, contig_id_list=None):
                contig_gc_content = dict()
                genome_obj_data = genome_obj['data']        
                #assembly_obj_data = gaAPI_get_assembly(genome_obj=genome_obj)
                for contig_id in genome_obj_data['contig_ids']:
                    if contig_id_list and contig_id not in contig_id_list:
                        continue
                    #contig_gc_content[contig_id] = assembly_obj_data['contigs'][contig_id]['gc_content']
                    contig_gc_content[contig_id] = 0.5 # DEBUG
                return contig_gc_content


        # Init SeqIO just for non-KBase use
        #
        if KBase_backend == None or not KBase_backend:
            from Bio import SeqIO


        # Caches
        #
        search_results = []
        search_done = []
        Contig_order = []
        Contig_order_lookup = {}
        Species_name_by_genome_ref = {}
        Global_KBase_Genomes = {}
        #Global_KBase_Assemblies = {}  # accessed through Genome now
        Global_Genbank_Genomes = []
        #Global_Features = []  # not using yet
        domain_family_desc = {}
        Global_Domains = dict()


        self.log(console,"GOT TO D")  # DEBUG
        
        # Build or append to GenomeSet_names
        #
        if KBase_backend:
            PivotFeatures_IDs = []
            GenomeSet_refs = []
            GenomeSet_names = dict()

            # Use GeneTree to build GenomeSet
            #
            if GeneTree_ref != None:
                GeneTree_ref_noVER = '/'.join(GeneTree_ref.split('/')[0:2])
                try:
                    geneTree_obj = ws.get_objects([{'ref':GeneTree_ref_noVER}])
                    geneTree_data = geneTree_obj[0]['data']
                    geneTree_info = geneTree_obj[0]['info']

                except Exception as e:
                    raise ValueError('Unable to fetch geneTree object from workspace: ' + str(e))
                    #to get the full stack trace: traceback.format_exc()

                GeneTree_obj_name = geneTree_info[NAME_I]
                feature_id_row_order = feature_id_order_from_genetree (geneTree_data['tree'])
                # DEBUG
                #for fid in feature_id_row_order:
                #    self.log(console, "FEATURE_ID_FROM_NEWICK: '{}'".format(fid))
                #for leaf_id in geneTree_data['leaf_list']:
                #    self.log(console, "LEAF_ID_FROM_OBJ: '{}'".format(leaf_id))

                # build PivotFeatures and GenomeSet_Refs
                def hex2ascii(matchobj):
                    hex_val = matchobj.group(0).lstrip('%')
                    return bytearray.fromhex(hex_val).decode()

                #for leaf_id in geneTree_data['leaf_list']:
                for leaf_id in feature_id_row_order:
                    leaf_id = re.sub (r'%[\da-f][\da-f]', hex2ascii, leaf_id)
                    #print ("LEAF_ID: '"+leaf_id+"'")  # DEBUG
                    [genome_ref, f_id] = leaf_id.split(genome_ref_feature_id_delim)
                    f_id = re.sub (r'^kb\;', 'kb|', f_id)
                    PivotFeatures_IDs.append(f_id)
                    GenomeSet_refs.append(genome_ref)
                    #GenomeSet_names[genome_ref] = genome_ref  # FIX
                    #print (genome_id+" "+genome_ref)

                self.log(console,"GOT TO E")  # DEBUG

            # or use FeatureSet to build GenomeSet
            #
            elif FeatureSet_ref is not None:
                FeatureSet_ref_noVER = '/'.join(FeatureSet_ref.split('/')[0:2])
                try:
                    featureSet_obj = ws.get_objects([{'ref':FeatureSet_ref_noVER}])
                    featureSet_data = featureSet_obj[0]['data']
                    featureSet_info = featureSet_obj[0]['info']

                except Exception as e:
                    raise ValueError('Unable to fetch featureSet object from workspace: ' + str(e))
                    #to get the full stack trace: traceback.format_exc()

                featureSet_id_list = []
                if 'element_order' in featureSet_data:
                    featureSet_id_list = featureSet_data['element_order']
                else:
                    featureSet_id_list = sorted(featureSet_data['elements'].keys())

                for f_id in featureSet_id_list:
                    PivotFeatures_IDs.append(f_id)
                    if len(featureSet_data['elements'][f_id]) == 0:
                        raise ValueError("missing genome reference for feature "+f_id+" in FeatureSet "+FeatureSet_ref)
                    for genome_ref in featureSet_data['elements'][f_id]:
                        #if genome_ref not in GenomeSet_refs:
                        GenomeSet_refs.append(genome_ref)
                            #GenomeSet_names[genome_ref] = genome_ref  # FIX
                            #print (genome_id+" "+genome_ref)

            # or use GenomeSet
            #
            elif GenomeSet_ref != None:
                GenomeSet_ref_noVER = '/'.join(GenomeSet_ref.split('/')[0:2])
                try:
                    genomeSet_obj = ws.get_objects([{'ref':GenomeSet_ref_noVER}])
                    genomeSet_data = genomeSet_obj[0]['data']
                    genomeSet_info = genomeSet_obj[0]['info']

                except Exception as e:
                    raise ValueError('Unable to fetch GenomeSet object from workspace: ' + str(e))
                    #to get the full stack trace: traceback.format_exc()

                for genome_id in genomeSet_data['elements'].keys():
                    if 'ref' not in genomeSet_data['elements'][genome_id] or \
                            genomeSet_data['elements'][genome_id]['ref'] == None:
                        raise ValueError("missing reference for genome "+genome_id+" in GenomeSet "+GenomeSet_ref)
                    genome_ref = genomeSet_data['elements'][genome_id]['ref']
                    if genome_ref not in GenomeSet_refs:
                        GenomeSet_refs.append(genome_ref)
                        #GenomeSet_names[genome_ref] = genome_id
                        #print (genome_id+" "+genome_ref)

            # Need GeneTree_ref, FeatureSet_ref, or GenomeSet_ref
            #
            else:
                raise ValueError ("KGB: GeneTree_ref, FeatureSet_ref, or GenomeSet_ref is required")


        # Instantiate GenomeAnnotationAPI and Build GenomeSet_names
        #
        try:
            are_there_genome_set_names = GenomeSet_names
        except:
            GenomeSet_names = dict()
        if KBase_backend:
            for genome_ref in GenomeSet_refs:

                print ("GETTING GENOME {}".format(genome_ref))
                
                genome_obj = gaAPI_get_genome (genome_ref=genome_ref)
                Global_KBase_Genomes[genome_ref] = genome_obj
                #Global_KBase_Assemblies[genome_ref] = ass = ga.get_assembly()
                #tax = ga.get_taxon()

                #        if ws_genome_id.count('/') == 2:
                #            [ws_id, genome_id, ver] = ws_genome_id.split('/')
                #        elif ws_genome_id.count('/') == 1:
                #            [ws_id, genome_id] = ws_genome_id.split('/')
                #            ver = 'auto'
                #        else:
                #            print ("badly formatted ws_genome_id")
                #            system.exit(-1)
                #        Species_name_by_genome_id[genome_id] = tax.get_scientific_name()
                Species_name_by_genome_ref[genome_ref] = gaAPI_get_scientific_name(genome_obj=genome_obj)
                #print ("Getting Contig IDs for Genome: "+ws_genome_id+": "+Species_name_by_genome_id[genome_id])  # DEBUG

                try:
                    is_there_a_genome_set_name = GenomeSet_names[genome_ref]
                except:
                    GenomeSet_names[genome_ref] = Species_name_by_genome_ref[genome_ref]


        # Build ContigSet_names from files or from KBase object
        #
        ContigSet_names = []
        genome_contig_id_delim = '/c:'
        Genome_ref_by_Contig_name = dict()
        if KBase_backend:
            for genome_i,genome_ref in enumerate(GenomeSet_refs):
                genome_obj = Global_KBase_Genomes[genome_ref]
                if FeatureSet_ref is None and GeneTree_ref is None:
                    for scaffold_id in gaAPI_get_contig_ids(genome_obj=genome_obj):            
                        contig_name = genome_ref+genome_contig_id_delim+scaffold_id
                        ContigSet_names.append(contig_name)
                        Genome_ref_by_Contig_name[contig_name] = genome_ref
                elif GeneTree_ref is not None or FeatureSet_ref is not None:
                    fid = PivotFeatures_IDs[genome_i]
                    #print ("THIS FID: '"+fid+"'")  # DEBUG
                    these_pivot_features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=[fid])
                    # HERE
                    #print ("THESE_PIVOT_FEATURES: '"+pformat(these_pivot_features)+"'")  # DEBUG
                    scaffold_id = these_pivot_features[fid]['location'][0][KB_LOC_CTG_I]
                    contig_name = genome_ref+genome_contig_id_delim+scaffold_id
                    ContigSet_names.append(contig_name)
                    Genome_ref_by_Contig_name[contig_name] = genome_ref
                else:
                    raise ValueError ("unknown target set type in Build ContigSet_names")

        elif genome_data_format == "Genbank":
            for genome_id in GenomeSet_names:
                genome_data_dir = genome_data_base_path+'/'+genome_id+genome_data_extra_subpath
                if not path.exists(genome_data_dir):
                    print ("DIRNOTFOUND: "+genome_data_dir)
                    sys.exit(-1)
                files = []
                for (dirpath, dirnames, filenames) in walk(genome_data_dir):
                    files.extend(filenames)
                    break
                for file in files:
                    if file.endswith("."+gbk_ext):
                        scaffold_id = file[0:file.index("."+gbk_ext)]
                        contig_name = genome_id+genome_contig_id_delim+scaffold_id
                        ContigSet_names.append(contig_name)
                        Genome_ref_by_Contig_name[contig_name] = genome_id
                        print("reading "+contig_name+" ...", flush=True)
        else:
            print ("unknown genome_data_format: '"+genome_data_format+"'")
            sys.exit(0)


        # Init variables / objects
        #
        tool_title = "KGB Genome Browser"
        color_namespace_names_disp = ['Annot', 'EC']
        color_namespace_names = ['annot', 'ec']
        if domain_data_exists == None:
            domain_data_exists = False
        if domain_data_exists:
            color_namespace_names_disp.extend(['COG', 'Pfam', 'Domains'])
            color_namespace_names.extend(['cog', 'pfam', 'domains'])

        #mode_names_disp = ['Contigs', 'Genome', 'Homologs', 'Tree', 'Strains']
        #mode_names = ['contigs', 'genome', 'homologs', 'tree', strains']
        mode_names_disp = ['Contigs', 'Genome']
        mode_names = ['contigs', 'genome']    
        def_genomebrowser_mode = "contigs"
        if len(ContigSet_names) == 1:
            def_genomebrowser_mode = "genome"
        if len(PivotFeatures_IDs) != 0:
            mode_names_disp.append('Homologs')
            mode_names.append('homologs')
            def_genomebrowser_mode = "homologs"
        if tree_data_file or GeneTree_ref:
            mode_names_disp.append('Tree')
            mode_names.append('tree')
            def_genomebrowser_mode = "tree"


        # Configuration
        #
        num_genomes = len(ContigSet_names)
        def_genome_mode_n_rows = 7
        max_rows = 100
        #max_rows = 10  # DEBUG
        if max_rows < def_genome_mode_n_rows:
            max_rows = def_genome_mode_n_rows
        if def_genome_mode_n_rows%2 == 0:
            def_genome_mode_n_rows += 1
        max_feature_disp_len = 8
            
        def_popup_zorder = 100
        def_genomebrowser_window_bp_width = 1000
        def_genomebrowser_zoom_tics = 7
        def_genomebrowser_zoom = 4  # zoom values are [0..zoom_tics] -> window width = is def_genomebrowser_window_bp_width * 2**i
        def_genomebrowser_xshift = 0
        genomebrowser_window_bp_width = def_genomebrowser_window_bp_width * 2**def_genomebrowser_zoom

        def_genomebrowser_color_namespace = "annot"
        #def_genomebrowser_color_namespace = "domains"
        #def_genomebrowser_color_namespace = "cog"
        genomebrowser_mode = def_genomebrowser_mode

        genomebrowser_color_namespace = def_genomebrowser_color_namespace

        KBASE_DOMAINHIT_GENE_ID_I        = 0
        #KBASE_DOMAINHIT_GENE_BEG_I       = 1  # not used
        #KBASE_DOMAINHIT_GENE_END_I       = 2  # not used
        #KBASE_DOMAINHIT_GENE_STRAND_I    = 3  # not used
        KBASE_DOMAINHIT_GENE_HITS_DICT_I = 4
        DOMHIT_BEG_I      = 0
        DOMHIT_END_I      = 1
        DOMHIT_EVALUE_I   = 2
        DOMHIT_BITSCORE_I = 3
        DOMHIT_ALNPERC_I  = 4 
        DOMHIT_DOMFAM_I   = 5
        KB_LOC_CTG_I    = 0
        KB_LOC_BEG_I    = 1
        KB_LOC_STRAND_I = 2
        KB_LOC_LEN_I    = 3
        #KB_LOC_CTG_I = 'contig_id'
        #KB_LOC_BEG_I = 'start'
        #KB_LOC_STR_I = 'strand'
        #KB_LOC_LEN_I = 'length'


        # State
        #
        Global_State = { \
                        "ContigSet_names":                   ContigSet_names, \
                        "PivotFeatures_IDs":                 PivotFeatures_IDs, \
                        "Dataset_names_list":                ["none"], \
                        "Dataset_species_list":              [], \
                        "PrimaryContig_GCavg":               0.0, \
                        "PrimaryContig_len":                 0, \
                        "Contig_lens":                       [], \
                        "PrimaryAnchor_leafId":              PrimaryAnchor_leafId, \
                        "PrimaryAnchor_pivot_pos":           0.0, \
                        "pivot_pos_list":                    [], \
                        "genome_mode_n_rows":                def_genome_mode_n_rows, \
                        "popup_zorder":                      def_popup_zorder, \
                        "genomebrowser_pivot_pos":           0.0, \
                        "def_genomebrowser_window_bp_width": def_genomebrowser_window_bp_width, \
                        "genomebrowser_window_bp_width":     def_genomebrowser_window_bp_width * 2**def_genomebrowser_zoom, \
                        "genomebrowser_zoom":                def_genomebrowser_zoom, \
                        "genomebrowser_xshift":              def_genomebrowser_xshift, \
                        "genomebrowser_mode":                def_genomebrowser_mode, \
                        "genomebrowser_color_namespace":     def_genomebrowser_color_namespace, \
                        "function_abundance_counts":         {}
                       }

        self.log(console,"GOT TO F")  # DEBUG

        # Set up canvas and arrow config
        #
        total_rows = num_genomes
        if def_genome_mode_n_rows > total_rows:
            total_rows = def_genome_mode_n_rows
        if total_rows > max_rows:
            total_rows = max_rows
        #figure_width = 12.0
        figure_width = 11.0
        figure_height_scaling = 0.75
        top_nav_height = 1.5
        #top_margin = 1.0/total_rows
        top_margin = 0.0
        left_margin = 0.05
        right_margin = 0.05
        row_delta = 1.0/total_rows
        #row_delta = 1.0/(total_rows+1)
        bottom_margin=0.75*row_delta
        arrow_w = 0.25/total_rows
        head_w = arrow_w
        #head_l = 0.4*head_w
        head_l = 0.015
        base_arrow_label_fontsize = 10
        foi_arrow_label_fontsize = 13
        text_padding = 0.03/total_rows
        text_yshift = text_padding+0.5*arrow_w
        #arrow_label_scaling = 125.0  # about appropriate value for fontsize 12 and figure width 15 (or maybe 20)
        #arrow_label_scaling = 95.0  # appropriate value for fontsize 10 and figure width 12
        arrow_label_scaling = 100.0  # appropriate value for fontsize 10 and figure width 12
        text_disp_window_bp_limit = 80000


        # Colors
        #
        nav_bg_color = "slategray"
        popup_bg_color = "aliceblue"
        #popup_bg_color = "cornsilk"
        #popup_bg_color = "lemonchiffon"
        search_color_names = [
            'magenta', \
            'orange', \
            'darkcyan', \
            'green', \
            'red', \
            'lightskyblue', \
            'crimson'
        ]
        # feature colors names list (left out excessively light colors)
        color_names = [
            #'aliceblue',
            'aqua',
            'aquamarine',
            #'azure',
            #'beige',
            #'bisque',
            #'blanchedalmond',
            'blue',
            'blueviolet',
            'brown',
            'burlywood',
            'cadetblue',
            'chartreuse',
            'chocolate',
            'coral',
            'cornflowerblue',
            #'cornsilk',
            'crimson',
            'cyan',
            'darkblue',
            'darkcyan',
            'darkgoldenrod',
            'darkgreen',
            'darkkhaki',
            'darkmagenta',
            'darkolivegreen',
            'darkorange',
            'darkorchid',
            'darkred',
            'darksalmon',
            'darkseagreen',
            'darkslateblue',
            #'darkslategray',
            'darkturquoise',
            'darkviolet',
            'deeppink',
            'deepskyblue',
            'dodgerblue',
            'firebrick',
            'forestgreen',
            'fuchsia',
            #'gainsboro',
            'gold',
            'goldenrod',
            'green',
            'greenyellow',
            #'honeydew',
            'hotpink',
            'indianred',
            'indigo',
            'khaki',
            #'lavender',
            #'lavenderblush',
            'lawngreen',
            #'lemonchiffon',
            'lightblue',
            #'lightcoral',
            #'lightcyan'
            #'lightgoldenrodyellow',
            'lightgreen',
            'lightpink',
            'lightsalmon',
            'lightseagreen',
            'lightskyblue',
            #'lightslategray',
            #'lightsteelblue',
            #'lightyellow',
            'lime',
            'limegreen',
            'magenta',
            'maroon',
            'mediumaquamarine',
            'mediumblue',
            'mediumorchid',
            'mediumpurple',
            'mediumseagreen',
            'mediumslateblue',
            'mediumspringgreen',
            'mediumturquoise',
            'mediumvioletred',
            'midnightblue',
            #'mintcream',
            #'mistyrose',
            #'moccasin',
            'navy',
            #'oldlace',
            'olive',
            'olivedrab',
            'orange',
            'orangered',
            'orchid',
            #'palegoldenrod',
            'palegreen',
            'paleturquoise',
            'palevioletred',
            #'papayawhip',
            'peachpuff',
            #'peru',
            'pink',
            'plum',
            'powderblue',
            'purple',
            'red',
            'rosybrown',
            'royalblue',
            'saddlebrown',
            'salmon',
            'sandybrown',
            'seagreen',
            #'seashell',
            'sienna',
            'skyblue',
            'slateblue',
            'springgreen',
            'steelblue',
            #'tan',
            'teal',
            #'thistle',
            'tomato',
            'turquoise',
            'violet',
            #'wheat',
            #'yellow',
            #'yellowgreen'
            ]

        ###############################################################################
        # Subs
        ###############################################################################

        # Utilities
        #
        def compute_GC (dna_seq):
            seq_len = len(dna_seq)
            GC_cnt = 0
            AT_cnt = 0
            other_cnt = 0
            dna_seq.upper()
            GC_cnt = dna_seq.count('G') + dna_seq.count('C')
            AT_cnt = dna_seq.count('A') + dna_seq.count('T') + dna_seq.count('U')
            other_cnt = seq_len - GC_cnt - AT_cnt
            seq_len -= other_cnt
            return GC_cnt / seq_len


        # KBase feature structure
        #
        def build_feature_rec_kbase (f, f_type='CDS', source_species='', contig_i=0, dna_seq='N'):
            feature_rec = {}
            functions = []
            id_delim = '.'

            # IDs and names
            gene_name = None
            name = f['id']
            name_split = name.split(id_delim)
            if len(name_split) > 2:
                name = id_delim.join(name_split[len(name_split)-2:])
                
            locus_tag = f['id']
            aliases = []
            if f.get('aliases'):
                #for alias in f['aliases'].keys():
                for alias in f['aliases']:
                    aliases.append(alias)
                    if isinstance(alias,str):
                        if locus_tag == f['id'] and 'IPR' not in alias:
                            locus_tag = alias  # fix this to match regexp \D+_?\d+ (but not IPR*), and stop assignment
                    elif isinstance(alias,tuple):
                        if alias[0] == 'gene':
                            name = alias[1]  # this is where we want the name from!
                            gene_name = alias[1]
                        if locus_tag == f['id'] and 'IPR' not in alias[0]:
                            locus_tag = alias[0]  # fix this to match regexp \D+_?\d+ (but not IPR*), and stop assignment
                    elif isinstance(alias,list):
                        if alias[0] == 'gene':
                            name = alias[1]  # this is where we want the name from!
                            gene_name = alias[1]
                    else:
                        raise ValueError ("unknown alias type "+str(alias))

            feature_ID = f['id']
            #print ("FEATURE_ID: '"+feature_ID+"'")  # DEBUG

            # coords
            strand = f['location'][0][KB_LOC_STRAND_I]
            f_len = f['location'][0][KB_LOC_LEN_I]
            if strand == '+':
                beg = f['location'][0][KB_LOC_BEG_I]
                end = beg + f_len - 1
            else:
                end = f['location'][0][KB_LOC_BEG_I]
                beg = end - f_len + 1    
            if end < beg:
                print ("WARNING: reversed gene: %s: %s (%d-%d)"%(source_species, locus_tag, beg, end))
                tmp_pos = end
                end = beg
                beg = tmp_pos

            # annotation
            if f.get('functions'):
                if isinstance(f['functions'],list):
                    annotation = "; ".join(f['functions'])
                    functions = f['functions']
                else:
                    annotation = f['functions']
                    functions = [f['functions']]
            elif f.get('function'):
                if isinstance(f['function'],list):
                    annotation = "; ".join(f['function'])
                    functions = f['function']
                else:
                    annotation = f['function']           
                    functions = [f['function']]
            else:
                annotation = ''
            EC_in_annotation = ''
            in_paren = False
            last_good_char = 0
            close_paren_pos = 0
            for i,c in enumerate(annotation[::-1]):
                if in_paren == False and c != ' ' and c != ')':
                    break
                if in_paren == False and c != ' ' and c == ')':
                    in_paren = True
                    close_paren_pos = i
                if in_paren and c == '(':
                    last_good_char = i+1
                    break
            if last_good_char > 0:
                candidate_EC = annotation[-1-last_good_char+2:-1-close_paren_pos]
                if candidate_EC[0:2].upper() == 'EC':
                    #print (name + " EC: '" + candidate_EC +"'")  # DEBUG
                    EC_in_annotation = candidate_EC[3:]

                annotation = annotation[0:-1-last_good_char]
            annotation = annotation.strip()
            annotation = annotation.lstrip()

            # add domain hits to annotation
            try:
                domain_hits = Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[contig_i]]][ContigSet_names[contig_i]][feature_ID]
                domfam_seen = {}
                if annotation != "":
                    annotation += "\n"
                for domhit in domain_hits:  # already reverse sorted by bitscore
                    domfam = domhit[DOMHIT_DOMFAM_I]
                    try:
                        seen = domfam_seen[domfam]
                        continue
                    except:
                        domfam_seen[domfam] = True
                        if domfam[0:3] == 'COG' or domfam[0:2] == 'PF' or domfam[0:4] == 'TIGR':
                            if domfam[0:2] == 'PF':
                                dom_namespace = 'PF'
                                [domfam, ver] = domfam.split('.')
                            elif domfam[0:3] == 'COG':
                                dom_namespace = 'COG'
                            elif domfam[0:4] == 'TIGR':
                                dom_namespace = 'TIGR'
                            line_break = "\n"
                            if annotation == "":
                                line_break = ""
                            try:
                                domfam_desc = ": "+domain_family_desc[dom_namespace][domfam]
                            except:
                                domfam_desc = ""
                            annotation += line_break+domfam+domfam_desc
            except:
                domain_hits = []

            # pull out E.C.    
            EC_number = ''
            if EC_in_annotation != '':
                EC_number = EC_in_annotation

            # count functions to determine winner when picking color
            for fxn in functions:
                if fxn not in Global_State["function_abundance_counts"]:
                    Global_State["function_abundance_counts"][fxn] = 0
                Global_State["function_abundance_counts"][fxn] += 1


            # create feature_rec
            feature_rec = {"source_species": source_species, \
                           "gene_name": gene_name,
                           "name": name, \
                           "locus_tag": locus_tag, \
                           "ID": feature_ID, # will probably need to make more unique \
                           "aliases": aliases, \
                           "type": f_type, \
                           "beg_pos": beg, \
                           "end_pos": end, \
                           "strand": strand, \
                           "annot": annotation, \
                           "functions": functions, \
                           "EC_number": EC_number, \
                           "dna_seq": dna_seq
                          }                         

            return feature_rec


        # Genbank feature structure
        #
        def build_feature_rec_genbank (f, f_type='CDS', source_species='', contig_i=0, dna_seq=None):
            feature_rec = {}

            id_delim = '.'

            # IDs and names
            if "gene" in f.qualifiers:
                name = f.qualifiers['gene'][0]
            elif "locus_tag" in f.qualifiers:
                name = f.qualifiers['locus_tag'][0]
            elif "db_xref" in f.qualifiers:
                name = f.qualifiers['db_xref'][0]
            else:
                print ("strangely formatted record.  SKIPPING")
                print (f)
                return {}
        #    name_split = name.split(id_delim)
        #    name = id_delim.join(name_split[1:])

            if "locus_tag" in f.qualifiers:
                locus_tag = f.qualifiers['locus_tag'][0]
            elif "gene" in f.qualifiers:
                locus_tag = f.qualifiers['gene'][0]
            elif "db_xref" in f.qualifiers:
                locus_tag = f.qualifiers['db_xref'][0]
            else:
                print ("strangely formatted record.  SKIPPING")
                print (f)
                return {}
            locus_tag_split = locus_tag.split(id_delim)
            locus_tag = id_delim.join(locus_tag_split[1:])

            # capture aliases
            aliases = []
            if "db_xref" in f.qualifiers:
                for alias in f.qualifiers['db_xref']:
                    aliases.append(alias)

            # For now, feature ID will be locus_tag
            feature_ID = locus_tag

            # coords
            beg = f.location.start + 1  # BioPython is cute and shifts just the begin pos by -1 
            end = f.location.end
            if end < beg:
                print ("WARNING: reversed gene: %s: %s (%d-%d)"%(source_species, locus_tag, beg, end))
                tmp_pos = end
                end = beg
                beg = tmp_pos
            if f.location.strand == -1:
                strand = '-'
            else :
                strand = '+'

            # annotation
            annotation = ''
            EC_in_annotation = ''
            if f_type == "CDS" or f_type == "rRNA" or f_type == "tRNA":    
                if genome_annotation_system == 'KBase':
                    annotation = f.qualifiers['function'][0]
                else:
                    annotation = f.qualifiers['product'][0]
                in_paren = False
                last_good_char = 0
                close_paren_pos = 0
                for i,c in enumerate(annotation[::-1]):
                    if in_paren == False and c != ' ' and c != ')':
                        break
                    if in_paren == False and c != ' ' and c == ')':
                        in_paren = True
                        close_paren_pos = i
                    if in_paren and c == '(':
                        last_good_char = i+1
                        break
                if last_good_char > 0:
                    candidate_EC = annotation[-1-last_good_char+2:-1-close_paren_pos]
                    if candidate_EC[0:2].upper() == 'EC':
                        #print (name + " EC: '" + candidate_EC +"'")  # DEBUG
                        EC_in_annotation = candidate_EC[3:]

                    annotation = annotation[0:-1-last_good_char]
                annotation = annotation.strip()
                annotation = annotation.lstrip()

            # add domain hits to annotation
            try:
                domain_hits = Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[contig_i]]][feature_ID]
                domfam_seen = {}
                if annotation != "":
                    annotation += "\n"
                for domhit in domain_hits:  # already reverse sorted by bitscore
                    domfam = domhit[DOMHIT_DOMFAM_I]
                    try:
                        seen = domfam_seen[domfam]
                        continue
                    except:
                        domfam_seen[domfam] = True
                        if domfam[0:3] == 'COG' or domfam[0:2] == 'PF' or domfam[0:4] == 'TIGR':
                            if domfam[0:2] == 'PF':
                                dom_namespace = 'PF'
                                [domfam, ver] = domfam.split('.')
                            elif domfam[0:3] == 'COG':
                                dom_namespace = 'COG'
                            elif domfam[0:4] == 'TIGR':
                                dom_namespace = 'TIGR'
                            line_break = "\n"
                            if annotation == "":
                                line_break = ""
                            try:
                                domfam_desc = ": "+domain_family_desc[dom_namespace][domfam]
                            except:
                                domfam_desc = ""
                            annotation += line_break+domfam+domfam_desc
            except:
                domain_hits = []

            # pull out E.C.    
            EC_number = ''
            if EC_in_annotation != '':
                EC_number = EC_in_annotation
            else:
                if f_type == "CDS" and 'EC_number' in f.qualifiers:
                    EC_number = f.qualifiers['EC_number'][0]

            # create feature_rec
            feature_rec = {"source_species": source_species, \
                           "name": name, \
                           "locus_tag": locus_tag, \
                           "ID": feature_ID, # will probably need to make more unique \
                           "aliases": aliases, \
                           "type": f_type, \
                           "beg_pos": beg, \
                           "end_pos": end, \
                           "strand": strand, \
                           "annot": annotation, \
                           "EC_number": EC_number, \
                           "dna_seq": dna_seq
                          }                         

            #if feature_rec['name'] == "DVU2932" or feature_rec['name'] == "Daud0298":
            #    arrow_color = color_names[sum([ord(c) for c in feature_rec['annot']]) % len(color_names)]
            #    print (feature_rec['source_species'] + " " + feature_rec['name']+" "+arrow_color)

            return feature_rec


        # Load domain family descriptions
        #
        def read_domain_family_desc ():
            #domain_family_desc = {}

            # read Pfam descs
            pfam_desc_file = domain_family_desc_base_path+'/Pfam-A.clans.tsv'
            if path.isfile(pfam_desc_file):
                domain_family_desc['PF'] = {}
                with open(pfam_desc_file,'r') as f:
                    reader=csv.reader(f,delimiter='\t')
                    for pfam,clan,clan_name,pfam_name,pfam_desc in reader:
                        domain_family_desc['PF'][pfam] = pfam_desc

            # read COG descs
            COG_desc_file = domain_family_desc_base_path+'/COG_2014.tsv'
            if path.isfile(COG_desc_file):
                domain_family_desc['COG'] = {}
                with open(COG_desc_file,'r') as f:
                    reader=csv.reader(f,delimiter='\t')
                    for COG,COG_funcat,COG_desc in reader:
                        domain_family_desc['COG'][COG] = COG_desc

            # read TIGR descs
            TIGR_desc_file = domain_family_desc_base_path+'/TIGRInfo.tsv'
            if path.isfile(TIGR_desc_file):
                domain_family_desc['TIGR'] = {}
                with open(TIGR_desc_file,'r') as f:
                    reader=csv.reader(f,delimiter='\t')
                    for id,tigrId,type,roleId,geneSymbol,ec,definition in reader:
                        domain_family_desc['TIGR'][tigrId] = definition
                        if ec != "":
                            domain_family_desc['TIGR'][tigrId] += " (EC "+ec+")"

            #return domain_family_desc


        # Check search terms
        #
        def search_term_match (feature_rec):    
            term_match = []
            for i,term in enumerate(Search_Terms):
                term_uc = term.upper()  # FIX: should sanitize search terms once at beginning
                term_match.append(False)

                for k in ['name', 'locus_tag', 'ID', 'aliases', 'annot', 'EC_number']:
                    info_uc_list = []

                    if k == 'aliases':
                        try:
                            for f in feature_rec[k]:
                                if isinstance(f,str):
                                    info_uc_list.append(f.upper())
                                elif isinstance(f,tuple):
                                    info_uc_list.append(f[0].upper())
                                elif isinstance(f,list):
                                    info_uc_list.append(f[1].upper())
                        except:
                            continue                
                    else:
                        #print ("CHECKING %s %s %s %s"%(term,k,feature_rec['name'], feature_rec[k]))  # DEBUG
                        try:
                            info_uc_list = [feature_rec[k].upper()]
                        except:
                            continue

                    search_words = term_uc.split(" ")  # FIX: should sanitize search terms once at beginning
                    logical = 'AND'
                    if " OR " in term_uc:
                        logical = 'OR'

                    this_info_match = False

                    for info_uc in info_uc_list:
                        if logical == 'AND':
                            this_info_match = True
                            for word in search_words:
                                #if word == 'AND' or word == 'OR':
                                if word == 'AND' or word == 'OR' or word == '':  # FIX: should sanitize search terms once at beginning
                                    continue
                                if word not in info_uc:
                                    this_info_match = False
                                    break
                        if logical == 'OR':
                            this_info_match = False
                            for word in search_words:
                                #if word == 'AND' or word == 'OR':
                                if word == 'AND' or word == 'OR' or word == '':  # FIX: should sanitize search terms once at beginning
                                    continue
                                if word in info_uc:
                                    this_info_match = True
                                    break

                        if this_info_match:
                            term_match[i] = True
                            break

                    if this_info_match:
                        term_match[i] = True
                        break

            return term_match


        # sort rules
        #
        def sort_by_beg_pos_key (feature):
            return feature['beg_pos']    

        def sort_by_bitscore_key (domhit):
            return domhit[DOMHIT_BITSCORE_I]    

        def sort_by_domhit_len (domhit):
            return domhit[DOMHIT_END_I]-domhit[DOMHIT_BEG_I]+1

        # Read KBase format domain hits
        #
        def getDomainHits (ContigSet_names, \
                           genomebrowser_mode=def_genomebrowser_mode,
                           domain_data_format=domain_data_format):

            #print ("GETTING DOMAINHITS")  # DEBUG

            for i,contig_name in enumerate(ContigSet_names):
                (genome_id,scaffold_id) = contig_name.split(genome_contig_id_delim)
                #print("GETTING DOMAINHITS FOR "+contig_name+" "+genome_id+" c:"+scaffold_id)  # DEBUG

                if Global_State['genomebrowser_mode'] == "genome" and i > 0:
                    break

                try:
                    Domain_Hits = Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[i]]]
                except:
                    Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[i]]] = dict()

                    if KBase_backend:
                        genome_ref = Genome_ref_by_Contig_name[contig_name]
                        domain_data = None
                        found_domain_data = False

                        #print ("CONTIG_NAME: '"+contig_name+"' GENOME_REF: '"+genome_ref+"'")  # DEBUG

                        if genome_ref.count('/') == 2:
                            (ws_id, ws_genome_id, ver) = genome_ref.split('/')
                        elif genome_ref.count('/') == 1:
                            (ws_id, ws_genome_id) = genome_ref.split('/')
                            ver = 'auto'

                        try:
                            domain_annotation_obj_list = ws.list_objects({'ids':[ws_id],'type':"KBaseGeneFamilies.DomainAnnotation"})
                        except Exception as e:
                            raise ValueError('Unable to list DomainAnnotation objects from workspace: ' + str(e))
                            #to get the full stack trace: traceback.format_exc()
                        for info in domain_annotation_obj_list:
                            # Object Info Contents
                            # 0 - obj_id objid
                            # 1 - obj_name name
                            # 2 - type_string type
                            # 3 - timestamp save_date
                            # 4 - int version
                            # 5 - username saved_by
                            # 6 - ws_id wsid
                            # 7 - ws_name workspace
                            # 8 - string chsum
                            # 9 - int size 
                            # 10 - usermeta meta
                            # absolute ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
                            domain_annotation_ref = str(info[6])+'/'+str(info[0])+'/'+str(info[4])
                            try:
                                domain_data = ws.get_objects([{'ref':domain_annotation_ref}])[0]['data']  
                            except Exception as e:
                                raise ValueError('Unable to fetch DomainAnnotation object from workspace: ' + str(e))
                            #to get the full stack trace: traceback.format_exc()

                            # we found the correct DomainAnnotation object
                            if domain_data['genome_ref'] == genome_ref:
                                #print ("FOUND DOM ANNOT "+str(domain_annotation_ref)+" FOR GENOME "+str(genome_ref))  # DEBUG
                                found_domain_data = True
                                break

                        if not found_domain_data:
                            continue

                        for scaffold_id_iter in domain_data['data'].keys():
                            contig_id_iter = genome_ref+genome_contig_id_delim+scaffold_id_iter
                            Global_Domains[genome_ref][contig_id_iter] = dict()

                            for CDS_domain_list in domain_data['data'][scaffold_id_iter]:
                                gene_ID   = CDS_domain_list[KBASE_DOMAINHIT_GENE_ID_I]
                                #gene_name = re.sub ('^'+genome_object_name+'.', '', gene_ID) 
                                gene_name = gene_ID
                                #(contig_name, gene_name) = (gene_ID[0:gene_ID.index(".")], gene_ID[gene_ID.index(".")+1:])
                                #print ("DOMAIN_HIT: "+contig_name+" "+gene_name)  # DEBUG
                                #print ("DOMAIN_HIT for gene: "+gene_name)  # DEBUG
                                #gene_beg       = CDS_domain_list[KBASE_DOMAINHIT_GENE_BEG_I]
                                #gene_end       = CDS_domain_list[KBASE_DOMAINHIT_GENE_END_I]
                                #gene_strand    = CDS_domain_list[KBASE_DOMAINHIT_GENE_STRAND_I]
                                gene_hits_dict = CDS_domain_list[KBASE_DOMAINHIT_GENE_HITS_DICT_I]
                                gene_hits_list = []
                                for domfam in gene_hits_dict.keys():
                                    # skip CD hits for now
                                    if domfam[0:2] != 'PF' and domfam[0:3] != 'COG' and domfam[0:4] != 'TIGR':
                                        continue
                                    #Global_Domains[i][gene_name] = gene_hits_dict
                                    for hit in gene_hits_dict[domfam]:
                                        list_format_hit = hit
                                        list_format_hit.append (domfam)
                                        #list_format_hit[DOMHIT_BEG_I]      = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_BEG_J]
                                        #list_format_hit[DOMHIT_END_I]      = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_END_J]
                                        #list_format_hit[DOMHIT_EVALUE_I]   = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_EVALUE_J]
                                        #list_format_hit[DOMHIT_BITSCORE_I] = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_BITSCORE_J]
                                        #list_format_hit[DOMHIT_ALNPERC_I]  = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_ALNPERC_J]
                                        gene_hits_list.append(list_format_hit)
                                        #print ("   DOMAIN_HIT: "+domfam)  # DEBUG
                                        #print ("%s\t%s\t%s\t%s\t%d\t%d\t%f"%(genome_id, scaffold_id, gene_name, domfam, hit_beg, hit_end, hit_evalue))
                                Global_Domains[genome_ref][contig_id_iter][gene_name] = sorted (gene_hits_list, key=sort_by_bitscore_key, reverse=True)


                    # Or running outside KBase
                    #
                    elif domain_data_format == "KBase_domains" and domain_data_exists:

                        domain_data_path = domain_data_base_path+'/'+genome_id+domain_data_extra_subpath+'/'+genome_id+"_Domain_annot"+'.json'
                        print ("reading "+domain_data_path+" ...")
                        with open(domain_data_path, 'r') as domain_file_handle:
                            kbase_domains = json.load(domain_file_handle)

                            for scaffold_id_iter in kbase_domains['data'].keys(): 
                                contig_id_iter = Genome_ref_by_Contig_name[ContigSet_names[i]]+genome_contig_id_delim+scaffold_id_iter
                                Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[i]]][contig_id_iter] = dict()
                                for CDS_domain_list in kbase_domains['data'][scaffold_id]:
                                    gene_ID        = CDS_domain_list[KBASE_DOMAINHIT_GENE_ID_I]
                                    (contig_name, gene_name) = (gene_ID[0:gene_ID.index(".")], gene_ID[gene_ID.index(".")+1:])
                                    #gene_beg       = CDS_domain_list[KBASE_DOMAINHIT_GENE_BEG_I]
                                    #gene_end       = CDS_domain_list[KBASE_DOMAINHIT_GENE_END_I]
                                    #gene_strand    = CDS_domain_list[KBASE_DOMAINHIT_GENE_STRAND_I]
                                    gene_hits_dict = CDS_domain_list[KBASE_DOMAINHIT_GENE_HITS_DICT_I]
                                    gene_hits_list = []
                                    for domfam in gene_hits_dict.keys():
                                        # skip CD hits for now
                                        if domfam[0:2] != 'PF' and domfam[0:3] != 'COG' and domfam[0:4] != 'TIGR':
                                            continue
                                        #Global_Domains[i][gene_name] = gene_hits_dict
                                        for hit in gene_hits_dict[domfam]:
                                            list_format_hit = hit
                                            list_format_hit.append (domfam)
                                            #list_format_hit[DOMHIT_BEG_I]      = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_BEG_J]
                                            #list_format_hit[DOMHIT_END_I]      = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_END_J]
                                            #list_format_hit[DOMHIT_EVALUE_I]   = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_EVALUE_J]
                                            #list_format_hit[DOMHIT_BITSCORE_I] = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_BITSCORE_J]
                                            #list_format_hit[DOMHIT_ALNPERC_I]  = hit[KBASE_DOMAINHIT_GENE_HITS_DICT_ALNPERC_J]
                                            gene_hits_list.append(list_format_hit)
                                            #print ("%s\t%s\t%s\t%s\t%d\t%d\t%f"%(genome_id, scaffold_id, gene_name, domfam, hit_beg, hit_end, hit_evalue))
                                    Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[i]]][contig_id_iter][gene_name] = sorted (gene_hits_list, key=sort_by_bitscore_key, reverse=True)


        # Workhorse for feature retrieval (KBase flavor)
        #   does more than window slice, so perhaps should split up functionality
        #
        def getFeatureSlicesKBase (ContigSet_names, \
                              PivotFeatures_IDs, \
                              genomebrowser_mode="contigs", \
                              genome_data_format="KBase", \
                              window_size=10000, \
                              genomebrowser_xshift=0):

            Feature_slices = []

            if genomebrowser_mode != "genome":

                if KBase_backend:  

                    for i,contig_name in enumerate(ContigSet_names):
                        if i >= max_rows:
                            break

                        (genome_id,scaffold_id) = contig_name.split(genome_contig_id_delim)
                        #print ("GENOME_ID: '"+genome_id+"' SCAFFOLD_ID: '"+scaffold_id+"' CONTIG_NAME: '"+contig_name+"'")  # DEBUG

                        genome_obj = Global_KBase_Genomes[genome_id]
                        #ass = Global_KBase_Assemblies[genome_id]

                        Feature_slice = []
                        Features_seen = set()
                        source = ""
                        pivot_pos = 0.0

                        try:
                            this_search_done = search_done[i]
                        except:
                            search_done.append([])                
                            search_done[i] = False
                            search_results.append([])
                            search_results[i] = []
                            for j,term in enumerate(Search_Terms):
                                search_results[i].append([])

                        # Get genome length
                        #if Global_State['PrimaryContig_len'] == 0 and i == 0:
                        contig_lens = gaAPI_get_contig_lengths(genome_obj=genome_obj)
                        if i == 0:
                            contig_GCs = gaAPI_get_contig_gc_content(genome_obj=genome_obj, contig_id_list=[scaffold_id])
                            Global_State['PrimaryContig_len'] = contig_lens[scaffold_id]                    
                            Global_State['PrimaryContig_GCavg'] = contig_GCs[scaffold_id]
                            Global_State['Contig_lens'] = []
                            Global_State['pivot_pos_list'] = []
                            #print ("%d"%Global_State['PrimaryContig_len'])
                        Global_State['Contig_lens'].append(contig_lens[scaffold_id])


                        # Find pivot feature and put in first position
                        #
                        contig_mode_xshift = 0
                        sci_name = gaAPI_get_scientific_name(genome_obj=genome_obj)
                        if sci_name == '' or sci_name.lower() == 'unknown':
                            source = genome_obj['info'][NAME_I]
                        else:
                            source = sci_name
                        try:
                            this_pivotfeature_id = PivotFeatures_IDs[i]
                        except:
                            PivotFeatures_IDs.append('')
                            Global_State['PivotFeatures_IDs'].append('')

                        if Global_State['genomebrowser_mode'] == 'contigs' \
                            or PivotFeatures_IDs[i] == '':

                            slice_beg = 1
                            slice_end = 100000
                            features = []
                            feature_slice_ids_list = []

        # RESTORE
        #                    feature_slice_ids = ga.get_feature_ids(group_by='region', filters={ "region_list": [{'contig_id':scaffold_id, 'strand':'?', 'start':slice_beg, 'length':slice_end-slice_beg+1}]})
                            #feature_slice_ids = ga.get_feature_ids(group_by='region')
                            feature_slice_ids = gaAPI_get_feature_ids(genome_obj=genome_obj)


                            # "by_region": dict<str contig_id, dict<str strand, dict<string range, list<string feature_id>>>>
                            for ctg_id in feature_slice_ids.keys():
                                if ctg_id != scaffold_id:  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                    continue
                                for strand in feature_slice_ids[ctg_id].keys():
                                    for f_range in feature_slice_ids[ctg_id][strand].keys():
                                        #print ("%s\t%s\t%s"%(ctg_id, strand, f_range))  # A
                                        [f_range_beg,f_range_end] = f_range.split('-')  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                        if int(f_range_beg) > slice_end or int(f_range_end) < slice_beg:
                                            continue
                                        feature_slice_ids_list.extend(feature_slice_ids[ctg_id][strand][f_range])
                            #features = ga.get_features(feature_id_list=feature_slice_ids_list)
                            features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=feature_slice_ids_list)

                            pivot_feature_rec = None
                            lowest_beg = 10000000000

                            for fid in features.keys():
                                #print (features[fid]['location'][0])  # DEBUG
                                strand = features[fid]['location'][0][KB_LOC_STRAND_I]
                                f_len = features[fid]['location'][0][KB_LOC_LEN_I]
                                if strand == '+':
                                    beg = features[fid]['location'][0][KB_LOC_BEG_I]
                                    end = beg + f_len - 1
                                else:
                                    end = features[fid]['location'][0][KB_LOC_BEG_I]
                                    beg = end - f_len + 1
                                #print ("%s\t%s\t%s\t%s\t%s\t%s"%(contig_id, ctg_id, fid, beg, end, strand))

                                if beg < lowest_beg:
                                    lowest_beg = beg
                                    f_type=features[fid]['type']
                                    pivot_feature_rec = build_feature_rec_kbase(features[fid], f_type=f_type, source_species=source, contig_i=i)

                            Feature_slice.append(pivot_feature_rec)
                            pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                            Global_State['pivot_pos_list'].append(pivot_pos)
                            if i == 0:
                                Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos

                            if Global_State['genomebrowser_mode'] == 'contigs':                    
                                contig_mode_xshift = 0.5*window_size - 0.5*(pivot_feature_rec['end_pos']-pivot_feature_rec['beg_pos'])

                        else:  # genomebrowser_mode != 'contigs' and != 'genome', so use PivotFeatures_IDs
                            #print ("FF "+str(i))  # DEBUG
                            fid = PivotFeatures_IDs[i]
                            #print ("GG "+str(fid))  # DEBUG
                            #features = ga.get_features(feature_id_list=[fid])
                            features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=[fid])
                            #print ("HH "+str(fid))  # DEBUG
                            f = features[fid]
                            #print ("II "+str(fid))  # DEBUG
                            f_type = f['type']
                            #print ("JJ "+str(fid))  # DEBUG

                            pivot_feature_rec = build_feature_rec_kbase(f, f_type=f_type, source_species=source, contig_i=i)
                            Feature_slice.append(pivot_feature_rec)
                            pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                            Global_State['pivot_pos_list'].append(pivot_pos)
                            if i == 0:
                                Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos


                        # Add in additional features within window.  Note: we want the duplicate pivot feature
                        #
                        slice_beg = pivot_pos + genomebrowser_xshift + contig_mode_xshift - 0.5*window_size
                        slice_end = pivot_pos + genomebrowser_xshift + contig_mode_xshift + 0.5*window_size
                        if slice_end < 1:
                            Feature_slices.append([Feature_slice[0]])
                            continue                 
                        if slice_beg < 1:
                            slice_beg = 1
                        features = []
                        feature_slice_ids_list = []
        # RESTORE
        #                feature_slice_ids = ga.get_feature_ids(group_by='region', filters={ "region_list": [{'contig_id':scaffold_id, 'strand':'?', 'start':slice_beg, 'length':slice_end-slice_beg+1}]})
                        #feature_slice_ids = ga.get_feature_ids(group_by='region')
                        feature_slice_ids = gaAPI_get_feature_ids(genome_obj=genome_obj)

                        #"by_region": dict<str contig_id, dict<str strand, dict<string range, list<string feature_id>>>>
                        for ctg_id in feature_slice_ids.keys():
                            if ctg_id != scaffold_id:  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                continue
                            for strand in feature_slice_ids[ctg_id].keys():
                                for f_range in feature_slice_ids[ctg_id][strand].keys():                    
                                    #print ("%s\t%s\t%s"%(ctg_id, strand, f_range)) # B
                                    [f_range_beg,f_range_end] = f_range.split('-')  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                    if int(f_range_beg) > slice_end or int(f_range_end) < slice_beg:
                                        continue
                                    feature_slice_ids_list.extend(feature_slice_ids[ctg_id][strand][f_range])
                        if len(feature_slice_ids_list) == 0:
                            Feature_slices.append([Feature_slice[0]])
                            continue  
                        #features = ga.get_features(feature_id_list=feature_slice_ids_list)
                        features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=feature_slice_ids_list)

                        for fid in features.keys():
                            f_type = features[fid]['type']
                            if features[fid].get('dna_sequence'):
                                dna_seq = features[fid]['dna_sequence']
                            else:
                                dna_seq = 'N'  # DEBUG
                            feature_rec = build_feature_rec_kbase(features[fid], f_type=f_type, source_species=source, contig_i=i, dna_seq=dna_seq)
                            Feature_slice.append(feature_rec)                        


                        # Check features against Search Terms
                        #
                        slice_beg = 1
                        slice_end = 10000000000
                        features = []
                        feature_slice_ids_list = []
        # RESTORE
        #                feature_slice_ids = ga.get_feature_ids(group_by='region', filters={ "region_list": [{'contig_id':scaffold_id, 'strand':'?', 'start':slice_beg, 'length':slice_end-slice_beg+1}]})
                        #feature_slice_ids = ga.get_feature_ids(group_by='region')
                        feature_slice_ids = gaAPI_get_feature_ids(genome_obj=genome_obj)

                        for ctg_id in feature_slice_ids.keys():
                            if ctg_id != scaffold_id:  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                continue
                            for strand in feature_slice_ids[ctg_id].keys():
                                for f_range in feature_slice_ids[ctg_id][strand].keys():                    
                                    #print ("%s\t%s\t%s"%(ctg_id, strand, f_range)) # C
                                    [f_range_beg,f_range_end] = f_range.split('-')  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                    if int(f_range_beg) > slice_end or int(f_range_end) < slice_beg:
                                        continue
                                    feature_slice_ids_list.extend(feature_slice_ids[ctg_id][strand][f_range])
                        if len(feature_slice_ids_list) == 0:
                            Feature_slices.append([Feature_slice[0]])
                            continue  
                        #features = ga.get_features(feature_id_list=feature_slice_ids_list)    
                        features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=feature_slice_ids_list)    

                        for fid in features.keys():
                            strand = features[fid]['location'][0][KB_LOC_STRAND_I]
                            f_len = features[fid]['location'][0][KB_LOC_LEN_I]
                            if strand == '+':
                                beg = features[fid]['location'][0][KB_LOC_BEG_I]
                                end = beg + f_len - 1
                            else:
                                end = features[fid]['location'][0][KB_LOC_BEG_I]
                                beg = end - f_len + 1                
                            pos_key = "%s,%d,%d"%(strand,beg,end)
                            if pos_key in Features_seen:
                                continue
                            else:
                                Features_seen.add(pos_key)                    

                            # check if a search hit
                            if not search_done[i]:
                                for j,match_flag in enumerate(search_term_match(feature_rec)):
                                    if match_flag:
                                        try:
                                            term_hit_in_genome = search_results[i][j][0]
                                        except:
                                            search_results[i][j] = []
                                        search_results[i][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                     'end_pos': feature_rec['end_pos'],
                                                                     'name': feature_rec['name'],
                                                                     'annot': feature_rec['annot']
                                                                    })

                        # sort results by position, store contig slice, and mark search as done on this contig
                        #
                        Sorted_Feature_slice = [Feature_slice[0]]  # retain pivot
                        for f in sorted (Feature_slice[1:], key=sort_by_beg_pos_key):
                            Sorted_Feature_slice.append (f)

                        Feature_slices.append(Sorted_Feature_slice)
                        search_done[i] = True   

                else:
                    print ("must use KBase_backend")

            # genomebrowser_mode == "genome"
            else:  
                if KBase_backend:  

                    for i in range (0,Global_State['genome_mode_n_rows']):
                        contig_name = ContigSet_names[0]
                        (genome_id,scaffold_id) = contig_name.split(genome_contig_id_delim)

                        genome_obj = Global_KBase_Genomes[genome_id]
                        #ass = Global_KBase_Assemblies[genome_id]

                        Feature_slice = []
                        Features_seen = set()
                        source = ""
                        pivot_pos = 0.0
                        track_xshift = window_size * (i - (Global_State['genome_mode_n_rows']-1)/2)  # e.g. n_rows=7 -> -3,-2,-1,0,1,2,3                

                        try:
                            this_search_done = search_done[0]
                        except:
                            search_done.append([])                
                            search_done[0] = False
                            search_results.append([])
                            search_results[0] = []
                            for j,term in enumerate(Search_Terms):
                                search_results[0].append([])

                        # Get genome length
                        #if Global_State['PrimaryContig_len'] == 0 and i == 0:
                        contig_lens = gaAPI_get_contig_lengths(genome_obj=genome_obj)
                        if Global_State['PrimaryContig_len'] == 0 and i == 0:
                            contig_GCs  = gaAPI_get_contig_gc_content(genome_obj=genome_obj, contig_id_list=[scaffold_id])
                            Global_State['PrimaryContig_len'] = contig_lens[scaffold_id]                    
                            Global_State['PrimaryContig_GCavg'] = contig_GCs[scaffold_id]
                            #print ("%d"%Global_State['PrimaryContig_len'])
                            Global_State['Contig_lens'].append(contig_lens[scaffold_id])


                        # Find pivot feature and put in first position
                        #
                        if i == 0:
                            try:
                                this_pivotfeature_id = PivotFeatures_IDs[0]
                            except:
                                PivotFeatures_IDs.append('')
                                Global_State['PivotFeatures_IDs'].append('')

                        if i > 0:
                            pivot_feature_rec = Feature_slices[0][0]  # use homolog feature from first contig only
                            pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                            Feature_slice.append(pivot_feature_rec)
                        else:
                            source = 'E. missingii'  # FIX THIS

                            if PivotFeatures_IDs[0] == '':
                                slice_beg = 1
                                slice_end = 100000
                                features = []
                                feature_slice_ids_list = []
        # RESTORE
        #                        feature_slice_ids = ga.get_feature_ids(group_by='region', filters={ "region_list": [{'contig_id':scaffold_id, 'strand':'?', 'start':slice_beg, 'length':slice_end-slice_beg+1}]})
                                #feature_slice_ids = ga.get_feature_ids(group_by='region')
                                feature_slice_ids = gaAPI_get_feature_ids(genome_obj=genome_obj)

                                #"by_region": dict<str contig_id, dict<str strand, dict<string range, list<string feature_id>>>>
                                for ctg_id in feature_slice_ids.keys():
                                    if ctg_id != scaffold_id:  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                        continue
                                    for strand in feature_slice_ids[ctg_id].keys():
                                        for f_range in feature_slice_ids[ctg_id][strand].keys():                    
                                            #print ("%s\t%s\t%s"%(ctg_id, strand, f_range)) # D
                                            [f_range_beg,f_range_end] = f_range.split('-')  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                            if int(f_range_beg) > slice_end or int(f_range_end) < slice_beg:
                                                continue
                                            feature_slice_ids_list.extend(feature_slice_ids[ctg_id][strand][f_range])
                                #features = ga.get_features(feature_id_list=feature_slice_ids_list)
                                features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=feature_slice_ids_list)

                                pivot_feature_rec = None
                                lowest_beg = 10000000000

                                for fid in features.keys():
                                    strand = features[fid]['location'][0][KB_LOC_STRAND_I]
                                    f_len = features[fid]['location'][0][KB_LOC_LEN_I]
                                    if strand == '+':
                                        beg = features[fid]['location'][0][KB_LOC_BEG_I]
                                        end = beg + f_len - 1
                                    else:
                                        end = features[fid]['location'][0][KB_LOC_BEG_I]
                                        beg = end - f_len + 1
                                    #print ("%s\t%s\t%s\t%s\t%s\t%s"%(contig_id, ctg_id, fid, beg, end, strand))

                                    if beg < lowest_beg:
                                        lowest_beg = beg
                                        f_type=features[fid]['type']
                                        pivot_feature_rec = build_feature_rec_kbase(features[fid], f_type=f_type, source_species=source, contig_i=i)

                                Feature_slice.append(pivot_feature_rec)
                                pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                                Global_State['pivot_pos_list'].append(pivot_pos)
                                if i == 0:
                                    Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos

                            else:  # PivotFeatures_IDs available
                                fid = PivotFeatures_IDs[0]
                                #features = ga.get_features(feature_id_list=[fid])
                                features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=[fid])
                                f = features[fid]
                                f_type = f['type']

                                pivot_feature_rec = build_feature_rec_kbase(f, f_type=f_type, source_species=source, contig_i=0)
                                Feature_slice.append(pivot_feature_rec)
                                pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                                Global_State['pivot_pos_list'].append(pivot_pos)
                                if i == 0:
                                    Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos


                        # Add in additional features within window.  Note: we want the duplicate pivot feature
                        #
                        slice_beg = int(round(pivot_pos + genomebrowser_xshift + track_xshift - 0.5*window_size))
                        slice_end = int(round(pivot_pos + genomebrowser_xshift + track_xshift + 0.5*window_size))
                        # DEBUG: temporary fix, but should make track_xshift correct for empty PivotFeatures_IDs
                        if slice_end < 1:
                            Feature_slices.append([Feature_slice[0]])
                            continue  
                        elif slice_beg < 1:
                            slice_beg = 1
                        #print ("%s\t%s\t%s\t%s"%(pivot_pos, genomebrowser_xshift, track_xshift, 0.5*window_size))
                        features = []
                        feature_slice_ids_list = []
        # RESTORE
        #                feature_slice_ids = ga.get_feature_ids(group_by='region', filters={ "region_list": [{'contig_id':scaffold_id, 'strand':'?', 'start':slice_beg, 'length':slice_end-slice_beg+1}]})
                        #feature_slice_ids = ga.get_feature_ids(group_by='region')
                        feature_slice_ids = gaAPI_get_feature_ids(genome_obj=genome_obj)

                        #"by_region": dict<str contig_id, dict<str strand, dict<string range, list<string feature_id>>>>
                        for ctg_id in feature_slice_ids.keys():
                            if ctg_id != scaffold_id:  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                continue
                            for strand in feature_slice_ids[ctg_id].keys():
                                for f_range in feature_slice_ids[ctg_id][strand].keys():                    
                                    #print ("%s\t%s\t%s"%(ctg_id, strand, f_range)) # E
                                    [f_range_beg,f_range_end] = f_range.split('-')  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                    if int(f_range_beg) > slice_end or int(f_range_end) < slice_beg:
                                        continue
                                    feature_slice_ids_list.extend(feature_slice_ids[ctg_id][strand][f_range])
                        if len(feature_slice_ids_list) == 0:
                            Feature_slices.append([Feature_slice[0]])
                            continue                    
                        #features = ga.get_features(feature_id_list=feature_slice_ids_list)
                        features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=feature_slice_ids_list)

                        for fid in features.keys():
                            f_type = features[fid]['type']
                            if features[fid].get('dna_sequence'):
                                dna_seq = features[fid]['dna_sequence']
                            else:
                                dna_seq = 'N'  # DEBUG
                            feature_rec = build_feature_rec_kbase(features[fid], f_type=f_type, source_species=source, contig_i=i, dna_seq=dna_seq)
                            Feature_slice.append(feature_rec)                        


                        # Check features against Search Terms
                        #
                        slice_beg = 1
                        slice_end = 100000000000
                        features = []
                        feature_slice_ids_list = []
        # RESTORE
        #                feature_slice_ids = ga.get_feature_ids(group_by='region', filters={ "region_list": [{'contig_id':scaffold_id, 'strand':'?', 'start':slice_beg, 'length':slice_end-slice_beg+1}]})
                        #feature_slice_ids = ga.get_feature_ids(group_by='region')
                        feature_slice_ids = gaAPI_get_feature_ids(genome_obj=genome_obj)

                        for ctg_id in feature_slice_ids.keys():
                            if ctg_id != scaffold_id:  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                continue
                            for strand in feature_slice_ids[ctg_id].keys():
                                for f_range in feature_slice_ids[ctg_id][strand].keys():                    
                                    #print ("%s\t%s\t%s"%(ctg_id, strand, f_range)) # F
                                    [f_range_beg,f_range_end] = f_range.split('-')  # SHOULDN'T BE NECESSARY IF get_feature_ids() WORKING
                                    if int(f_range_beg) > slice_end or int(f_range_end) < slice_beg:
                                        continue
                                    feature_slice_ids_list.extend(feature_slice_ids[ctg_id][strand][f_range])
                        if len(feature_slice_ids_list) == 0:
                            Feature_slices.append([Feature_slice[0]])
                            continue
                        #features = ga.get_features(feature_id_list=feature_slice_ids_list)    
                        features = gaAPI_get_features(genome_obj=genome_obj, feature_id_list=feature_slice_ids_list)    

                        for fid in features.keys():
                            strand = features[fid]['location'][0][KB_LOC_STRAND_I]
                            f_len = features[fid]['location'][0][KB_LOC_LEN_I]
                            if strand == '+':
                                beg = features[fid]['location'][0][KB_LOC_BEG_I]
                                end = beg + f_len - 1
                            else:
                                end = features[fid]['location'][0][KB_LOC_BEG_I]
                                beg = end - f_len + 1                
                            pos_key = "%s,%d,%d"%(strand,beg,end)
                            if pos_key in Features_seen:
                                continue
                            else:
                                Features_seen.add(pos_key)                    

                            # check if a search hit
                            if not search_done[0]:
                                for j,match_flag in enumerate(search_term_match(feature_rec)):
                                    if match_flag:
                                        try:
                                            term_hit_in_genome = search_results[0][j][0]
                                        except:
                                            search_results[0][j] = []
                                        search_results[0][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                     'end_pos': feature_rec['end_pos'],
                                                                     'name': feature_rec['name'],
                                                                     'annot': feature_rec['annot']
                                                                    })

                        # sort results by position, store contig slice, and mark search as done on this contig
                        #
                        Sorted_Feature_slice = [Feature_slice[0]]  # retain pivot
                        for f in sorted (Feature_slice[1:], key=sort_by_beg_pos_key):
                            Sorted_Feature_slice.append (f)

                        Feature_slices.append(Sorted_Feature_slice)
                        search_done[0] = True   

                else:
                    print ("must use KBase_backend")


            return Feature_slices


        # Workhorse for feature retrieval (Genbank flavor)
        #   does more than window slice, so perhaps should split up functionality
        #
        def getFeatureSlicesGenbank (ContigSet_names, \
                              PivotFeatures_IDs, \
                              genomebrowser_mode="contigs", \
                              genome_data_format="Genbank", \
                              window_size=10000, \
                              genomebrowser_xshift=0):

            Feature_slices = []

            if genomebrowser_mode != "genome":

                if genome_data_format == "Genbank":
                    for i,contig_name in enumerate(ContigSet_names):
                        if i >= max_rows:
                            break

                        try:
                            t = Global_Genbank_Genomes[i]
                        except:
                            (genome_id,scaffold_id) = contig_name.split(genome_contig_id_delim)
                            #print ("reading " + contig_name + " ...")
                            #genome_data_path = 'data/'+contig_name+'.'+gbk_ext
                            genome_data_path = genome_data_base_path+'/'+genome_id+genome_data_extra_subpath+'/'+scaffold_id+'.'+gbk_ext
                            print ("%d "%i+'reading '+genome_data_path)
                            Global_Genbank_Genomes.append (SeqIO.read(genome_data_path, 'genbank'))
                            #print ("%d"%len(Global_Genbank_Genomes))
                            #contig_seq = Global_Genbank_Genomes[i].seq
                            #print (contig_seq[-100:])
                            #print ("%d"%len(contig_seq))

        #                try:
        #                    feature_cache = Global_Features[i]
        #                except:
        #                    Global_Features[i] = {}

                        Feature_slice = []
                        Features_seen = set()
                        source = ""
                        pivot_pos = 0.0

                        try:
                            this_search_done = search_done[i]
                        except:
                            search_done.append([])                
                            search_done[i] = False
                            search_results.append([])
                            search_results[i] = []
                            for j,term in enumerate(Search_Terms):
                                search_results[i].append([])

                        # Get genome length
                        #if Global_State['PrimaryContig_len'] == 0 and i == 0:
                        if i == 0:
                            Global_State['PrimaryContig_len'] = len(Global_Genbank_Genomes[i].seq)
                            Global_State['PrimaryContig_GCavg'] = compute_GC (Global_Genbank_Genomes[i].seq)
                            Global_State['Contig_lens'] = []
                            Global_State['pivot_pos_list'] = []
                            #print ("%d"%Global_State['PrimaryContig_len'])
                        Global_State['Contig_lens'].append(len(Global_Genbank_Genomes[i].seq))

                        # Find pivot feature and put in first position
                        contig_mode_xshift = 0
                        for f in Global_Genbank_Genomes[i].features:                        
                            if f.type == "source":
                                source = f.qualifiers['organism'][0]
                                break

                        try:
                            this_pivotfeature_id = PivotFeatures_IDs[i]
                        except:
                            PivotFeatures_IDs.append('')
                            Global_State['PivotFeatures_IDs'].append('')

                        if Global_State['genomebrowser_mode'] == 'contigs' \
                            or PivotFeatures_IDs[i] == '':

                            pivot_feature_rec = None
                            lowest_beg = 10000000000
                            for f in Global_Genbank_Genomes[i].features:                        
                                #if f.type == "CDS" and "locus_tag" in f.qualifiers:                    
        #                        if (f.type == "CDS" and "gene" in f.qualifiers):
                                if (genome_annotation_system == 'KBase' and f.type == "CDS" and "gene" in f.qualifiers) \
                                    or (genome_annotation_system != 'KBase' and f.type == "CDS" and "locus_tag" in f.qualifiers) \
                                    or f.type == 'rRNA' \
                                    or f.type == 'tRNA' \
                                    or (f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR") \
                                    or (f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR spacer") \
                                    or (f.type == 'gene' and "pseudo" in f.qualifiers):

                                    if f.location.start == None or f.location.start <= 0:
                                        raise ValueError ("Bad Feature.  gene:'"+str(f.qualifiers.gene)+"' locus_tag:'"+str(f.qualifiers.locus_tag)+"'")

                                    if f.location.start < lowest_beg:
                                        lowest_beg = f.location.start
                                        pivot_feature_rec = build_feature_rec_genbank(f, f_type=f.type, source_species=source, contig_i=i)

                            Feature_slice.append(pivot_feature_rec)
                            pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                            Global_State['pivot_pos_list'].append(pivot_pos)
                            if i == 0:
                                Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos

                            if Global_State['genomebrowser_mode'] == 'contigs':                    
                                contig_mode_xshift = 0.5*window_size - 0.5*(pivot_feature_rec['end_pos']-pivot_feature_rec['beg_pos'])

                        else:             
                            for f in Global_Genbank_Genomes[i].features:        

                                #if f.type == "CDS" and "locus_tag" in f.qualifiers:                    
                                if (genome_annotation_system == 'KBase' and f.type == "CDS" and "gene" in f.qualifiers) \
                                    or (genome_annotation_system != 'KBase' and f.type == "CDS" and "locus_tag" in f.qualifiers):  # should it permit non-CDS anchor?             

                                    #elif f.qualifiers['locus_tag'][0] == PivotFeatures_IDs[i]:
                                    if (genome_annotation_system == 'KBase' and f.qualifiers['gene'][0] == PivotFeatures_IDs[i]) \
                                        or (genome_annotation_system != 'KBase' and f.qualifiers['locus_tag'][0] == PivotFeatures_IDs[i]):  
                                        pivot_feature_rec = build_feature_rec_genbank(f, f_type='CDS', source_species=source, contig_i=i)
                                        Feature_slice.append(pivot_feature_rec)
                                        pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                                        Global_State['pivot_pos_list'].append(pivot_pos)
                                        if i == 0:
                                            Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos
                                        break

                        # Add in additional features within window.  Note: we want the duplicate pivot feature
                        for f in Global_Genbank_Genomes[i].features:

                            if f.type == "CDS":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                # fast enough for now to build recs for all features, but might want to only do for window
                                feature_rec = build_feature_rec_genbank(f, f_type='CDS', source_species=source, contig_i=i)

                                # check if a search hit
                                if not search_done[i]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[i][j][0]
                                            except:
                                                search_results[i][j] = []
                                            search_results[i][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })

                                # add feature to slice if in viewable window
                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + contig_mode_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + contig_mode_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == "rRNA":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                feature_rec = build_feature_rec_genbank(f, f_type='rRNA', source_species=source, contig_i=i)

                                # check if a search hit
                                if not search_done[i]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[i][j][0]
                                            except:
                                                search_results[i][j] = []
                                            search_results[i][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + contig_mode_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + contig_mode_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == "tRNA":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)                        
                                feature_rec = build_feature_rec_genbank(f, f_type='tRNA', source_species=source, contig_i=i)

                                # check if a search hit
                                if not search_done[i]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[i][j][0]
                                            except:
                                                search_results[i][j] = []
                                            search_results[i][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + contig_mode_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + contig_mode_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)                        
                                fwd_dna_seq = Global_Genbank_Genomes[i].seq[f.location.start:f.location.end-1]  # BioPython shifts start by -1 but not end?
                                feature_rec = build_feature_rec_genbank(f, f_type='CRISPR', source_species=source, contig_i=i, dna_seq=fwd_dna_seq)

                                # check if a search hit
                                if not search_done[i]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[i][j][0]
                                            except:
                                                search_results[i][j] = []
                                            search_results[i][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + contig_mode_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + contig_mode_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR spacer":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                fwd_dna_seq = Global_Genbank_Genomes[i].seq[f.location.start:f.location.end-1]  # BioPython shifts start by -1 but not end?
                                feature_rec = build_feature_rec_genbank(f, f_type='CRISPR spacer', source_species=source, contig_i=i, dna_seq=fwd_dna_seq)

                                # check if a search hit
                                if not search_done[i]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[i][j][0]
                                            except:
                                                search_results[i][j] = []
                                            search_results[i][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + contig_mode_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + contig_mode_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == 'gene' and "pseudo" in f.qualifiers:
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)                        
                                feature_rec = build_feature_rec_genbank(f, f_type='pseudogene', source_species=source, contig_i=i)

                                # check if a search hit
                                if not search_done[i]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[i][j][0]
                                            except:
                                                search_results[i][j] = []
                                            search_results[i][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + contig_mode_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + contig_mode_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                        # sort results by position, store contig slice, and mark search as done on this contig
                        #
                        Sorted_Feature_slice = [Feature_slice[0]]  # retain pivot
                        for f in sorted (Feature_slice[1:], key=sort_by_beg_pos_key):
                            Sorted_Feature_slice.append (f)

                        Feature_slices.append(Sorted_Feature_slice)
                        search_done[i] = True   

                else:
                    print ("unknown data format for Genomes")

            # genomebrowser_mode == "genome"
            else:  

                if genome_data_format == "Genbank":
                    #for i in range (0,total_rows):
                    for i in range (0,Global_State['genome_mode_n_rows']):
                        contig_name = ContigSet_names[0]
                        try:
                            t=Global_Genbank_Genomes[0]
                        except:
                            (genome_id,scaffold_id) = contig_name.split(genome_contig_id_delim)
                            #print ("reading " + contig_name + " ...")
                            #genome_data_path = 'data/'+contig_name+'.'+gbk_ext
                            genome_data_path = genome_data_base_path+'/'+genome_id+genome_data_extra_subpath+'/'+scaffold_id+'.'+gbk_ext
                            print ("%d "%i+'reading '+genome_data_path)
                            Global_Genbank_Genomes.append (SeqIO.read(genome_data_path, 'genbank'))

                        Feature_slice = []
                        Features_seen = set()
                        source = ""
                        pivot_pos = 0.0
                        #track_xshift = window_size * (i - (total_rows-1)/2)  # e.g. total_rows=7 -> -3,-2,-1,0,1,2,3
                        track_xshift = window_size * (i - (Global_State['genome_mode_n_rows']-1)/2)  # e.g. n_rows=7 -> -3,-2,-1,0,1,2,3
                        #print ("%d %f"%(i,track_xshift))

                        try:
                            this_search_done = search_done[0]
                        except:
                            search_done.append([])                
                            search_done[0] = False
                            search_results.append([])
                            search_results[0] = []
                            for j,term in enumerate(Search_Terms):
                                search_results[0].append([])                

                        # Get genome length
                        if Global_State['PrimaryContig_len'] == 0 and i == 0:
                            Global_State['PrimaryContig_len'] = len(Global_Genbank_Genomes[i].seq)
                            Global_State['PrimaryContig_GCavg'] = compute_GC (Global_Genbank_Genomes[i].seq)
                            #print ("%d"%Global_State['PrimaryContig_len'])
                            Global_State['Contig_lens'].append(len(Global_Genbank_Genomes[i].seq))

                        # Find pivot feature and put in first position
                        if i == 0:
                            try:
                                this_pivotfeature_id = PivotFeatures_IDs[0]
                            except:
                                PivotFeatures_IDs.append('')
                                Global_State['PivotFeatures_IDs'].append('')              

                        if i > 0:
                            pivot_feature_rec = Feature_slices[0][0]  # use homolog feature from first contig only
                            pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                            Feature_slice.append(pivot_feature_rec)
                        else:
                            for f in Global_Genbank_Genomes[0].features:  
                                if f.type == "source":
                                    source = f.qualifiers['organism'][0]
                                    break

                            if PivotFeatures_IDs[0] == '':      
                                pivot_feature_rec = None
                                lowest_beg = 10000000000
                                for f in Global_Genbank_Genomes[0].features:                        
                                    #if f.type == "CDS" and "locus_tag" in f.qualifiers:                    
            #                        if (f.type == "CDS" and "gene" in f.qualifiers):
                                    if (genome_annotation_system == 'KBase' and f.type == "CDS" and "gene" in f.qualifiers) \
                                        or (genome_annotation_system != 'KBase' and f.type == "CDS" and "locus_tag" in f.qualifiers) \
                                        or f.type == 'rRNA' \
                                        or f.type == 'tRNA' \
                                        or (f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR") \
                                        or (f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR spacer") \
                                        or (f.type == 'gene' and "pseudo" in f.qualifiers):

                                        if f.location.start < lowest_beg:
                                            lowest_beg = f.location.start
                                            pivot_feature_rec = build_feature_rec_genbank(f, f_type=f.type, source_species=source, contig_i=0)

                                Feature_slice.append(pivot_feature_rec)
                                pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                                Global_State['pivot_pos_list'].append(pivot_pos)
                                if i == 0:
                                    Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos  

                            else:
                                for f in Global_Genbank_Genomes[0].features:  # should it permit non-CDS anchor?    
                                    #if f.type == "CDS" and "locus_tag" in f.qualifiers and f.qualifiers['locus_tag'][0] == PivotFeatures_IDs[i]:
                                    if (genome_annotation_system == 'KBase' and f.type == "CDS" and "gene" in f.qualifiers) \
                                        or (genome_annotation_system != 'KBase' and f.type == "CDS" and "locus_tag" in f.qualifiers):  
                                        #if f.qualifiers['gene'][0] == PivotFeatures_IDs[i]:
                                        if (genome_annotation_system == 'KBase' and f.qualifiers['gene'][0] == PivotFeatures_IDs[0]) \
                                            or (genome_annotation_system != 'KBase' and f.qualifiers['locus_tag'][0] == PivotFeatures_IDs[0]):
                                            pivot_feature_rec = build_feature_rec_genbank(f, f_type='CDS', source_species=source, contig_i=0)
                                            Feature_slice.append(pivot_feature_rec)
                                            pivot_pos = 0.5 * (pivot_feature_rec['beg_pos']+pivot_feature_rec['end_pos'])
                                            if i == 0:
                                                Global_State['PrimaryAnchor_pivot_pos'] = pivot_pos
                                            break

                        # Add in additional features within window.  Note: we want the duplicate pivot feature
                        for f in Global_Genbank_Genomes[0].features:
                            if f.type == "CDS":  # we'll add other types later when we have icons for them
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                # fast enough for now to build recs for all features, but might want to only do for window
                                feature_rec = build_feature_rec_genbank(f, f_type='CDS', source_species=source, contig_i=0)

                                # check if a search hit
                                if not search_done[0]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[0][j][0]
                                            except:
                                                search_results[0][j] = []
                                            search_results[0][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })                        

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + track_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + track_xshift + 0.5*window_size:

                                    #print ("%d %d %d"%(i, feature_rec['beg_pos'], feature_rec['end_pos']))
                                    Feature_slice.append(feature_rec)

                            elif f.type == "rRNA":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                feature_rec = build_feature_rec_genbank(f, f_type='rRNA', source_species=source, contig_i=0)

                                # check if a search hit
                                if not search_done[0]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[0][j][0]
                                            except:
                                                search_results[0][j] = []
                                            search_results[0][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })                          

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + track_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + track_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == "tRNA":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                feature_rec = build_feature_rec_genbank(f, f_type='tRNA', source_species=source, contig_i=0)

                                # check if a search hit
                                if not search_done[0]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[0][j][0]
                                            except:
                                                search_results[0][j] = []
                                            search_results[0][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })                          

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + track_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + track_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)                    

                            elif f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                fwd_dna_seq = Global_Genbank_Genomes[0].seq[f.location.start:f.location.end-1]  # BioPython shifts start by -1 but not end?
                                feature_rec = build_feature_rec_genbank(f, f_type='CRISPR', source_species=source, contig_i=0, dna_seq=fwd_dna_seq)

                                # check if a search hit
                                if not search_done[0]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[0][j][0]
                                            except:
                                                search_results[0][j] = []
                                            search_results[0][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })                          

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + track_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + track_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == 'gene' and "comment" in f.qualifiers and f.qualifiers['comment'][0] == "CRISPR spacer":
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                fwd_dna_seq = Global_Genbank_Genomes[0].seq[f.location.start:f.location.end-1]  # BioPython shifts start by -1 but not end?
                                feature_rec = build_feature_rec_genbank(f, f_type='CRISPR spacer', source_species=source, contig_i=0, dna_seq=fwd_dna_seq)

                                # check if a search hit
                                if not search_done[0]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[0][j][0]
                                            except:
                                                search_results[0][j] = []
                                            search_results[0][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })                          

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + track_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + track_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)

                            elif f.type == 'gene' and "pseudo" in f.qualifiers:
                                pos_key = "%d,%d,%d"%(f.location.strand,f.location.start+1,f.location.end)
                                if pos_key in Features_seen:
                                    continue
                                else:
                                    Features_seen.add(pos_key)
                                feature_rec = build_feature_rec_genbank(f, f_type='pseudogene', source_species=source, contig_i=0)

                                # check if a search hit
                                if not search_done[0]:
                                    for j,match_flag in enumerate(search_term_match(feature_rec)):
                                        if match_flag:
                                            try:
                                                term_hit_in_genome = search_results[0][j][0]
                                            except:
                                                search_results[0][j] = []
                                            search_results[0][j].append({'beg_pos': feature_rec['beg_pos'],
                                                                         'end_pos': feature_rec['end_pos'],
                                                                         'name': feature_rec['name'],
                                                                         'annot': feature_rec['annot']
                                                                        })                          

                                if feature_rec['end_pos'] >= pivot_pos + genomebrowser_xshift + track_xshift - 0.5*window_size and \
                                    feature_rec['beg_pos'] <= pivot_pos + genomebrowser_xshift + track_xshift + 0.5*window_size:

                                    Feature_slice.append(feature_rec)


                        # sort results by position, store contig slice, and mark search as done on this contig
                        #
                        Sorted_Feature_slice = [Feature_slice[0]]  # retain pivot
                        for f in sorted (Feature_slice[1:], key=sort_by_beg_pos_key):
                            Sorted_Feature_slice.append (f)

                        Feature_slices.append(Sorted_Feature_slice)
                        search_done[0] = True

                else:
                    print ("unknown data format for Genomes")

            return Feature_slices


        # Slider events
        #
        class DraggableRectangleSlider:
            def __init__(self, rect, func='', move='x', min_pos=0.0, max_pos=1.0, step=0.1, tol=0.01):
                self.rect = rect
                self.press = None
                self.func = func
                self.move = move.lower()
                self.min_pos = min_pos
                self.max_pos = max_pos
                self.step = step
                self.tol = tol
                self.n_steps = int(round((max_pos-min_pos)/step))
                self.last_pos = 0
                if self.move != 'x' and self.move != 'y':
                    print ("bad DraggableRectSlider move: '"+ self.move +"' for func: '"+ self.func +"'")

            def connect(self):
                'connect to all the events we need'
                self.cidpress = self.rect.figure.canvas.mpl_connect('button_press_event', self.on_press)
                self.cidrelease = self.rect.figure.canvas.mpl_connect('button_release_event', self.on_release)
                self.cidmotion = self.rect.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

            def on_press(self, event):
                'on button press we will see if the mouse is over us and store some data'
                if event.inaxes != self.rect.axes: return

                contains, attrd = self.rect.contains(event)
                if not contains: return
                #print ('event contains', self.rect.xy)
                x0, y0 = self.rect.xy
                self.press = x0, y0, event.xdata, event.ydata
                if self.move == 'x':
                    self.last_pos = x0
                elif self.move == 'y':
                    self.last_pos = y0

            def on_motion(self, event):
                'on motion we will move the rect if the mouse is over us'
                update_draggable = False
                if self.press is None: return
                if event.inaxes != self.rect.axes: return
                x0, y0, xpress, ypress = self.press
                #last_pos = self.last_pos
                step = self.step
                tol = self.tol
                min_pos = self.min_pos
                max_pos = self.max_pos
                dx = event.xdata - xpress
                dy = event.ydata - ypress
                #print ('x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f, last_pos, step'%(x0, xpress, event.xdata, dx, x0+dx, self.last_pos, self.step))
                if self.move == 'x':
                    if x0+dx <= min_pos + tol:
                        if self.last_pos != min_pos:
                            self.last_pos = min_pos
                            nearest_step_i = 0
                            update_draggable = True
                    elif x0+dx >= max_pos - tol:
                        if self.last_pos != max_pos:
                            self.last_pos = max_pos
                            nearest_step_i = int(round((max_pos-min_pos)/step))
                            update_draggable = True
                    else:
                        nearest_step_i = int(round((x0+dx-min_pos)/step))
                        new_pos = min_pos + nearest_step_i*step
                        if new_pos != self.last_pos:
                            self.last_pos = new_pos
                            update_draggable = True

                    if update_draggable:
                        self.rect.set_x(self.last_pos)
                        #self.rect.set_x(x0+dx)
                        #self.rect.set_y(y0)
                elif self.move == 'y':
                    if y0+dy <= min_pos + tol:
                        if self.last_pos != min_pos:
                            self.last_pos = min_pos
                            nearest_step_i = 0
                            update_draggable = True
                    elif y0+dy >= max_pos - tol:
                        if self.last_pos != max_pos:
                            self.last_pos = max_pos
                            nearest_step_i = int(round((max_pos-min_pos)/step))
                            update_draggable = True
                    else:
                        nearest_step_i = int(round((y0+dy-min_pos)/step))
                        new_pos = min_pos + nearest_step_i*step
                        if new_pos != self.last_pos:
                            self.last_pos = new_pos
                            update_draggable = True

                    if update_draggable:
                        self.rect.set_y(self.last_pos)
                        #self.rect.set_x(x0)
                        #self.rect.set_y(y0+dy)

                if update_draggable:
                    self.rect.figure.canvas.draw()
                    if self.func == 'zoom':
                        Global_State['genomebrowser_zoom'] = int(round((self.last_pos-min_pos)/step))
                        Global_State['genomebrowser_window_bp_width'] = Global_State['def_genomebrowser_window_bp_width']*(2**Global_State['genomebrowser_zoom'])
                    elif self.func == 'pan':
                        Global_State['genomebrowser_xshift'] = Global_State['PrimaryContig_len']*(self.last_pos-min_pos)/(max_pos-min_pos) - Global_State['PrimaryAnchor_pivot_pos']
                    elif self.func == 'mode':
                        Global_State['genomebrowser_mode'] = mode_names[nearest_step_i]
                    elif self.func == 'color':
                        Global_State['genomebrowser_color_namespace'] = color_namespace_names[nearest_step_i]

                    update_genomebrowser_panel (ax_center)  # This should happen first so Features and Domains are updated
                    update_sidenav_panel (ax_left)
                    update_mode_panel (ax_top_left)

            def on_release(self, event):
                'on release we reset the press data'
                self.press = None
                self.rect.figure.canvas.draw()

            def disconnect(self):
                'disconnect all the stored connection ids'
                self.rect.figure.canvas.mpl_disconnect(self.cidpress)
                self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
                self.rect.figure.canvas.mpl_disconnect(self.cidmotion)


        # handle contig nav events
        #
        def onpick_contig_nav (event):
            artist = event.artist
            mouseevent = event.mouseevent
            ax = artist.axes
            if artist.get_label() != "circular_contig_nav" and artist.get_label() != "linear_contig_nav":
                return
            artist_info = artist.gid.split(",")

            repaint = False
            if artist.get_label() == "circular_contig_nav":
                genome_i = int(artist_info[0])
                x0 = float(artist_info[1])
                y0 = float(artist_info[2])
                y_scale = float(artist_info[3])
                dx = mouseevent.xdata - x0
                dy = (mouseevent.ydata - y0) / y_scale

                cosine = dx / sqrt(dx*dx+dy*dy)
                # Math.acos and Math.asin are in radians, and acos returns 0-pi
                theta = acos(cosine)
                if dy < 0:
                    theta = 2*pi - theta            
                if theta < pi/2.0:
                    phi = pi/2.0 - theta
                else:
                    phi = 2.0*pi - theta + pi/2.0
                genome_frac_pos = phi / (2.0*pi)

                new_genome_pos = int(round(genome_frac_pos * Global_State['PrimaryContig_len']))

                if (new_genome_pos <= 0):
                    new_genome_pos = 1
                elif new_genome_pos > Global_State['PrimaryContig_len']:
                    new_genome_pos = Global_State['PrimaryContig_len']

                Global_State['genomebrowser_xshift'] = new_genome_pos - Global_State['PrimaryAnchor_pivot_pos']
                repaint = True

            elif artist.get_label() == "linear_contig_nav":
                genome_i = int(artist_info[0])
                x0 = float(artist_info[1])
                w = float(artist_info[2])
                dx = mouseevent.xdata - x0
                #dy = mouseevent.ydata - y0

                new_genome_pos = int(round(Global_State['Contig_lens'][genome_i] * dx/w))

                if (new_genome_pos <= 0):
                    new_genome_pos = 1
                elif new_genome_pos > Global_State['Contig_lens'][genome_i]:
                    new_genome_pos = Global_State['Contig_lens'][genome_i]

                Global_State['genomebrowser_xshift'] = new_genome_pos - Global_State['pivot_pos_list'][genome_i]        
                repaint = True

            if repaint:
                #update_control_panel(ax_top_center)  # for some reason this causes things to drag.  memory leak for event listeners?
                # FIX THIS
                #pos_knob_patch.set_x(self.last_pos)  # global ref (total kludge).  not working yet.        
                update_genomebrowser_panel (ax_center)
                update_sidenav_panel (ax_left)


        # handle pan events
        #
        def onpick_pan (event):
            artist = event.artist
            #mouseevent = event.mouseevent
            #ax = artist.axes
            if artist.get_label() != "pan_button":
                return
            control = artist.gid

            repaint = False
            if control == "reset_pan":
                Global_State['genomebrowser_xshift'] = 0
                repaint = True
            elif control == "left_pan":
                if Global_State['genomebrowser_mode'] == "genome":
                    Global_State['genomebrowser_xshift'] -= 1.0*Global_State['genomebrowser_window_bp_width']
                else:
                    Global_State['genomebrowser_xshift'] -= 0.4*Global_State['genomebrowser_window_bp_width']
                repaint = True
            elif control == "right_pan":
                if Global_State['genomebrowser_mode'] == "genome":
                    Global_State['genomebrowser_xshift'] += 1.0*Global_State['genomebrowser_window_bp_width']
                else:
                    Global_State['genomebrowser_xshift'] += 0.4*Global_State['genomebrowser_window_bp_width']
                repaint = True

            if repaint:
                #update_control_panel(ax_top_center)  # for some reason this causes things to drag.  memory leak for event listeners?
                # FIX THIS
                #pos_knob_patch.set_x(self.last_pos)  # global ref (total kludge).  not working yet.        
                update_genomebrowser_panel (ax_center)
                update_sidenav_panel (ax_left)


        # Feature events
        #
        """
        # Something about the object solution here really makes ipython drag.
        #   Implemented solution with artist labels, feature_artist_label_to_popup_box and
        #   feature_artist_label_to_feature dicts and labels to link artists ugly but effective
        class InteractiveGenomeFeature:
            def __init__(self, arrow, feature):
                self.arrow = arrow
                self.feature = feature
                self.dialog_box = None
                self.box = None

            def connect(self):
                'connect to all the events we need'
                self.cidpick = self.arrow.figure.canvas.mpl_connect('pick_event', self.on_pick)

            def on_pick (self, event):
                artist = event.artist
                mouseevent = event.mouseevent
                arrow = artist
                ax = arrow.axes
                #arrow.set_color("black")
                self.dialog_box = Rectangle((mouseevent.xdata, mouseevent.ydata), 0.1, 0.1, \
                                            facecolor="aliceblue", edgecolor="gray", alpha=1.0, zorder=4)
                ax.add_patch(self.dialog_box)

                #arrow.remove()
                arrow.figure.canvas.draw()

            def disconnect(self):
                'disconnect all the stored connection ids'
                self.arrow.figure.canvas.mpl_disconnect(self.cidpick)

        """
        # handle feature events
        #
        def onpick_feature (event):
            artist = event.artist
            mouseevent = event.mouseevent
            ax = artist.axes
            if artist.get_label() != "feature" and artist.get_label() != "popup_destroy":
                return
            feature_index = artist.gid

            # popup config
            annot_word_limit = 5
            popup_x_margin = 0.015
            popup_y_margin = 0.020 * 7/(total_rows+1)
            close_button_w = .02
            close_button_h = .03 * 7/(total_rows+1)
            close_button_marker_y_shift = 0.0037 * 7/(total_rows+1)
            popup_text_w_scaling = 0.0085                     # appropriate value for fontsize 10
            popup_text_h_scaling = 0.0300 * 7/(total_rows+1)  # appropriate value for fontsize 10

            #popup_zorder += 2   # why doesn't this work even when popup_zorder is a global?
            Global_State['popup_zorder'] += 4
            popup_zorder = Global_State['popup_zorder']

            if not feature_index in feature_artist_label_to_popup_box or \
                not feature_artist_label_to_popup_box[feature_index]:

                feature = feature_artist_label_to_feature[feature_index]

                # build info strings
                word_join = " "
                row_join = "\n"

                # popup_info_1 (name)
                popup_info_str1_row_split = [feature['name']]

                # popup_info_2 (annot)
                annot_processed = False
                if not "\n" in feature['annot']:
                    popup_info_str2_word_split = feature['annot'].split(" ")
                    word_count = len(popup_info_str2_word_split)
                    if word_count <= annot_word_limit:
                        if feature['annot'] != '':
                            popup_info_str2_row_split = [feature['annot']]
                        annot_processed = True
                if not annot_processed:
                    popup_info_str2_row_split = []
                    for row in feature['annot'].split("\n"):
                        #if row == '':
                        #    continue
                        popup_info_str2_word_split = row.split(" ")
                        word_count = len(popup_info_str2_word_split)
                        if word_count <= annot_word_limit:
                            popup_info_str2_row_split.append(row)
                        else:
                            imax = int(round(word_count/annot_word_limit))
                            indent = ''
                            for i in range(0,imax):
                                row_slice = popup_info_str2_word_split[i*annot_word_limit:(i+1)*annot_word_limit]
                                if i > 0:
                                    indent = '  '
                                popup_info_str2_row_split.append(indent+word_join.join(row_slice))
                            remainder = len(popup_info_str2_word_split)%annot_word_limit
                            if remainder > 0:
                                row_slice = popup_info_str2_word_split[imax*annot_word_limit:]
                                indent = '  '
                                popup_info_str2_row_split.append(indent+word_join.join(row_slice))

                # popup_info_3 (rest of info)
                popup_info_str3_row_split = ["locus_tag: %s"%feature['locus_tag'], \
                                             "id(s): %s"%feature['ID'], \
                                             "E.C.: %s"%feature['EC_number'], \
                                             "type: %s"%feature['type'], \
                                             "position: %d-%d (%s)"%(feature['beg_pos'],feature['end_pos'],feature['strand']), \
                                             feature['source_species']
                                            ]

                # get dimensions and shift popup if hitting walls  
                #
                popup_info_row_cnt = len(popup_info_str1_row_split) + \
                                     len(popup_info_str2_row_split) + \
                                     1 + \
                                     len(popup_info_str3_row_split)
                popup_info_h = popup_text_h_scaling * popup_info_row_cnt

                longest_info_str = 0
                for popup_info_str_row_split in [popup_info_str1_row_split, popup_info_str2_row_split, popup_info_str3_row_split]:
                #for popup_info_str_row_split in [popup_info_str3_row_split]:
                    for popup_info_str in popup_info_str_row_split:
                        row_len = len(popup_info_str)
                        if row_len > longest_info_str:
                            longest_info_str = row_len
                popup_info_w = popup_text_w_scaling * longest_info_str

                popup_box_w = popup_info_w + 2*popup_x_margin + 1.5*close_button_w
                popup_box_h = popup_info_h + 2*popup_y_margin + 1.0*close_button_h
                popup_x_shift = 0
                popup_y_shift = 0
                if mouseevent.xdata+popup_box_w > 1.0:
                    popup_x_shift = 1.0 - (mouseevent.xdata+popup_box_w) - 0.01  # extra .01 to not hit wall
                if mouseevent.ydata+popup_box_h > 1.0:
                    popup_y_shift = 1.0 - (mouseevent.ydata+popup_box_h) - 0.01  # extra .01 to not hit wall


                # Add Box
                #
                popup_box = Rectangle((mouseevent.xdata+popup_x_shift, mouseevent.ydata+popup_y_shift), \
                                      popup_box_w, \
                                      popup_box_h, \
                                      facecolor=popup_bg_color, edgecolor="gray", \
                                      alpha=1.0, zorder=popup_zorder)
                popup_box.gid = feature_index  # if we want to make it and contents draggable later
                popup_box.set_path_effects([path_effects.PathPatchEffect(offset=(1, -1), facecolor='black', alpha=0.75),
                                            path_effects.PathPatchEffect(facecolor=popup_bg_color,
                                                                         edgecolor="gray")
                                           ])
                ax.add_patch(popup_box)

                # Add Info
                #
                text_y_shift = popup_text_h_scaling*len(popup_info_str3_row_split) + \
                               popup_text_h_scaling*1.0 + \
                               popup_text_h_scaling*len(popup_info_str2_row_split) + \
                               popup_text_h_scaling*0.5
                popup_info_1 = ax.text(mouseevent.xdata+popup_x_shift+popup_x_margin, \
                                       mouseevent.ydata+popup_y_shift+popup_y_margin + text_y_shift, \
                                        row_join.join(popup_info_str1_row_split), \
                                        verticalalignment="bottom", \
                                        horizontalalignment="left", \
                                        transform=ax.transAxes, \
                                        color="black", \
                                        fontsize=12, \
                                        zorder=popup_zorder+1)
                #popup_info_1.set_label(feature_index)
                #ax.add_patch(popup_info_1)   # text cannot be added as patch to axes, or else things get wonky

                text_y_shift = popup_text_h_scaling*len(popup_info_str3_row_split) + \
                               popup_text_h_scaling*1.0
                popup_info_2 = ax.text(mouseevent.xdata+popup_x_shift+popup_x_margin, \
                                       mouseevent.ydata+popup_y_shift+popup_y_margin + text_y_shift, \
                                        row_join.join(popup_info_str2_row_split), \
                                        verticalalignment="bottom", \
                                        horizontalalignment="left", \
                                        transform=ax.transAxes, \
                                        color="red", \
                                        fontsize=10, \
                                        zorder=popup_zorder+1)
                #popup_info_2.set_label(feature_index)

                text_y_shift = 0
                popup_info_3 = ax.text(mouseevent.xdata+popup_x_shift+popup_x_margin, \
                                       mouseevent.ydata+popup_y_shift+popup_y_margin + text_y_shift, \
                                       row_join.join(popup_info_str3_row_split), \
                                       verticalalignment="bottom", \
                                       horizontalalignment="left", \
                                       transform=ax.transAxes, \
                                       color="black", \
                                       fontsize=10, \
                                       zorder=popup_zorder+1)
                #popup_info_3.set_label(feature_index)


                # Add Close button
                #
                close_button = Rectangle((mouseevent.xdata+popup_x_shift+popup_box_w-popup_x_margin-close_button_w,
                                          mouseevent.ydata+popup_y_shift+popup_box_h-popup_y_margin-close_button_h), \
                                          close_button_w, \
                                          close_button_h, \
                                          facecolor="white", edgecolor="gray", \
                                          picker=True, \
                                          alpha=1.0, zorder=popup_zorder+2)
                close_button.gid = feature_index
                close_button.set_label("popup_destroy")
                ax.add_patch(close_button)
                #close_button.figure.canvas.mpl_connect ('pick_event', onpick_feature)  # already listening
                close_button_x = ax.text(mouseevent.xdata+popup_x_shift+popup_box_w-popup_x_margin-0.5*close_button_w,
                                         mouseevent.ydata+popup_y_shift+popup_box_h-popup_y_margin-0.5*close_button_h-close_button_marker_y_shift,
                                         'x', \
                                         verticalalignment="center", \
                                         horizontalalignment="center", \
                                         transform=ax.transAxes, \
                                         color="gray", \
                                         fontsize=10, \
                                         zorder=popup_zorder+3)
                #close_button_x.set_label(feature_index)

                # Register artists for moving or later destruction
                #
                feature_artist_label_to_popup_box[feature_index] = []
                feature_artist_label_to_popup_box[feature_index].append(popup_box)
                feature_artist_label_to_popup_box[feature_index].append(popup_info_1)
                feature_artist_label_to_popup_box[feature_index].append(popup_info_2)
                feature_artist_label_to_popup_box[feature_index].append(popup_info_3)
                feature_artist_label_to_popup_box[feature_index].append(close_button)
                feature_artist_label_to_popup_box[feature_index].append(close_button_x)

            # Close popup (conflicts with close button)   
            #
            else: 
                popup_box = feature_artist_label_to_popup_box[feature_index][0]
                popup_box.remove()
                popup_info_1 = feature_artist_label_to_popup_box[feature_index][1]
                popup_info_1.remove()
                popup_info_2 = feature_artist_label_to_popup_box[feature_index][2]
                popup_info_2.remove()
                popup_info_3 = feature_artist_label_to_popup_box[feature_index][3]
                popup_info_3.remove()
                close_button = feature_artist_label_to_popup_box[feature_index][4]
                close_button.remove()
                close_button_x = feature_artist_label_to_popup_box[feature_index][5]
                close_button_x.remove()

                feature_artist_label_to_popup_box[feature_index] = None
                feature_index  = None

            artist.figure.canvas.draw()


        # Diplay coordinate transform
        #
        def disp_coord_transform (x, pivot_pos, pivot_strand, genomebrowser_window_bp_width, genomebrowser_xshift, track_xshift, contig_mode_xshift):

            track_disp_scale_factor = (1.0-right_margin-left_margin) / genomebrowser_window_bp_width
            #transformed_x = track_disp_scale_factor * (x - (pivot_pos-0.5*genomebrowser_window_bp_width))

            if pivot_strand == '-' and Global_State['genomebrowser_mode'] == 'tree':
                transformed_x = track_disp_scale_factor * ((pivot_pos - x) - (genomebrowser_xshift+track_xshift+contig_mode_xshift-0.5*genomebrowser_window_bp_width))
            else:
                # ORIGINAL WORKING FORM
                #transformed_x = track_disp_scale_factor * (x - (pivot_pos+genomebrowser_xshift+track_xshift+contig_mode_xshift-0.5*genomebrowser_window_bp_width))
                transformed_x = track_disp_scale_factor * ((x - pivot_pos) - (genomebrowser_xshift+track_xshift+contig_mode_xshift - 0.5*genomebrowser_window_bp_width))
                
            return left_margin + transformed_x


        # draw feature element
        #
        def draw_feature_element (ax, \
                             row_n, \
                             max_row_n, \
                             feature_i, \
                             feature_j, \
                             feature, \
                             pivot_pos, \
                             pivot_strand, \
                             track_xshift, \
                             contig_mode_xshift, \
                             pivot_feature_flag, \
                             annot_repeat, \
                             seq_repeat, \
                             this_label_show_top, this_label_show_bot):

            left_cropped = False
            right_cropped = False
            feature_oi_flag = False
            pivot_feature_flag = False
            if feature['beg_pos'] <= pivot_pos and feature['end_pos'] >= pivot_pos:
                pivot_feature_flag = True

            gene_id = feature['ID']
            gene_name = feature['name']
            if feature.get('gene_name'):
                gene_name = feature['gene_name']
            
            # Determine if a search result feature
            try:
                if Global_State['genomebrowser_mode'] == "genome":
                    contig_search_results = search_results[0]
                else:
                    contig_search_results = search_results[feature_i]
            except:
                contig_search_results = []
            for j,result_list in enumerate(contig_search_results):
                for result in contig_search_results[j]:
                    if result['beg_pos'] == feature['beg_pos'] and result['end_pos'] == feature['end_pos']:
                        feature_oi_flag = True
                        highlight_color = search_color_names[j % len(search_color_names)]
                        break

            # window coords
            if Global_State['genomebrowser_mode'] == 'tree' and pivot_strand == '-':
                genome_beg_pos = feature['end_pos']
                genome_end_pos = feature['beg_pos']
            else:
                genome_beg_pos = feature['beg_pos']
                genome_end_pos = feature['end_pos']
                
            window_beg_pos = disp_coord_transform (feature['beg_pos'], \
                                                   pivot_pos, \
                                                   pivot_strand, \
                                                   Global_State['genomebrowser_window_bp_width'], \
                                                   Global_State['genomebrowser_xshift'], \
                                                   track_xshift, \
                                                   contig_mode_xshift)
            window_end_pos = disp_coord_transform (feature['end_pos'],
                                                   pivot_pos, \
                                                   pivot_strand, \
                                                   Global_State['genomebrowser_window_bp_width'], \
                                                   Global_State['genomebrowser_xshift'], \
                                                   track_xshift, \
                                                   contig_mode_xshift)

            # handle tiny number problems  (FIX: NOT QUITE WORKING)
            if window_end_pos < window_beg_pos:
                window_end_pos = window_beg_pos + .00000001
            if window_beg_pos > 1.0-right_margin:
                window_beg_pos = 1.0-right_margin
            if window_end_pos < left_margin:
                window_end_pos = left_margin

            # handle features that exceed margins
            if window_beg_pos < left_margin:
                window_beg_pos = left_margin
                left_cropped = True
            if window_end_pos > 1.0-right_margin:
                window_end_pos = 1.0-right_margin
                right_cropped = True     

            # deal with short features
            if feature['type'] != "CRISPR" and feature['type'] != "CRISPR spacer":
                if window_end_pos - window_beg_pos < head_l:
                    feature_element_base_length = .00000001
                    arrow_params = {'head_width':head_w, 'head_length':window_end_pos-window_beg_pos}
                else:
                    feature_element_base_length = window_end_pos - window_beg_pos - head_l
                    arrow_params = {'head_width':head_w, 'head_length':head_l}
            elif feature['type'] == "CRISPR":
                feature_element_base_length = .00000001
                arrow_params = {'head_width':0.75*head_w, 'head_length':0.5*(window_end_pos-window_beg_pos)}
            elif feature['type'] == "CRISPR spacer":
                feature_element_base_length = window_end_pos - window_beg_pos


            # not going to change window_beg_pos again, so can set x_pos and y_pos now
            x_pos = window_beg_pos   
            # FIX
            y_pos_base = bottom_margin + (1.0-bottom_margin-top_margin)*row_delta*(total_rows-row_n)
            #y_pos_base = bottom_margin + (1.0-bottom_margin-top_margin)*row_delta*(max_row_n-row_n)
            #y_pos_base = 1.0-0.02*row_n  # DEBUG

            
            # shift up to match tree due to scale bar at bottom
            # Strangely, this has no effect, so removing it
            #if Global_State['genomebrowser_mode'] == 'tree':
            #    y_pos_base += .5 / float(max_row_n)

            # direction
            if Global_State['genomebrowser_mode'] == 'tree':
                if pivot_strand == '-':
                    if feature['strand'] == '+':
                        direction = 'rev'
                        feature_element_start_pos = x_pos+feature_element_base_length+head_l
                        #label_start_pos = x_pos+head_l
                        label_start_pos = x_pos
                        dir_adj = -1.0
                        drop_shadow_xshift = 1
                    else:
                        direction = 'fwd'
                        feature_element_start_pos = x_pos
                        label_start_pos = x_pos
                        dir_adj = 1.0
                        drop_shadow_xshift = 0
                else:
                    if feature['strand'] == '-':
                        direction = 'rev'
                        feature_element_start_pos = x_pos+feature_element_base_length+head_l
                        #label_start_pos = x_pos+head_l
                        label_start_pos = x_pos
                        dir_adj = -1.0
                        drop_shadow_xshift = 1
                    else:
                        direction = 'fwd'
                        feature_element_start_pos = x_pos
                        label_start_pos = x_pos
                        dir_adj = 1.0
                        drop_shadow_xshift = 0
            else:
                if feature['strand'] == '-':
                    direction = 'rev'
                    feature_element_start_pos = x_pos+feature_element_base_length+head_l
                    #label_start_pos = x_pos+head_l
                    label_start_pos = x_pos
                    dir_adj = -1.0
                    drop_shadow_xshift = 1
                else:
                    direction = 'fwd'
                    feature_element_start_pos = x_pos
                    label_start_pos = x_pos
                    dir_adj = 1.0
                    drop_shadow_xshift = 0


            # color
            #
            feature_element_color = "lightgray"
            if feature['type'] == "pseudogene":
                feature_element_color = "gray"

            elif feature['type'] == "rRNA" or feature['type'] == "tRNA":
                feature_element_color = "black"

            elif feature['type'] == "CRISPR":
                feature_element_color = "black"        
                if feature['dna_seq'] in seq_repeat:
                    feature_element_color = color_names[sum([ord(c) for c in feature['dna_seq']]) % len(color_names)]           
            elif feature['type'] == "CRISPR spacer":
                feature_element_color = "lightgray"
                if feature['dna_seq'] in seq_repeat:
                    feature_element_color = color_names[sum([ord(c) for c in feature['dna_seq']]) % len(color_names)]

            elif feature['type'] == "CDS":

                if Global_State['genomebrowser_color_namespace'] == "annot":

                    if feature.get('gene_name'):
                        feature_element_color = color_names[sum([ord(c) for c in feature['gene_name']]) % len(color_names)]

                    elif 'annot' in feature:
                        if feature['annot'] != '' and \
                            (feature['annot'].lower() == "hypothetical protein" \
                            or feature['annot'].lower() == "conserved hypothetical protein" \
                            or feature['annot'].lower() == "orf, hypothetical protein" \
                            or feature['annot'].lower() == "orf, conserved hypothetical protein" \
                            or feature['annot'].lower() == "hypothetical protein; orphan" \
                            or feature['annot'].lower() == "hypothetical protein; genus orphan"):

                            feature_element_color = "darkgray"
                        elif feature['annot'] == '':
                            feature_element_color = "lightgray"
                        else:
                            # use most common function for the color
                            #if Global_State['genomebrowser_mode'] == "tree" \
                            #    or Global_State['genomebrowser_mode'] == "homologs":
                            #    if feature['annot'] in annot_repeat:
                            #        feature_element_color = color_names[sum([ord(c) for c in feature['annot']]) % len(color_names)]
                            #else:
                            #    feature_element_color = color_names[sum([ord(c) for c in feature['annot']]) % len(color_names)]

                            # color by most abundant function
                            most_abundant_fxn = None
                            most_abundant_cnt = 0
                            for fxn in sorted(feature['functions']): # must sort to avoid ties
                                #print ("checking fxn {}".format(fxn))  # DEBUG
                                if Global_State["function_abundance_counts"][fxn] > most_abundant_cnt:
                                    most_abundant_count = Global_State["function_abundance_counts"][fxn]
                                    most_abundant_fxn = fxn
                            #print ("most abundant fxn {}".format(most_abundant_fxn))  # DEBUG
                            feature_element_color = color_names[sum([ord(c) for c in most_abundant_fxn]) % len(color_names)]

                elif Global_State['genomebrowser_color_namespace'] == "ec":
                    if 'EC_number' in feature and feature['EC_number'] != '':
                        if Global_State['genomebrowser_mode'] == "tree" \
                            or Global_State['genomebrowser_mode'] == "homologs":
                            if feature['EC_number'] in annot_repeat:
                                feature_element_color = color_names[sum([ord(c) for c in feature['EC_number']]) % len(color_names)]
                        else:
                            feature_element_color = color_names[sum([ord(c) for c in feature['EC_number']]) % len(color_names)]

                elif Global_State['genomebrowser_color_namespace'] == 'cog' \
                    or Global_State['genomebrowser_color_namespace'] == 'pfam':

                    if Global_State['genomebrowser_mode'] == "genome":
                        contig_i = 0
                    else:
                        contig_i = feature_i
                    try:
                        domain_hits = Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[contig_i]]][ContigSet_names[contig_i]][feature['ID']]
                    except:
                        domain_hits = []

                    top_domfam = None                
                    for domhit in domain_hits:  # already reverse sorted by bitscore
                        domfam = domhit[DOMHIT_DOMFAM_I]
                        if (Global_State['genomebrowser_color_namespace'] == 'cog' and domfam[0:3] == 'COG') \
                            or (Global_State['genomebrowser_color_namespace'] == 'pfam' and domfam[0:2] == 'PF'):

                            if (Global_State['genomebrowser_mode'] == "tree" or Global_State['genomebrowser_mode'] == "homologs") \
                                  and not domfam in annot_repeat:
                                continue                            
                            top_domfam = domfam
                            break

                    if top_domfam != None:
                        if Global_State['genomebrowser_mode'] == "tree" \
                            or Global_State['genomebrowser_mode'] == "homologs":
                            if top_domfam in annot_repeat:  # already checked, but safer to repeat in case logic above changes
                                feature_element_color = color_names[sum([ord(c) for c in top_domfam]) % len(color_names)]
                        else:
                            feature_element_color = color_names[sum([ord(c) for c in top_domfam]) % len(color_names)]

                elif Global_State['genomebrowser_color_namespace'] == 'domains':
                    #feature_element_color = "none"
                    feature_element_color = "lightgray"
                    # individual domains will be drawn separately with rects below


            # edge color
            feature_element_edge_lw = 1
            feature_element_edge_color = "none"
        #    if (Global_State['genomebrowser_color_namespace'] == 'domains'):
        #        feature_element_edge_color = "black"
            if pivot_feature_flag:
                if Global_State['genomebrowser_mode'] == 'contigs':
                    feature_element_edge_color = "black"
                    feature_element_edge_lw = 1
                else:
                    feature_element_edge_color = "black"
                    feature_element_edge_lw = 2
            if feature_oi_flag:
                feature_element_edge_color = "black"
                feature_element_edge_lw = 2        


            # Draw Rectangle (patch) or Arrow (from pylab) for features
            #
            if feature['type'] == "CRISPR":
                new_arrow_l = ax.arrow (feature_element_start_pos+0.5*(window_end_pos-window_beg_pos), \
                                  y_pos_base, \
                                  -feature_element_base_length, \
                                  0.0, \
                                  width=0.75*arrow_w, \
                                  color=feature_element_color, \
                                  ec="none",
                                  #lw=feature_element_edge_lw, \
                                  alpha=0.75, \
                                  zorder=1, \
                                  picker=False, \
                                  **arrow_params)
                new_arrow_r = ax.arrow (feature_element_start_pos+0.5*(window_end_pos-window_beg_pos), \
                                  y_pos_base, \
                                  feature_element_base_length, \
                                  0.0, \
                                  width=0.75*arrow_w, \
                                  color=feature_element_color, \
                                  ec="none", \
                                  #lw=feature_element_edge_lw, \
                                  alpha=0.75, \
                                  zorder=1, \
                                  picker=False, \
                                  **arrow_params)
                # FIX: listeners do work, but seems excessive so turning off for now
                #feature_index = "%d,%d"%(feature_i,feature_j)
                #feature_artist_label_to_feature[feature_index] = feature
                #new_arrow_l.figure.canvas.mpl_connect ('pick_event', onpick_feature)
                #new_arrow_l.gid = feature_index
                #new_arrow_l.set_label("feature")
                #new_arrow_r.figure.canvas.mpl_connect ('pick_event', onpick_feature)  
                #new_arrow_r.gid = feature_index  # not sure I can have 2 elements with the same gid?
                new_arrow_r.set_label("feature")

            elif feature['type'] == "CRISPR spacer":
                rect_feature = Rectangle((x_pos, y_pos_base-0.325*arrow_w), \
                                         feature_element_base_length, 0.75*arrow_w, \
                                         facecolor=feature_element_color, ec="none", alpha=0.75, zorder=2, \
                                         picker=False)
                # FIX: listeners don't work as written below, but rect picker works elsewhere so should be doable
                #feature_index = "%d,%d"%(feature_i,feature_j)
                #feature_artist_label_to_feature[feature_index] = feature
                #rect_feature.figure.canvas.mpl_connect ('pick_event', onpick_feature)
                #rect_feature.gid = feature_index
                #rect_feature.set_label("feature")
                ax.add_patch(rect_feature)

            else:
                new_arrow = ax.arrow (feature_element_start_pos, \
                                  y_pos_base, \
                                  dir_adj*feature_element_base_length, \
                                  0.0, \
                                  width=arrow_w, \
                                  #color=feature_element_color, \
                                  #ec=feature_element_edge_color, \
                                  #lw=feature_element_edge_lw, \
                                  #alpha=0.75, \
                                  zorder=1, \
                                  picker=True, \
                                  **arrow_params)

                # Drop shadow baby!
                new_arrow.set_path_effects([path_effects.PathPatchEffect(offset=(drop_shadow_xshift, -1),
                                                                         facecolor='black', alpha=0.4),
                                            path_effects.PathPatchEffect(facecolor=feature_element_color,
                                                                         edgecolor=feature_element_edge_color,
                                                                         lw=feature_element_edge_lw)])


                # add popup onclick
                #
                # the object approach causes iPython to hang.  Use artist label to index into hashes instead
                #ax.add_patch (new_arrow)
                #active_arrow = InteractiveGenomeFeature (new_arrow, feature)
                #active_arrow.connect()
                #
                feature_index = "%d,%d"%(feature_i,feature_j)
                feature_artist_label_to_feature[feature_index] = feature    
                new_arrow.figure.canvas.mpl_connect ('pick_event', onpick_feature)
                new_arrow.gid = feature_index
                new_arrow.set_label("feature")


            # Add Domains
            #
            if Global_State['genomebrowser_color_namespace'] == 'domains':
                if Global_State['genomebrowser_mode'] == "genome":
                    contig_i = 0
                else:
                    contig_i = feature_i
                try:
                    domain_hits = Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[contig_i]]][ContigSet_names[contig_i]][feature['ID']]
                except:
                    domain_hits = []

                uncovered_frac_thresh = 0.5
                domhits_disp = []
                # needs to extend beyond by 1 because protein numbering starting from 1
                covered_mask = []
                for i in range(0,int(round(feature['end_pos']-feature['beg_pos']+1)/3)+2):
                    covered_mask.append(False)

                for dom_namespace in ['TIGR', 'PF', 'COG']:
        #        for dom_namespace in ['TIGR']:
        #        for dom_namespace in ['PF']:
        #        for dom_namespace in ['COG']:
                    for domhit in domain_hits:  # already reverse sorted by bitscore
                        domfam = domhit[DOMHIT_DOMFAM_I]
                        if domfam[0:len(dom_namespace)] != dom_namespace:
                            continue
                        if (Global_State['genomebrowser_mode'] == "tree" or Global_State['genomebrowser_mode'] == "homologs") \
                              and not domfam in annot_repeat:
                            continue

                        occluded = False
                        uncovered_pos_cnt = 0
                        for i in range(domhit[DOMHIT_BEG_I],domhit[DOMHIT_END_I]+1):
                            if not covered_mask[i]:
                                uncovered_pos_cnt += 1

                        if uncovered_pos_cnt/(domhit[DOMHIT_END_I]-domhit[DOMHIT_BEG_I]+1) < uncovered_frac_thresh:
                            occluded = True

                        if not occluded:
                            domhits_disp.append(domhit)
                            for i in range(domhit[DOMHIT_BEG_I],domhit[DOMHIT_END_I]+1):
                                covered_mask[i] = True

                dom_zorder = 3                    
                for domhit in sorted(domhits_disp, key=sort_by_domhit_len, reverse=True):
                    top_domfam = domhit[DOMHIT_DOMFAM_I]

                    if Global_State['genomebrowser_mode'] == "tree" \
                          or Global_State['genomebrowser_mode'] == "homologs":
                        if top_domfam in annot_repeat:  # already checked, but safer to repeat in case logic above changes
                            dom_element_color = color_names[sum([ord(c) for c in top_domfam]) % len(color_names)]
                    else:
                        dom_element_color = color_names[sum([ord(c) for c in top_domfam]) % len(color_names)]

                    # get coords
                    if direction == 'fwd':
                        dom_beg_pos = feature['beg_pos']+3*(domhit[DOMHIT_BEG_I]-1)
                        dom_end_pos = feature['beg_pos']+3*(domhit[DOMHIT_END_I]-1)
                    else:
                        dom_beg_pos = feature['end_pos']-3*(domhit[DOMHIT_END_I]-1)
                        dom_end_pos = feature['end_pos']-3*(domhit[DOMHIT_BEG_I]-1)

                    dom_window_beg_pos = disp_coord_transform (dom_beg_pos, \
                                                   pivot_pos, \
                                                   pivot_strand, \
                                                   Global_State['genomebrowser_window_bp_width'], \
                                                   Global_State['genomebrowser_xshift'], \
                                                   track_xshift, \
                                                   contig_mode_xshift)
                    dom_window_end_pos = disp_coord_transform (dom_end_pos, \
                                                   pivot_pos, \
                                                   pivot_strand, \
                                                   Global_State['genomebrowser_window_bp_width'], \
                                                   Global_State['genomebrowser_xshift'], \
                                                   track_xshift, \
                                                   contig_mode_xshift)

                    if dom_window_end_pos <= left_margin or dom_window_beg_pos >= 1.0-right_margin:
                        continue
                    if dom_window_beg_pos < left_margin:
                        dom_window_beg_pos = left_margin
                    if dom_window_end_pos > 1.0-right_margin:
                        dom_window_end_pos = 1.0-right_margin

                    # rect representation  # DEBUG
                    """
                    ##y_shift_tweak = -0.0005  # rectangles and arrows don't necessarily wind up on same pixel for close calls
                    y_shift_tweak = 0.0  # rectangles and arrows don't necessarily wind up on same pixel for close calls
                    dom_x_pos = dom_window_beg_pos
                    dom_y_pos = y_pos_base - 0.5*arrow_w + y_shift_tweak
                    dom_w = dom_window_end_pos - dom_window_beg_pos + .00000001
                    dom_h = arrow_w

                    rect_feature = Rectangle((dom_x_pos, dom_y_pos), dom_w, dom_h, \
                                             facecolor=dom_element_color, ec="black", lw=1, alpha=1.0, zorder=dom_zorder, \
                                             picker=False)  # Probably will want picker
                    ax.add_patch(rect_feature)
                    """

                    # arrow representation of domains
                    #
                    dom_y_scale = 0.9
                    #dom_y_scale = 1.0  # DEBUG
                    dom_arrow_w = dom_y_scale * arrow_w
                    dom_element_edge_color = 'none'
                    dom_element_edge_lw = 0
                    add_tip_crop = False

                    # direction and position            
                    if direction == 'fwd':

                        # domain completely in arrow head
                        if window_end_pos - dom_window_beg_pos <= head_l:
                            dom_element_base_length = .00000001
                            dom_head_l = window_end_pos - dom_window_beg_pos
                            if window_end_pos - window_beg_pos <= head_l:
                                this_head_l = window_end_pos - window_beg_pos
                            else:
                                this_head_l = head_l
                            dom_head_w = dom_y_scale * dom_head_l * (head_w/this_head_l)
                            dom_arrow_params = {'head_width':dom_head_w, 'head_length':dom_head_l}
                            if dom_window_end_pos < window_end_pos:
                                add_tip_crop = True
                                tip_start_pos = dom_window_end_pos
                                tip_base_length = .00000001
                                tip_head_l = window_end_pos - dom_window_end_pos
                                tip_head_w = tip_head_l * (head_w/this_head_l)
                                tip_arrow_params = {'head_width':tip_head_w, 'head_length':tip_head_l}

                        # domain completely *not* in arrow head
                        elif dom_window_end_pos <= window_end_pos - head_l:
                            dom_element_base_length = dom_window_end_pos - dom_window_beg_pos
                            dom_head_l = 0
                            dom_head_w = 0
                            dom_arrow_params = {'head_width':dom_head_w, 'head_length':dom_head_l}

                        # domain in both regions
                        else:  # dom_window_beg_pos <= window_end_pos - head_l and dom_window_end_pos > window_end_pos - head_l
                            dom_element_base_length = window_end_pos - dom_window_beg_pos - head_l
                            if dom_element_base_length < 0:
                                dom_element_base_length = .00000001
                            dom_head_l = head_l
                            dom_head_w = dom_y_scale * head_w
                            dom_arrow_params = {'head_width':dom_head_w, 'head_length':dom_head_l}
                            if dom_window_end_pos < window_end_pos:
                                add_tip_crop = True
                                tip_start_pos = dom_window_end_pos
                                tip_base_length = .00000001
                                tip_head_l = window_end_pos - dom_window_end_pos
                                tip_head_w = tip_head_l * (head_w/head_l)
                                tip_arrow_params = {'head_width':tip_head_w, 'head_length':tip_head_l}

                        dom_element_start_pos = dom_window_beg_pos
                        dir_adj = 1.0

                    else:  # direction == "rev"  

                        # domain completely in arrow head
                        if dom_window_end_pos - window_beg_pos <= head_l:
                            dom_element_base_length = .00000001
                            dom_head_l = dom_window_end_pos - window_beg_pos
                            if window_end_pos - window_beg_pos <= head_l:
                                this_head_l = window_end_pos - window_beg_pos
                            else:
                                this_head_l = head_l
                            dom_head_w = dom_y_scale * dom_head_l * (head_w/this_head_l)
                            dom_arrow_params = {'head_width':dom_head_w, 'head_length':dom_head_l}
                            if dom_window_beg_pos > window_beg_pos:
                                add_tip_crop = True
                                tip_start_pos = dom_window_beg_pos
                                tip_base_length = .00000001
                                tip_head_l = dom_window_beg_pos - window_beg_pos
                                tip_head_w = tip_head_l * (head_w/this_head_l)
                                tip_arrow_params = {'head_width':tip_head_w, 'head_length':tip_head_l}

                        # domain completely *not* in arrow head
                        elif dom_window_beg_pos >= window_beg_pos + head_l:
                            dom_element_base_length = dom_window_end_pos - dom_window_beg_pos
                            dom_head_l = head_l
                            dom_head_w = 0
                            dom_arrow_params = {'head_width':dom_head_w, 'head_length':0}

                        # domain in both regions
                        else:  # dom_window_beg_pos <= window_end_pos - head_l and dom_window_end_pos > window_end_pos - head_l
                            dom_element_base_length = dom_window_end_pos - window_beg_pos - head_l
                            if dom_element_base_length < 0:
                                dom_element_base_length = .00000001
                            dom_head_l = head_l
                            dom_head_w = dom_y_scale * head_w
                            dom_arrow_params = {'head_width':dom_head_w, 'head_length':dom_head_l}
                            if dom_window_beg_pos > window_beg_pos:
                                add_tip_crop = True
                                tip_start_pos = dom_window_beg_pos
                                tip_base_length = .00000001
                                tip_head_l = dom_window_beg_pos - window_beg_pos
                                tip_head_w = tip_head_l * (head_w/head_l)
                                tip_arrow_params = {'head_width':tip_head_w, 'head_length':tip_head_l}

                        dom_element_start_pos = window_beg_pos+dom_element_base_length+dom_head_l
                        dir_adj = -1.0 


                    # draw domain element 
                    #
                    dom_arrow = ax.arrow (dom_element_start_pos, \
                                  y_pos_base, \
                                  dir_adj*dom_element_base_length, \
                                  0.0, \
                                  width=dom_arrow_w, \
                                  color=dom_element_color, \
                                  ec=dom_element_edge_color, \
                                  #lw=dom_element_edge_lw, \
                                  alpha=1.0, \
                                  zorder=dom_zorder, \
                                  picker=False, \
                                  **dom_arrow_params)
                    if add_tip_crop:
                        tip_arrow = ax.arrow (tip_start_pos, \
                                  y_pos_base, \
                                  dir_adj*tip_base_length, \
                                  0.0, \
                                  width=dom_arrow_w, \
                                  color=feature_element_color, \
                                  ec='none', \
                                  #lw=dom_element_edge_lw, \
                                  alpha=1.0, \
                                  zorder=dom_zorder+1, \
                                  picker=False, \
                                  **tip_arrow_params)

                    dom_zorder += 2

                # restore border for highlighted features           
                if pivot_feature_flag or feature_oi_flag:
                    border_arrow = ax.arrow (feature_element_start_pos, \
                                  y_pos_base, \
                                  dir_adj*feature_element_base_length, \
                                  0.0, \
                                  width=arrow_w, \
                                  color='none', \
                                  ec=feature_element_edge_color, \
                                  lw=feature_element_edge_lw, \
                                  alpha=1.0, \
                                  zorder=dom_zorder+1, \
                                  picker=False, \
                                  **arrow_params)

            # Feature label
            #
            if Global_State['genomebrowser_window_bp_width'] <= text_disp_window_bp_limit or pivot_feature_flag:

                label_color = "black"
                arrow_label_fontsize = base_arrow_label_fontsize        

                # indicate pivot feature
                if pivot_feature_flag and Global_State['genomebrowser_mode'] != "contigs":
                    #label_color = 'red'
                    label_color = 'black'            
                    arrow_label_fontsize = foi_arrow_label_fontsize
                    if this_label_show_top == False and this_label_show_bot == False:
                        this_label_show_top = True

                if feature_oi_flag:
                    label_color = highlight_color            
                    arrow_label_fontsize = foi_arrow_label_fontsize
                    if this_label_show_top == False and this_label_show_bot == False:
                        this_label_show_top = True


                        
                # label position

                # old label pos
                """
                #label_y_pos = bottom_margin + (1.0-bottom_margin-top_margin)*row_delta*(total_rows-row_n)

                if Global_State['genomebrowser_window_bp_width'] > text_disp_window_bp_limit or this_label_show_top:
                    vert_align = "bottom"
                    label_y_pos = y_pos_base + 0.9*text_yshift
                    if direction == "rev":
                        label_start_pos += .0015
                    #else:
                    #    label_start_pos += .0005
                else:
                    vert_align = "top"
                    label_y_pos = y_pos_base - 1.25*text_yshift
                    if direction == "rev":
                        label_start_pos += .0050
                    #else:
                    #    label_start_pos += .0005
                """

                # this is what we drew the arrow with, so use these coords
                """
                new_arrow = ax.arrow (feature_element_start_pos, \
                                  y_pos_base, \
                                  dir_adj*feature_element_base_length, \
                                  0.0, \
                                  width=arrow_w, \
                                  #color=feature_element_color, \
                                  #ec=feature_element_edge_color, \
                                  #lw=feature_element_edge_lw, \
                                  #alpha=0.75, \
                                  zorder=1, \
                                  picker=True, \
                                  **arrow_params)
                """

                # new label pos
                label_start_pos = x_pos
                label_y_pos = y_pos_base
                print ("ROW: {} y_pos: {}".format(row_n, label_y_pos))
                label_y_pos += 0.07 * float(row_n) / float(max_row_n)
                vert_align = "center"   # top,bottom,center,baseline,center_baseline
                
#                if Global_State['genomebrowser_window_bp_width'] > text_disp_window_bp_limit or this_label_show_top:
#                    vert_align = "bottom"
#                    label_y_pos += 0.9*text_yshift
#                    if direction == "rev":
#                        label_start_pos += .0015
#                else:
#                    vert_align = "top"
#                    label_y_pos -= 1.25*text_yshift
#                    if direction == "rev":
#                        label_start_pos += .0050

                # draw label
                #
                if pivot_feature_flag or \
                   feature_oi_flag or \
                   this_label_show_top or \
                   this_label_show_bot:
                    
                    feature_name_disp = gene_name  # maybe add more options here
                    if row_n == 0 or row_n == 1:
                    #if len(feature_name_disp) <= max_feature_disp_len:
                    
                        # label
                        feature_label = ax.text(label_start_pos, \
                                                label_y_pos, \
                                                feature_name_disp, \
                                                #row_n, \  # DEBUG
                                                verticalalignment=vert_align, \
                                                horizontalalignment="left", \
                                                transform=ax.transAxes, \
                                                color=label_color, \
                                                fontsize=arrow_label_fontsize, \
                                                zorder=1)

                #print (feature['source_species']+" "+feature['name']+" %f %f %f"%(arrow_start_pos, arrow_base_length, y_pos_base))

            # add ellipsis to indicate cropped feature
            #
            if left_cropped:
                ellipsis_spacing = 0.2 * left_margin
                x_center_base = left_margin - 2*ellipsis_spacing
                y_center_base = y_pos_base
                diameter = 0.20 * arrow_w
                ellipse_to_circle_scaling = 0.75*figure_width / (figure_height_scaling*(total_rows+1))        
                x_diameter = 1.0 * diameter
                y_diameter = ellipse_to_circle_scaling * diameter
                for i in range(0,3):
                    ellipse_center = (x_center_base-i*ellipsis_spacing, y_center_base)
                    ellipsis_dot_patch = Ellipse (ellipse_center, x_diameter, y_diameter, \
                                                 facecolor=feature_element_color, edgecolor=feature_element_color, \
                                                 lw=1, alpha=1.0, zorder=1)

                    ax.add_patch(ellipsis_dot_patch)

            if right_cropped:
                ellipsis_spacing = 0.2 * left_margin
                x_center_base = 1.0-right_margin + 2*ellipsis_spacing
                y_center_base = y_pos_base
                diameter = 0.20 * arrow_w
                ellipse_to_circle_scaling = 0.75*figure_width / (figure_height_scaling*(total_rows+1))        
                x_diameter = 1.0 * diameter
                y_diameter = ellipse_to_circle_scaling * diameter
                for i in range(0,3):
                    ellipse_center = (x_center_base+i*ellipsis_spacing, y_center_base)
                    ellipsis_dot_patch = Ellipse (ellipse_center, x_diameter, y_diameter, \
                                                 facecolor=feature_element_color, edgecolor=feature_element_color, \
                                                 lw=1, alpha=1.0, zorder=1)

                    ax.add_patch(ellipsis_dot_patch)


        # Paint methods
        #
        def draw_mode_panel (ax, contig_name, genomebrowser_mode, data_set_name):

            # Config
            field_name_fontsize=12
            field_val_fontsize=11
            field_color="white"
            value_color="lightsalmon"
            #shadow_color="magenta"
            shadow_color="blue"
            #shadow_color="gray"
            base_zorder=5

            # Draw background with rounded corners
            #
            bb = mtransforms.Bbox([[0.06, 0.10], [0.94, 0.89]])
            bb_shadow = mtransforms.Bbox([[0.07, 0.08], [0.95, 0.87]])

            p_fancy = FancyBboxPatch((bb.xmin, bb.ymin),
                                     abs(bb.width), abs(bb.height),
                                     boxstyle="round,pad=0.04",
                                     mutation_aspect=0.25*figure_width/top_nav_height,
                                     fc=nav_bg_color,
                                     ec='white',
                                     lw=2,
                                     zorder=base_zorder)
            ax.add_patch(p_fancy)
            """
            p_fancy_back = FancyBboxPatch((bb.xmin, bb.ymin),
                                     abs(bb.width), abs(bb.height),
                                     boxstyle="round,pad=0.06",
                                     mutation_aspect=0.25*figure_width/top_nav_height,
                                     fc='white',
                                     ec='black',
                                     lw=1,
                                     zorder=base_zorder-1)
            ax.add_patch(p_fancy_back)
            """
            p_fancy_shadow = FancyBboxPatch((bb_shadow.xmin, bb_shadow.ymin),
                                     abs(bb_shadow.width), abs(bb_shadow.height),
                                     boxstyle="round,pad=0.04",
                                     mutation_aspect=0.25*figure_width/top_nav_height,
                                     fc=shadow_color,
                                     ec=shadow_color,
                                     lw=1,
                                     alpha=0.5,
                                     zorder=base_zorder-2)
            ax.add_patch(p_fancy_shadow)

            # Tool title
            ax.text(0.05, \
                    0.9, \
                    tool_title, \
                    verticalalignment="top", \
                    horizontalalignment="left", \
                    #transform=ax_top_left.transAxes, \
                    color="white", \
                    style="normal", # style="oblique", \
                    fontsize=14,
                    zorder=base_zorder+1)

            # Genome
            name_disp = contig_name
            if KBase_backend:
                #[ws_id, genome_contig_id] = contig_name.split('/')
                #[genome_id, contig_id] = genome_contig_id.split(genome_contig_id_delim)
                [genome_ref, scaffold_id] = contig_name.split(genome_contig_id_delim)
                name_disp = Species_name_by_genome_ref[genome_ref]
            ax.text(0.30, 0.5, "Genome", verticalalignment="bottom", horizontalalignment="right", color=field_color, fontsize=field_name_fontsize, zorder=base_zorder+1)
            ax.text(0.35, 0.5, name_disp, verticalalignment="bottom", horizontalalignment="left", color=value_color, fontsize=field_val_fontsize, zorder=base_zorder+1)

            # Mode
            ax.text(0.30, 0.3, "Mode", verticalalignment="bottom", horizontalalignment="right", color=field_color, fontsize=field_name_fontsize, zorder=base_zorder+1)
            ax.text(0.35, 0.3, genomebrowser_mode, verticalalignment="bottom", horizontalalignment="left", color=value_color, fontsize=field_val_fontsize, zorder=base_zorder+1)

            # Dataset
            #ax.text(0.30, 0.1, "Data", verticalalignment="bottom", horizontalalignment="right", color=field_color, fontsize=field_name_fontsize, zorder=base_zorder+1)
            #ax.text(0.35, 0.1, data_set_name, verticalalignment="bottom", horizontalalignment="left", color=value_color, fontsize=field_val_fontsize, zorder=base_zorder+1)


        def draw_control_panel (ax, genomebrowser_zoom, genomebrowser_xshift, genomebrowser_color_namespace):

            # Config
            text_color = "white"
            base_zorder = 5
            #shadow_color = 'magenta'
            shadow_color = 'blue'
            #shadow_color = 'gray'

            # Draw background with rounded corners
            #
            bb = mtransforms.Bbox([[0.02, 0.11], [0.98, 0.88]])
            bb_shadow = mtransforms.Bbox([[0.0225, 0.09], [0.9825, 0.86]])

            p_fancy = FancyBboxPatch((bb.xmin, bb.ymin),
                                     abs(bb.width), abs(bb.height),
                                     boxstyle="round,pad=0.015",
                                     mutation_aspect=0.75*figure_width/top_nav_height,
                                     fc=nav_bg_color,
                                     ec='white',
                                     lw=2,
                                     zorder=base_zorder)
            ax.add_patch(p_fancy)
            """
            p_fancy_back = FancyBboxPatch((bb.xmin, bb.ymin),
                                     abs(bb.width), abs(bb.height),
                                     boxstyle="round,pad=0.020",
                                     mutation_aspect=0.75*figure_width/top_nav_height,
                                     fc='white',
                                     ec='black',
                                     lw=1,
                                     zorder=base_zorder-1)
            ax.add_patch(p_fancy_back)
            """
            p_fancy_shadow = FancyBboxPatch((bb_shadow.xmin, bb_shadow.ymin),
                                     abs(bb_shadow.width), abs(bb_shadow.height),
                                     boxstyle="round,pad=0.015",
                                     mutation_aspect=0.75*figure_width/top_nav_height,
                                     fc=shadow_color,
                                     ec=shadow_color,
                                     lw=1,
                                     alpha=0.5,
                                     zorder=base_zorder-2)
            ax.add_patch(p_fancy_shadow)


            # Color slider
            #
            slider_name = "Color"
            slider_name_h_align = "right"
            tic_label_v_align="bottom"
            #slider_x = 0.64
            slider_x = 0.13
            slider_y = 0.55
            #track_w = 0.25
            #track_w = 0.45
            track_w = 0.40
            slider_h = 0.2
            track_h  = 0.05
            knob_w = 0.025
            knob_h = 0.25
            slider_margin = 0.02
            tic_margin = 0.05
            tic_w = 0.004

            num_tics = len(color_namespace_names)
            step = track_w / (num_tics-1)
            tol = 0.1 * step
            knob_initial_pos = 0
            for i in range(0,num_tics):
                knob_initial_pos = i
                if Global_State['genomebrowser_color_namespace'] == color_namespace_names[i]:
                    break
            tic_labels = []
            for i in range (0,num_tics):
                tic_labels.append(color_namespace_names_disp[i])

            ax.text(slider_x-1.5*slider_margin, slider_y+0.5*slider_h, slider_name, verticalalignment="center", horizontalalignment=slider_name_h_align, color=text_color, fontsize=12, zorder=base_zorder+1)
            ax.add_patch(Rectangle((slider_x-slider_margin, slider_y), track_w+2*slider_margin, slider_h, facecolor="gray", alpha=0.5, zorder=base_zorder+2))
            ax.add_patch(Rectangle((slider_x, slider_y+0.5*slider_h-0.5*track_h), track_w, track_h, facecolor="black", alpha=1.0, zorder=base_zorder+4))
            for i in range (0,num_tics):
                tic_shift = 0
                if i == 0:
                    tic_shift = 0.5*tic_w
                elif i == num_tics-1:
                    tic_shift = -0.5*tic_w
                ax.add_patch(Rectangle((slider_x+i*step-0.5*tic_w+tic_shift, slider_y+tic_margin), tic_w, slider_h-2.0*tic_margin, facecolor="lightgray", edgecolor="lightgray", alpha=0.6, zorder=base_zorder+3))
                ax.text(slider_x+i*step+tic_shift, slider_y+slider_h+1.00*tic_margin, tic_labels[i], verticalalignment=tic_label_v_align, horizontalalignment="center", color=text_color, fontsize=10, zorder=base_zorder+1)

            #color_knob_patch = Rectangle((slider_x+0.5*track_w-0.5*knob_w, slider_y-0.5*(knob_h-slider_h)), knob_w, knob_h, facecolor="lightgray", alpha=1.0, zorder=4)
            color_knob_patch = Rectangle((slider_x+knob_initial_pos*step-0.5*knob_w, slider_y-0.5*(knob_h-slider_h)), knob_w, knob_h, facecolor="lightgray", alpha=1.0, zorder=base_zorder+5)
            ax.add_patch (color_knob_patch)
            color_knob = DraggableRectangleSlider(color_knob_patch, func='color', move='x', min_pos=slider_x-0.5*knob_w, max_pos=slider_x+track_w-0.5*knob_w, step=step, tol=tol)
            color_knob.connect()


            # Mode slider  
            #
            slider_name = "Mode"
            slider_name_h_align = "left"
            tic_label_v_align="bottom"
            #slider_x = 0.13
            slider_x = 0.64
            slider_y = 0.55
            #track_w = 0.45
            track_w = 0.25
            slider_h = 0.2
            track_h  = 0.05
            knob_w = 0.025
            knob_h = 0.25
            slider_margin = 0.02
            tic_margin = 0.05
            tic_w = 0.004
            num_tics = len(mode_names)
            step = track_w / (num_tics-1)
            tol = 0.1 * step
            knob_initial_pos = 0
            for i in range(0,num_tics):
                knob_initial_pos = i
                if Global_State['genomebrowser_mode'] == mode_names[i]:
                    break
            tic_labels = []
            for i in range (0,num_tics):
                tic_labels.append(mode_names_disp[i])

            ax.text(slider_x+track_w+1.5*slider_margin, slider_y+0.5*slider_h, slider_name, verticalalignment="center", horizontalalignment=slider_name_h_align, color=text_color, fontsize=12, zorder=base_zorder+1)
            ax.add_patch(Rectangle((slider_x-slider_margin, slider_y), track_w+2*slider_margin, slider_h, facecolor="gray", alpha=0.5, zorder=base_zorder+2))
            ax.add_patch(Rectangle((slider_x, slider_y+0.5*slider_h-0.5*track_h), track_w, track_h, facecolor="black", alpha=1.0, zorder=base_zorder+4))
            for i in range (0,num_tics):
                tic_shift = 0
                if i == 0:
                    tic_shift = 0.5*tic_w
                elif i == num_tics-1:
                    tic_shift = -0.5*tic_w
                ax.add_patch(Rectangle((slider_x+i*step-0.5*tic_w+tic_shift, slider_y+tic_margin), tic_w, slider_h-2.0*tic_margin, facecolor="lightgray", edgecolor="lightgray", alpha=0.6, zorder=base_zorder+3))
                ax.text(slider_x+i*step+tic_shift, slider_y+slider_h+1.00*tic_margin, tic_labels[i], verticalalignment=tic_label_v_align, horizontalalignment="center", color=text_color, fontsize=10, zorder=base_zorder+1)

            mode_knob_patch = Rectangle((slider_x+knob_initial_pos*step-0.5*knob_w, slider_y-0.5*(knob_h-slider_h)), knob_w, knob_h, facecolor="lightgray", alpha=1.0, zorder=base_zorder+5, transform=ax.transAxes)
            ax.add_patch (mode_knob_patch)
            mode_knob = DraggableRectangleSlider(mode_knob_patch, func='mode', move='x', min_pos=slider_x-0.5*knob_w, max_pos=slider_x+track_w-0.5*knob_w, step=step, tol=tol)
            mode_knob.connect()


            # Position slider  
            #
            slider_name = "Position"
            slider_name_h_align = "right"
            tic_label_v_align="top"
            pan_button_w = 0.025
            pan_button_h = 0.2
            pan_button_margin = 0.01
            slider_margin = 0.02
            slider_x = 0.13 + 3*(pan_button_w+pan_button_margin)
            slider_y = 0.25
            #track_w = 0.45 - 3*(pan_button_w+pan_button_margin)
            track_w = 0.40 - 3*(pan_button_w+pan_button_margin)
            slider_h = 0.2
            track_h  = 0.05
            knob_w = 0.025
            knob_h = 0.25
            tic_margin = 0.05
            tic_w = 0.004
            num_tics = 11
            step = track_w / (num_tics-1)
            tol = 0.1 * step
            label_bp_disp = 1000000
            genome_step = int(round(Global_State['PrimaryContig_len'] / (num_tics-1)))
            knob_initial_pos = (num_tics-1)*(1.0-(Global_State['PrimaryContig_len']-(Global_State['PrimaryAnchor_pivot_pos']+Global_State['genomebrowser_xshift']))/Global_State['PrimaryContig_len'])
            tic_labels = []
            for i in range (0,num_tics):
                if i == num_tics-1:
                    tic_val = "%5.1f"%(Global_State['PrimaryContig_len']/label_bp_disp)
                else:
                #    tic_val = i*genome_step/label_bp_disp      
                #tic_labels.append("%4.1f"%tic_val)
                     tic_val = "%5.1f"%(i*genome_step/label_bp_disp)        
                tic_labels.append(tic_val.strip())

            ax.text(slider_x-3*(pan_button_w+pan_button_margin)-1.5*slider_margin, slider_y+0.5*slider_h, slider_name, verticalalignment="center", horizontalalignment=slider_name_h_align, color=text_color, fontsize=12, zorder=base_zorder+1)
            ax.add_patch(Rectangle((slider_x-slider_margin, slider_y), track_w+2*slider_margin, slider_h, facecolor="gray", alpha=0.5, zorder=base_zorder+2))
            ax.add_patch(Rectangle((slider_x, slider_y+0.5*slider_h-0.5*track_h), track_w, track_h, facecolor="black", alpha=1.0, zorder=base_zorder+4))
            for i in range (0,num_tics):
                tic_shift = 0
                if i == 0:
                    tic_shift = 0.5*tic_w
                elif i == num_tics-1:
                    tic_shift = -0.5*tic_w
                ax.add_patch(Rectangle((slider_x+i*step-0.5*tic_w+tic_shift, slider_y+tic_margin), tic_w, slider_h-2.0*tic_margin, facecolor="lightgray", edgecolor="lightgray", alpha=0.6, zorder=base_zorder+3))
                ax.text(slider_x+i*step+tic_shift, slider_y-1.25*tic_margin, tic_labels[i], verticalalignment=tic_label_v_align, horizontalalignment="center", color=text_color, fontsize=10, zorder=base_zorder+1)

            pos_knob_patch = Rectangle((slider_x+knob_initial_pos*step-0.5*knob_w, slider_y-0.5*(knob_h-slider_h)), knob_w, knob_h, facecolor="lightgray", alpha=1.0, zorder=base_zorder+5, transform=ax.transAxes)
            ax.add_patch (pos_knob_patch)
            pos_knob = DraggableRectangleSlider(pos_knob_patch, func='pan', move='x', min_pos=slider_x-0.5*knob_w, max_pos=slider_x+track_w-0.5*knob_w, step=step, tol=tol)
            pos_knob.connect()

            # pan buttons
            pan_button_left = Rectangle((slider_x-3*(pan_button_w+pan_button_margin)-slider_margin, slider_y), pan_button_w, pan_button_h, facecolor="black", alpha=1.0, zorder=base_zorder+3, transform=ax.transAxes, picker=True)
            ax.text(slider_x-3*(pan_button_w+pan_button_margin)-slider_margin+0.5*pan_button_w, slider_y+0.5*pan_button_h, "<", verticalalignment="center", horizontalalignment="center", color="white", fontsize=12, zorder=base_zorder+4)
            pan_button_left.set_label("pan_button")
            pan_button_left.gid = "left_pan"
            ax.add_patch(pan_button_left)

            pan_button_reset = Rectangle((slider_x-2*(pan_button_w+pan_button_margin)-slider_margin, slider_y), pan_button_w, pan_button_h, facecolor="black", alpha=1.0, zorder=base_zorder+3, transform=ax.transAxes, picker=True)
            ax.text(slider_x-2*(pan_button_w+pan_button_margin)-slider_margin+0.5*pan_button_w, slider_y+0.5*pan_button_h, "o", verticalalignment="center", horizontalalignment="center", color="white", fontsize=14, zorder=base_zorder+4)
            pan_button_reset.set_label("pan_button")
            pan_button_reset.gid = "reset_pan"
            ax.add_patch(pan_button_reset)

            pan_button_right = Rectangle((slider_x-1*(pan_button_w+pan_button_margin)-slider_margin, slider_y), pan_button_w, pan_button_h, facecolor="black", alpha=1.0, zorder=base_zorder+3, transform=ax.transAxes, picker=True)
            ax.text(slider_x-1*(pan_button_w+pan_button_margin)-slider_margin+0.5*pan_button_w, slider_y+0.5*pan_button_h, ">", verticalalignment="center", horizontalalignment="center", color="white", fontsize=12, zorder=base_zorder+4)
            pan_button_right.set_label("pan_button")
            pan_button_right.gid = "right_pan"
            ax.add_patch(pan_button_right)

            pan_button_reset.figure.canvas.mpl_connect ('pick_event', onpick_pan)


            # Zoom slider
            #
            slider_name = "Zoom"
            slider_name_h_align = "left"
            tic_label_v_align="top"
            slider_x = 0.64
            slider_y = 0.25
            track_w = 0.25
            slider_h = 0.2
            track_h  = 0.05
            knob_w = 0.025
            knob_h = 0.25
            slider_margin = 0.02
            tic_margin = 0.05
            tic_w = 0.004
            num_tics = def_genomebrowser_zoom_tics
            step = track_w / (num_tics-1)
            tol = 0.1 * step
            knob_initial_pos = Global_State['genomebrowser_zoom']
            tic_labels = []
            for i in range (0,num_tics):
                tic_val = int(round(0.001*Global_State['def_genomebrowser_window_bp_width']*(2**i)))
                tic_labels.append("%d"%tic_val)

            ax.text(slider_x+track_w+1.5*slider_margin, slider_y+0.5*slider_h, slider_name, verticalalignment="center", horizontalalignment=slider_name_h_align, color=text_color, fontsize=12, zorder=base_zorder+1)
            ax.add_patch(Rectangle((slider_x-slider_margin, slider_y), track_w+2*slider_margin, slider_h, facecolor="gray", alpha=0.5, zorder=base_zorder+2))
            ax.add_patch(Rectangle((slider_x, slider_y+0.5*slider_h-0.5*track_h), track_w, track_h, facecolor="black", alpha=1.0, zorder=base_zorder+4))
            for i in range (0,num_tics):
                tic_shift = 0
                if i == 0:
                    tic_shift = 0.5*tic_w
                elif i == num_tics-1:
                    tic_shift = -0.5*tic_w
                ax.add_patch(Rectangle((slider_x+i*step-0.5*tic_w+tic_shift, slider_y+tic_margin), tic_w, slider_h-2.0*tic_margin, facecolor="lightgray", edgecolor="lightgray", alpha=0.6, zorder=base_zorder+3))
                ax.text(slider_x+i*step+tic_shift, slider_y-1.25*tic_margin, tic_labels[i], verticalalignment=tic_label_v_align, horizontalalignment="center", color=text_color, fontsize=10, zorder=base_zorder+1)

            zoom_knob_patch = Rectangle((slider_x+knob_initial_pos*step-0.5*knob_w, slider_y-0.5*(knob_h-slider_h)), knob_w, knob_h, facecolor="lightgray", alpha=1.0, zorder=base_zorder+5)
            ax.add_patch (zoom_knob_patch)
            zoom_knob = DraggableRectangleSlider(zoom_knob_patch, func='zoom', move='x', min_pos=slider_x-0.5*knob_w, max_pos=slider_x+track_w-0.5*knob_w, step=step, tol=tol)
            zoom_knob.connect()


        def draw_sidenav_panel (ax, genomebrowser_mode):

            # Tree representation
            #
            if genomebrowser_mode == "tree":

                """
                # BioPython Phylo version screws up height scaling
                treeObj = Phylo.read(tree_data_base_path+'/'+tree_data_file, tree_data_format)
                #print (treeObj)
                treeObj.ladderize()   # Flip branches so deeper clades are displayed at top
                # TODO: reorder genome tracks following species tree?  Must put primary contig on top.
                Phylo.draw(treeObj, axes=ax)
                """

                # Draw Tree with ETE3
                if GeneTree_ref:
                    tree_img_path = os.path.join (output_dir, GeneTree_obj_name+'.png')
                    newick_string = geneTree_data['tree']
                else:
                    tree_img_path = path.join('.','tree.png')
                    newick_string = ''
                    tree_path = path.join(tree_data_base_path,tree_data_file)
                    #tree_img_path = path.join('.','tree.pdf')  # can't put pdf into plot
                    #with open (tree_path, 'r', 0) as tree_file_handle:   # can't have unbuffered text I/O
                    with open (tree_path, 'r') as tree_file_handle:
                        for tree_line in tree_file_handle:
                            newick_string += tree_line
                self.log(console,"LOADED NEWICK")  # DEBUG

                # ETE3 customization
                import ete3
                self.log(console,"INSTANTIATE TREE")  # DEBUG
                treeObj = ete3.Tree(newick_string)
                treeObj.ladderize()  # read row order from leaves?
                ts = ete3.TreeStyle()
                self.log(console,"SET TREE STYLE")  # DEBUG
                #ts.show_leaf_name = True
                ts.show_leaf_name = False
                ts.show_branch_length = False
                ts.show_branch_support = True
                #ts.scale = 50 # 50 pixels per branch length unit
                ts.branch_vertical_margin = 70 # pixels between adjacent branches
                #ts.title.add_face(ete3.TextFace(params['output_name']+": "+params['desc'], fsize=10), column=0)
                node_style = ete3.NodeStyle()
                node_style["fgcolor"] = "#606060"  # for node balls
                node_style["size"] = 10  # for node balls (gets reset based on support)
                node_style["vt_line_color"] = "#606060"
                node_style["hz_line_color"] = "#606060"
                node_style["vt_line_width"] = 4
                node_style["hz_line_width"] = 4
                node_style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
                node_style["hz_line_type"] = 0

                leaf_style = ete3.NodeStyle()
                leaf_style["fgcolor"] = "#ffffff"  # for node balls
                leaf_style["size"] = 0  # for node balls (we're not using it to add space)
                leaf_style["vt_line_color"] = "#606060"  # unecessary
                leaf_style["hz_line_color"] = "#606060"
                leaf_style["vt_line_width"] = 4
                leaf_style["hz_line_width"] = 4
                leaf_style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
                leaf_style["hz_line_type"] = 0

                self.log(console,"TRAVERSING TREE")  # DEBUG
                for n in treeObj.traverse():
                    if n.is_leaf():
                        style = leaf_style
                        #n.name = genome_sci_name_by_id[genome_id]
                        #n.name = None
                        #n.add_face(ete3.TextFace(leaf_name_disp, fsize=10), column=0, position="branch-right")
                    else:
                        style = ete3.NodeStyle()
                        for k in node_style.keys():
                            style[k] = node_style[k]

                        if n.support > 0.95:
                            style["size"] = 12
                        elif n.support > 0.90:
                            style["size"] = 10
                        elif n.support > 0.80:
                            style["size"] = 8
                        else:
                            style["size"] = 4

                    n.set_style(style)

                # write tree to image file

                # this is what I had before
                """
                dpi = 600
                img_units = "in"
                img_pix_width = 200
                img_in_width = round(float(img_pix_width)/float(dpi), 1)
                #img_in_height = img_in_width * total_rows / 5.0
                #treeObj.render(tree_img_path, w=img_in_width, h=img_in_height, units=img_units, dpi=dpi, tree_style=ts)
                self.log(console,"SAVING TREE IMAGE")  # DEBUG
                treeObj.render(tree_img_path, w=img_in_width, units=img_units, dpi=dpi, tree_style=ts)
                """

                # new 
                dpi = 600
                img_units = "in"
                img_pix_width = 200
                img_in_width = round(float(img_pix_width)/float(dpi), 1)
                img_pix_height = 40 * total_rows
                #img_pix_height = 10 * total_rows
                img_in_height = round(float(img_pix_height)/float(dpi), 1)
                self.log(console,"SAVING TREE IMAGE")  # DEBUG
                treeObj.render(tree_img_path, w=img_in_width, h=img_in_height, units=img_units, dpi=dpi, tree_style=ts)
                #treeObj.render(tree_img_path, w=img_in_width, units=img_units, dpi=dpi, tree_style=ts)

                self.log(console,"TREE IMAGE SAVED")  # DEBUG

                # load and display tree
                from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
                xy = [0.5, 0.51]
                arr_img = pyplot.imread(tree_img_path, format='png')
                #ax.imshow(arr_img)
                self.log(console,"PLOT LOADED TREE")  # DEBUG

                imagebox = OffsetImage(arr_img, zoom=1.0)
                imagebox.image.axes = ax
                self.log(console,"AXES SET")  # DEBUG

                ab = AnnotationBbox(imagebox, xy, frameon=False)
                ax.add_artist(ab)
                self.log(console,"IMAGE BOX ADDED TO AXES")  # DEBUG


            # Line representation
            #
            elif genomebrowser_mode == "contigs" \
                or genomebrowser_mode == "homologs":

                contig_label_fontsize = 12
                contig_label_x_margin = 0.0075
                contig_label_y_margin = 0.06 / (figure_height_scaling*(total_rows+1))
                contig_color = "steelblue"
                contig_disp_margin = 0.10
                contig_h           = 0.17 / (figure_height_scaling*(total_rows+1))
                clickable_y_margin = 0.75 * contig_h


                longest_contig_len = 0
                for i,contig_len in enumerate(Global_State['Contig_lens']):
                    if i >= max_rows:
                        break
                    if longest_contig_len < contig_len:
                        longest_contig_len = contig_len

                for i,contig_len in enumerate(Global_State['Contig_lens']):
                    if i >= max_rows:
                        break
                    pivot_pos = Global_State['pivot_pos_list'][i]

                    # Contig coords
                    #
                    x0 = contig_disp_margin
                    y0 = bottom_margin + (1.0-bottom_margin-top_margin)*row_delta*(total_rows-(i+1)) - 0.5*contig_h
                    w = (1.0 - 2*contig_disp_margin) * contig_len / longest_contig_len
                    h = contig_h

                    # Add clickable base
                    #
                    x0_click = 0.5*contig_disp_margin
                    y0_click = y0 - 0.75*clickable_y_margin
                    w_click = (1.0 - 2*contig_disp_margin) * contig_len / longest_contig_len + contig_disp_margin
                    h_click = contig_h + 3.0*clickable_y_margin
                    clickable_contig_rect = Rectangle((x0_click,y0_click), w_click, h_click, \
                                                          facecolor="none", edgecolor="none", \
                                                          alpha=0.0, zorder=0, \
                                                          picker=True)         
                    clickable_contig_rect.set_label("linear_contig_nav")
                    clickable_contig_rect.gid = "%d,%f,%f"%(i,x0,w)
                    ax.add_patch(clickable_contig_rect)
                    clickable_contig_rect.figure.canvas.mpl_connect ('pick_event', onpick_contig_nav)

                    # Show contig
                    #
                    contig_rect = Rectangle((x0,y0), w, h, \
                                            facecolor=contig_color, edgecolor=contig_color, \
                                                alpha=1.0, zorder=1)                 
                    # drop shadow
                    contig_rect.set_path_effects([path_effects.PathPatchEffect(offset=(2, -2),
                                                                               facecolor='magenta', edgecolor='none', alpha=0.5),
                                                  path_effects.PathPatchEffect(facecolor=contig_color, edgecolor=contig_color)]) 
                    ax.add_patch (contig_rect)

                    # Add label
                    #
                    #print ("CONTIGSET_NAME: '"+ContigSet_names[i]+"'\n")  # DEBUG
                    disp_name = ContigSet_names[i]
                    if KBase_backend:
                        #full_name_split = ContigSet_names[i].split('/')
                        #[genome_id, scaffold_id] = full_name_split[1].split(genome_contig_id_delim)
                        [genome_id, scaffold_id] = ContigSet_names[i].split(genome_contig_id_delim)
                        scaffold_id_disp = scaffold_id.split('.')[-1]
                        disp_name = Species_name_by_genome_ref[genome_id]+' c:'+scaffold_id_disp
                    ax.text(x0+contig_label_x_margin, y0-contig_label_y_margin,
                            disp_name, verticalalignment="top", horizontalalignment="left",
                            color="black", fontsize=contig_label_fontsize, zorder=1)

                    # Mark pivot
                    #
                    if Global_State['genomebrowser_mode'] == 'homologs':
                        mark_h = 1.25 * contig_h
                        mark_w = 0.01
                        mark_margin = 0.02 * mark_h
                        disp_mark_x0 = x0 + w * (1-(contig_len-pivot_pos)/contig_len) - 0.5*mark_w
                        disp_mark_y0 = y0 + contig_h + mark_margin
                        pivot_feature_rect = Rectangle((disp_mark_x0,disp_mark_y0), mark_w, mark_h,
                                                       facecolor="black", edgecolor="none", alpha=1.0, zorder=2) 
                        ax.add_patch (pivot_feature_rect)

                    # Mark search results
                    #
                    try:
                        contig_search_results = search_results[i]
                    except:
                        contig_search_results = []
                    for j,result_list in enumerate(contig_search_results):
                        mark_color = search_color_names[j % len(search_color_names)]
                        for result in contig_search_results[j]:
                            mid_pos = 0.5 * (result['beg_pos'] + result['end_pos'])
                            mark_h = 1.25 * contig_h
                            mark_w = 0.01
                            mark_margin = 0.02 * mark_h
                            disp_mark_x0 = x0 + w * (1-(contig_len-mid_pos)/contig_len) - 0.5*mark_w
                            disp_mark_y0 = y0 + contig_h + mark_margin
                            search_feature_rect = Rectangle((disp_mark_x0,disp_mark_y0), mark_w, mark_h,
                                                            facecolor=mark_color, edgecolor="none", alpha=1.0, zorder=3) 
                            ax.add_patch (search_feature_rect)

                    # Indicate window viewed by track panel
                    #
                    contig_mode_xshift = 0.0
                    if Global_State['genomebrowser_mode'] == 'contigs':
                            #contig_mode_xshift = 0.5*Global_State['genomebrowser_window_bp_width'] - 0.5*(Feature_slices[i][0]['end_pos']-Feature_slices[i][0]['beg_pos'])
                        contig_mode_xshift = 0.5*Global_State['genomebrowser_window_bp_width']

                    #if (pivot_pos+Global_State['genomebrowser_xshift']) >= -.0000001 \
                    #    and (pivot_pos+Global_State['genomebrowser_xshift']) <= (contig_len + .0000001):
                    track_window_h = 1.6 * contig_h
                    track_window_w = w * Global_State['genomebrowser_window_bp_width'] / contig_len
                    track_window_x0 = x0 - 0.5*track_window_w + w * ((pivot_pos+Global_State['genomebrowser_xshift']+contig_mode_xshift) / contig_len)
                    track_window_y0 = y0 + 0.5*h - 0.5*track_window_h
                    if track_window_x0 > x0+w or track_window_x0 + track_window_w < x0:
                        continue                            
                    track_window_rect = Rectangle((track_window_x0,track_window_y0), track_window_w, track_window_h, \
                                                      facecolor="red", edgecolor="red", alpha=0.4, zorder=4) 
                    ax.add_patch (track_window_rect)


            # Circle representation
            #
            elif genomebrowser_mode == "genome":

                # set up circle genome representation
                contig_label_fontsize = 12
                contig_label_y_margin = 0.5 / (figure_height_scaling*(total_rows+1))
                contig_color = "steelblue"
                diameter = 0.75
                clickable_diameter = 1.25 * diameter
                genome_lw = 10
                marker_lw = 15
                track_window_lw = 20
                marker_len = 0.1
                GC_window_bp = 5000
                GC_bar_high_lw = 20
                GC_bar_high_perc = 0.20
                GC_base_diameter = 0.55
                GC_lw_to_coord_scale = 0.01
                ellipse_center_x = 0.50
                #ellipse_center_y = 0.70
                ellipse_center_y = (total_rows+1 - 2.5) / (total_rows+1)
                ellipse_center = (ellipse_center_x, ellipse_center_y)
                ellipse_to_circle_scaling = 0.25*figure_width / (figure_height_scaling*(total_rows+1))        


                # Add clickable base
                #
                x_diameter = 1.0 * clickable_diameter
                y_diameter = ellipse_to_circle_scaling * clickable_diameter
                #Ellipse(xy, width, height, angle=0.0, **kwargs)
                clickable_contig_image = Ellipse (ellipse_center, x_diameter, y_diameter, \
                                                  facecolor="none", edgecolor="none", alpha=0.0, zorder=0,
                                                  picker=True)
                clickable_contig_image.set_label("circular_contig_nav")
                clickable_contig_image.gid = "%d,%f,%f,%f"%(0,ellipse_center_x,ellipse_center_y,ellipse_to_circle_scaling)
                ax.add_patch(clickable_contig_image)
                clickable_contig_image.figure.canvas.mpl_connect ('pick_event', onpick_contig_nav)

                # Draw circle genome representation
                #
                x_diameter = 1.0 * diameter
                y_diameter = ellipse_to_circle_scaling * diameter
                #Ellipse(xy, width, height, angle=0.0, **kwargs)
                contig_image = Ellipse (ellipse_center, x_diameter, y_diameter, \
                                        facecolor="white", edgecolor=contig_color, lw=genome_lw, alpha=1.0, zorder=1,
                                       picker=False)
                ax.add_patch(contig_image)

                # Add label
                #
                disp_name = ContigSet_names[0]
                scaffold_id = None
                if KBase_backend:
                    #full_name_split = ContigSet_names[0].split('/')
                    #[genome_id, scaffold_id] = full_name_split[1].split(genome_contig_id_delim)
                    [genome_id, scaffold_id] = ContigSet_names[i].split(genome_contig_id_delim)
                    scaffold_id_disp = scaffold_id.split('.')[-1]
                    disp_name = Species_name_by_genome_ref[genome_id] + ' c:'+scaffold_id_disp
                ax.text(ellipse_center_x, ellipse_center_y - 0.5*y_diameter - contig_label_y_margin,
                        disp_name, verticalalignment="top", horizontalalignment="center",
                        color="black", fontsize=contig_label_fontsize)
                if scaffold_id != None:
                    ax.text(ellipse_center_x, ellipse_center_y - 0.5*y_diameter - 1.5*contig_label_y_margin,
                            scaffold_id, verticalalignment="top", horizontalalignment="center",
                            color="black", fontsize=contig_label_fontsize)            

                # Mark contig edge
                #
                mark_width = 1.0
                arc_end = 90 + 0.5*mark_width
                arc_beg = 90 - 0.5*mark_width
                contig_edge_arc = Arc (ellipse_center, x_diameter, y_diameter, \
                                        theta1=arc_beg, theta2=arc_end, \
                                        edgecolor="white", lw=genome_lw, alpha=1.0, zorder=2)
                ax.add_patch (contig_edge_arc)        

                # Mark pivot
                #
                mark_width = 1.0
                mark_pad = 0.125
                mark_x_diameter = x_diameter + 1.0 * mark_pad
                mark_y_diameter = y_diameter + ellipse_to_circle_scaling * mark_pad
                # pivot feature
                arc_end = 90 + 0.5*mark_width - 360 * (Global_State['PrimaryAnchor_pivot_pos'] / Global_State['PrimaryContig_len'])
                arc_beg = 90 - 0.5*mark_width - 360 * (Global_State['PrimaryAnchor_pivot_pos'] / Global_State['PrimaryContig_len'])
                pivot_feature_arc = Arc (ellipse_center, mark_x_diameter, mark_y_diameter, \
                                            theta1=arc_beg, theta2=arc_end, \
                                            edgecolor="black", lw=marker_lw, alpha=1.0, zorder=3)
                ax.add_patch (pivot_feature_arc)

                # Mark search features
                #
                try:
                    contig_search_results = search_results[0]
                except:
                    contig_search_results = []
                for j,result_list in enumerate(contig_search_results):
                    mark_color = search_color_names[j % len(search_color_names)]
                    for result in contig_search_results[j]:
                        mid_pos = 0.5 * (result['beg_pos'] + result['end_pos'])
                        arc_end = 90 + 0.5*mark_width - 360 * (mid_pos / Global_State['PrimaryContig_len'])
                        arc_beg = 90 - 0.5*mark_width - 360 * (mid_pos / Global_State['PrimaryContig_len'])
                        search_result_arc = Arc (ellipse_center, mark_x_diameter, mark_y_diameter, \
                                            theta1=arc_beg, theta2=arc_end, \
                                            edgecolor=mark_color, lw=marker_lw, alpha=1.0, zorder=4)
                        ax.add_patch (search_result_arc)

                # Indicate window viewed by track panel
                #
                track_extra = 0.5*(total_rows-1)*Global_State['genomebrowser_window_bp_width']
                arc_end = 90 - 360 * ((Global_State['PrimaryAnchor_pivot_pos']+Global_State['genomebrowser_xshift']-track_extra-0.5*Global_State['genomebrowser_window_bp_width']) / Global_State['PrimaryContig_len'])
                arc_beg = 90 - 360 * ((Global_State['PrimaryAnchor_pivot_pos']+Global_State['genomebrowser_xshift']+track_extra+0.5*Global_State['genomebrowser_window_bp_width']) / Global_State['PrimaryContig_len']) 
                if arc_end > 90:
                    arc_end = 90 - 0.5
                if arc_beg < -270:
                    arc_beg = -270 + 0.5

                track_window_arc = Arc (ellipse_center, x_diameter, y_diameter, \
                                        theta1=arc_beg, theta2=arc_end, \
                                        edgecolor="red", lw=track_window_lw, alpha=0.40, zorder=5)  # facecolor does nothing (no fill for Arc)
                ax.add_patch (track_window_arc)

                # Add GC
                #       
                GC_bar_scale = GC_bar_high_lw / GC_bar_high_perc
                GC_steps = int(round(Global_State['PrimaryContig_len'] / GC_window_bp))
                for GC_i in range(1,GC_steps-2):
                    arc_end = 90 - 360 * (GC_i*GC_window_bp / Global_State['PrimaryContig_len'])
                    arc_beg = 90 - 360 * ((GC_i+1)*GC_window_bp / Global_State['PrimaryContig_len'])
                    GC_in_window = compute_GC (Global_Genbank_Genomes[0].seq[(GC_i-1)*GC_window_bp:(GC_i+2)*GC_window_bp])
                    #GC_in_window = 0.5  # DEBUG
                    GC_delta = GC_in_window - Global_State['PrimaryContig_GCavg']
                    if (GC_delta > 0):
                        GC_color = "blue"
                        GC_bar_lw = int(round(GC_bar_scale * GC_delta))
                        GC_bar_diameter = GC_base_diameter + 0.5*GC_bar_lw*GC_lw_to_coord_scale
                        GC_x_diameter = 1.0 * GC_bar_diameter
                        GC_y_diameter = ellipse_to_circle_scaling * GC_bar_diameter
                        GC_arc = Arc (ellipse_center, GC_x_diameter, GC_y_diameter, \
                                        theta1=arc_beg, theta2=arc_end, \
                                        edgecolor=GC_color, lw=GC_bar_lw, alpha=1.0, zorder=1)  # facecolor does nothing (no fill for Arc)
                        ax.add_patch (GC_arc)
                    else:
                        GC_color = "red"
                        GC_bar_lw = int(round(GC_bar_scale * -GC_delta))
                        GC_bar_diameter = GC_base_diameter - 0.5*GC_bar_lw*GC_lw_to_coord_scale
                        GC_x_diameter = 1.0 * GC_bar_diameter
                        GC_y_diameter = ellipse_to_circle_scaling * GC_bar_diameter
                        GC_arc = Arc (ellipse_center, GC_x_diameter, GC_y_diameter, \
                                        theta1=arc_beg, theta2=arc_end, \
                                        edgecolor=GC_color, lw=GC_bar_lw, alpha=1.0, zorder=1)  # facecolor does nothing (no fill for Arc)
                        ax.add_patch (GC_arc)


        def draw_genomebrowser_panel (ax, \
                                       ContigSet_names, \
                                       PivotFeatures_IDs, \
                                       genomebrowser_window_bp_width, \
                                       genomebrowser_xshift, \
                                       genomebrowser_mode):                  

            # Read Domains
            #
            if Global_State['genomebrowser_color_namespace'] == 'domains' \
                or Global_State['genomebrowser_color_namespace'] == 'cog' \
                or Global_State['genomebrowser_color_namespace'] == 'pfam':

                # update Global_Domains, but really should be method of Domains object
                getDomainHits (ContigSet_names, \
                               genomebrowser_mode=genomebrowser_mode, \
                               domain_data_format=domain_data_format)


            # Read Features
            #
            if KBase_backend:
                Feature_slices = getFeatureSlicesKBase (ContigSet_names, \
                                                            PivotFeatures_IDs, \
                                                            genomebrowser_mode=genomebrowser_mode, \
                                                            genome_data_format=genome_data_format, \
                                                            window_size=genomebrowser_window_bp_width, \
                                                            genomebrowser_xshift=genomebrowser_xshift)
            elif genome_data_format == 'Genbank':
                Feature_slices = getFeatureSlicesGenbank (ContigSet_names, \
                                                            PivotFeatures_IDs, \
                                                            genomebrowser_mode=genomebrowser_mode, \
                                                            genome_data_format=genome_data_format, \
                                                            window_size=genomebrowser_window_bp_width, \
                                                            genomebrowser_xshift=genomebrowser_xshift)


            # Add features elements to figure
            #
            if genomebrowser_mode == "contigs" \
                or genomebrowser_mode == "genome" \
                or genomebrowser_mode == "homologs" \
                or genomebrowser_mode == "tree":

                # Determine repeats
                #
                annot_repeat = set()
                annot_seen   = set()
                seq_repeat   = set()
                seq_seen     = set()
                for i in range(0,len(Feature_slices)):
                    if i >= max_rows:
                        break

                    for j in range(1,len(Feature_slices[i])):  # skip [0]. it's the pivot and if in window will be repeated
                        if Feature_slices[i][j]['type'] == "CDS":
                            if Global_State['genomebrowser_color_namespace'] == 'annot' \
                                and 'annot' in Feature_slices[i][j] \
                                and Feature_slices[i][j]['annot'] != '':
                                if Feature_slices[i][j]['annot'] in annot_seen:
                                    annot_repeat.add(Feature_slices[i][j]['annot'])
                                else:
                                    annot_seen.add(Feature_slices[i][j]['annot'])

                            elif Global_State['genomebrowser_color_namespace'] == 'ec' \
                                and 'EC_number' in Feature_slices[i][j] \
                                and Feature_slices[i][j]['EC_number'] != '':
                                if Feature_slices[i][j]['EC_number'] in annot_seen:
                                    annot_repeat.add(Feature_slices[i][j]['EC_number'])
                                else:
                                    annot_seen.add(Feature_slices[i][j]['EC_number'])

                            elif Global_State['genomebrowser_color_namespace'] == 'cog' \
                                or Global_State['genomebrowser_color_namespace'] == 'pfam' \
                                or Global_State['genomebrowser_color_namespace'] == 'domains':

                                if Global_State['genomebrowser_mode'] == "genome":
                                    contig_i = 0
                                else:
                                    contig_i = i
                                try:
                                    domain_hits = Global_Domains[Genome_ref_by_Contig_name[ContigSet_names[contig_i]]][ContigSet_names[contig_i]][Feature_slices[i][j]['ID']]
                                except:
                                    domain_hits = []
                                for domhit in domain_hits:
                                    domfam = domhit[DOMHIT_DOMFAM_I]
                                    if domfam in annot_seen:
                                        annot_repeat.add(domfam)
                                    else:
                                        annot_seen.add(domfam)

                        elif (Feature_slices[i][j]['type'] == "CRISPR" or Feature_slices[i][j]['type'] == "CRISPR spacer") \
                            and 'dna_seq' in Feature_slices[i][j]:
                            if Feature_slices[i][j]['dna_seq'] in seq_seen:
                                seq_repeat.add(Feature_slices[i][j]['dna_seq'])
                            else:
                                seq_seen.add(Feature_slices[i][j]['dna_seq'])


                # Draw feature elements
                #
                for i in range(0,len(Feature_slices)):
                    if i >= max_rows:
                        break

                    pivot_pos = 0.5 * (Feature_slices[i][0]['beg_pos']+Feature_slices[i][0]['end_pos'])
                    pivot_strand = Feature_slices[i][0]['strand']
                    track_xshift = 0
                    if genomebrowser_mode == "genome":
                        #track_xshift = genomebrowser_window_bp_width * (i - (total_rows-1)/2)  # e.g. total_rows=7 -> -3,-2,-1,0,1,2,3
                        track_xshift = genomebrowser_window_bp_width * (i - (Global_State['genome_mode_n_rows']-1)/2)  # e.g. n_rows=7 -> -3,-2,-1,0,1,2,3
                        #print ("%d %f"%(i,track_xshift))
                    contig_mode_xshift = 0
                    if genomebrowser_mode == "contigs":
                        contig_mode_xshift = 0.5*genomebrowser_window_bp_width - 0.5*(Feature_slices[i][0]['end_pos']-Feature_slices[i][0]['beg_pos'])
                    prev_top_beg_pos = 0
                    prev_bot_beg_pos = 0
                    prev_top_name = ''
                    prev_bot_name = ''

                    for j in range(1,len(Feature_slices[i])):  # skip [0]. it's the pivot and if in window will be repeated
                        pivot_feature_flag = False
                        if Feature_slices[i][j]['beg_pos'] == Feature_slices[i][0]['beg_pos'] and \
                           Feature_slices[i][j]['end_pos'] == Feature_slices[i][0]['end_pos']:
                            pivot_feature_flag = True

                        # show label on top or bottom, if it fits
                        this_label_show_top = True
                        this_label_show_bot = False
                        this_beg_pos = disp_coord_transform (Feature_slices[i][j]['beg_pos'], \
                                                             pivot_pos, \
                                                             pivot_strand, \
                                                             genomebrowser_window_bp_width, \
                                                             genomebrowser_xshift, \
                                                             track_xshift, \
                                                             contig_mode_xshift)
                        #if j == 1:
                        if this_beg_pos < left_margin:
                            this_beg_pos = left_margin
                        elif len(prev_top_name) > arrow_label_scaling*(this_beg_pos-prev_top_beg_pos):
                            this_label_show_top = False
                        if this_label_show_top == False and len(prev_bot_name) <= arrow_label_scaling*(this_beg_pos-prev_bot_beg_pos):
                            this_label_show_bot = True

                        # FIX
                        max_row_n = total_rows  
                        if num_genomes < total_rows:
                            max_row_n = num_genomes

                        # finally draw arrow!  (or other shape)
                        draw_feature_element (ax, \
                                         i+1, \
                                         max_row_n, \
                                         i, \
                                         j, \
                                         Feature_slices[i][j], \
                                         pivot_pos, \
                                         pivot_strand, \
                                         track_xshift, \
                                         contig_mode_xshift, \
                                         pivot_feature_flag, \
                                         annot_repeat, \
                                         seq_repeat, \
                                         this_label_show_top, this_label_show_bot)

                        if this_label_show_top:
                            prev_top_beg_pos = this_beg_pos
                            prev_top_name = Feature_slices[i][j]['name']
                        if this_label_show_bot:
                            prev_bot_beg_pos = this_beg_pos
                            prev_bot_name = Feature_slices[i][j]['name']



        # Panel updates
        #
        def update_mode_panel (ax):
            ax.cla()
            draw_mode_panel (ax, \
                             Global_State['ContigSet_names'][0], \
                             Global_State['genomebrowser_mode'], \
                             Global_State['Dataset_names_list'][0])    
            if Display_plot:
                pyplot.show()

        def update_control_panel (ax):
            ax.cla()
            draw_control_panel (ax, \
                                Global_State['genomebrowser_zoom'], \
                                Global_State['genomebrowser_xshift'], \
                                Global_State['genomebrowser_color_namespace'])    
            if Display_plot:
                pyplot.show()

        def update_sidenav_panel (ax):
            ax.cla()
            draw_sidenav_panel (ax, \
                                Global_State['genomebrowser_mode'])
            if Display_plot:
                pyplot.show()

        def update_genomebrowser_panel (ax):

            # clear panel
            ax.cla()

            # Update contig order
            #   might be used elsewhere later, likely when track ordering changes
            #
            Contig_order = []
            Contig_order_lookup = {}
            for i,contig_name in enumerate(ContigSet_names):
                (genome_id,scaffold_id) = contig_name.split(genome_contig_id_delim)        
                try:
                    col = Contig_order_lookup['genome_id']
                except:
                    Contig_order_lookup['genome_id'] = {}
                Contig_order_lookup['genome_id']['scaffold_id'] = i
                Contig_order.append(contig_name)

            # Reset feature and popup_box storage
            #
            feature_artist_label_to_feature = {}
            feature_artist_label_to_popup_box = {}
            Global_State['popup_zorder'] = def_popup_zorder

            self.log(console, "GOT TO G.1")  # DEBUG

            # Add Genome tracks
            #
            draw_genomebrowser_panel (ax, \
                                      ContigSet_names = Global_State['ContigSet_names'], \
                                      PivotFeatures_IDs = Global_State['PivotFeatures_IDs'], \
                                      genomebrowser_window_bp_width = Global_State['genomebrowser_window_bp_width'], \
                                      genomebrowser_xshift = Global_State['genomebrowser_xshift'], \
                                      genomebrowser_mode = Global_State['genomebrowser_mode'])
            self.log(console, "GOT TO G.2")  # DEBUG
            if Display_plot:
                pyplot.show()
            self.log(console, "GOT TO G.3")  # DEBUG


        ###############################################################################                
        # Main
        ###############################################################################

        # Instantiate fig_KGB
        #   
        if not Hide_controls:
            figure_height = top_nav_height + figure_height_scaling*(total_rows+1)
        else:
            figure_height = figure_height_scaling*(total_rows+1)
        fig_KGB = pyplot.figure(1, (figure_width, figure_height))
        #fig_KGB.suptitle(tool_title, fontsize=20)

        if not Hide_controls:
            # want to maintain constant top nav height, so adjust proportion according to number of rows
            # e.g. with 6 tracks, top nav should be 2 units, so top rowspan=2 when bottom rowspan = 6 and grid is 8
            top_nav_rows = int(round(1000*top_nav_height/(top_nav_height+figure_height_scaling*(total_rows+1))))

            # subplot2grid(shape, loc, rowspan=1, colspan=1)
            ax_top_left   = pyplot.subplot2grid((1000,4),            (0,0), rowspan=top_nav_rows,        colspan=1)
            ax_top_center = pyplot.subplot2grid((1000,4),            (0,1), rowspan=top_nav_rows,        colspan=3)
            ax_left       = pyplot.subplot2grid((1000,4), (top_nav_rows,0), rowspan=(1000-top_nav_rows), colspan=1)
            ax_center     = pyplot.subplot2grid((1000,4), (top_nav_rows,1), rowspan=(1000-top_nav_rows), colspan=3)
        else:
            ax_left       = pyplot.subplot2grid((1000,4), (0,0), rowspan=1000, colspan=1)
            ax_center     = pyplot.subplot2grid((1000,4), (0,1), rowspan=1000, colspan=3)   

        # Let's turn off visibility of all tick labels and boxes here
        for ax in fig_KGB.axes:
            ax.xaxis.set_visible(False)  # remove axis labels and ticks
            ax.yaxis.set_visible(False)
            for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tick labels
                t.set_visible(False)
        if not Hide_controls:
            all_axes = [ax_top_left, ax_top_center, ax_left, ax_center]
            #for ax in [ax_top_left, ax_top_center]:
            #    ax.set_axis_bgcolor(nav_bg_color)
        else:
            all_axes = [ax_left, ax_center]
            
        for ax in all_axes:
            ax.spines['top'].set_visible(False) # Get rid of top axis line
            ax.spines['bottom'].set_visible(False) #  bottom axis line
            ax.spines['left'].set_visible(False) #  Get rid of bottom axis line
            ax.spines['right'].set_visible(False) #  Get rid of bottom axis line

        fig_KGB.tight_layout()  # left justify and space subplots reasonably.  Must follow axis creation


        # Load domain family descriptions
        #
        #domain_family_desc = read_domain_family_desc ()  # these need to be global (kludge)
        if domain_family_desc_base_path != None:
            read_domain_family_desc ()  # these need to be global (kludge)


        # Draw genome browser panel (must happen first to load Global_State)
        #
        self.log(console,"GOT TO G")  # DEBUG
        feature_artist_label_to_feature = {}     # these need to be global (kludge)
        feature_artist_label_to_popup_box = {}   # these need to be global (kludge)
        pos_knob_patch = None                    # these need to be global (kludge)
        update_genomebrowser_panel (ax_center)


        # Draw sidenav (tree, circles, etc.)
        #
        self.log(console,"GOT TO H")  # DEBUG
        update_sidenav_panel (ax_left)

        self.log(console,"GOT TO I")  # DEBUG

        # Draw mode panel
        #
        if not Hide_controls:
            update_mode_panel (ax_top_left)

        # Draw control panel (currently, must occur after primary genome read in update_genomebrowser_panel)
        if not Hide_controls:
            update_control_panel (ax_top_center)

        ###############################################################################
        # End KGB
        ###############################################################################


        # Save image
        #
        self.log(console, "SAVING REPORT IMAGES")
        img_dpi = 200
        #pyplot.show()
        self.log(console, "SAVING IMAGE")
        out_file_basename = GeneTree_obj_name+'-genome_context'
        png_file = out_file_basename+'.png'
        pdf_file = out_file_basename+'.pdf'
        output_png_file_path = os.path.join(html_dir, png_file);
        output_pdf_file_path = os.path.join(html_dir, pdf_file);
        fig_KGB.savefig(output_png_file_path, dpi=img_dpi)
        fig_KGB.savefig(output_pdf_file_path, format='pdf')


        # Create HTML
        #
        self.log(console, "SAVING REPORT HTML")
        img_html_width=1000
        html_file = 'index.html'
        output_html_file_path = os.path.join(html_dir, html_file);
        html_report_lines = []
        html_report_lines += ['<html>']
        html_report_lines += ['<head><title>KBase Gene Tree Genome Context</title></head>']
        html_report_lines += ['<body bgcolor="white">']
        html_report_lines += ['<img width=' + str(img_html_width) + ' src="' + png_file + '">']
        html_report_lines += ['</body>']
        html_report_lines += ['</html>']

        html_report_str = "\n".join(html_report_lines)
        with open(output_html_file_path, 'w') as html_handle:
            html_handle.write(html_report_str)
            
        # upload images and html
        dfu = DFUClient(self.callback_url, token=session_token)
        try:
            png_upload_ret = dfu.file_to_shock({'file_path': output_png_file_path,
                                                #'pack': 'zip'})
                                                'make_handle': 0})
        except:
            raise ValueError('error uploading png file to shock')
        try:
            pdf_upload_ret = dfu.file_to_shock({'file_path': output_pdf_file_path,
                                                #'pack': 'zip'})
                                                'make_handle': 0})
        except:
            raise ValueError('error uploading pdf file to shock')
        try:
            html_upload_ret = dfu.file_to_shock({'file_path': html_dir,
                                                 'make_handle': 0,
                                                 'pack': 'zip'})
        except:
            raise ValueError('error uploading png file to shock')
        

        # Create report obj
        #
        self.log(console, "BUILDING REPORT OBJECT")
        reportName = GeneTree_obj_name+'-genome_context_report_' + str(uuid.uuid4())
        #report += output_newick_buf+"\n"
        reportObj = {'objects_created': [],
                     'direct_html_link_index': 0,
                     'file_links': [],
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }
        reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                    'name': html_file,
                                    'label': GeneTree_obj_name+' genome context' + ' HTML'
                                    }
                                   ]        
        reportObj['file_links'].extend([{'shock_id': png_upload_ret['shock_id'],
                                         'name': GeneTree_obj_name+'-genome_context' + '.png',
                                         'label': GeneTree_obj_name+' genome context' + ' PNG'
                                         },
                                        {'shock_id': pdf_upload_ret['shock_id'],
                                         'name': GeneTree_obj_name+'-genome_context' + '.pdf',
                                         'label': GeneTree_obj_name+' genome context' + ' PDF'
                                         }])        

        reportClient = KBaseReport(self.callback_url, token=ctx['token'])
        report_info = reportClient.create_extended_report(reportObj)

        # Done
        #
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref']
                  }
        self.log(console, "run_genetree_genome_context() DONE")
        #END run_genetree_genome_context

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_genetree_genome_context return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
