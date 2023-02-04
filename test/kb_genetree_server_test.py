# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from kb_genetree.kb_genetreeImpl import kb_genetree
from kb_genetree.kb_genetreeServer import MethodContext
from kb_genetree.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class kb_genetreeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_genetree'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_genetree',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_genetree(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')


    ##############
    # UNIT TESTS #
    ##############


    #### Gene Tree Genome Context
    ##
    # HIDE @unittest.skip("skipped test_01_run_genetree_genome_context()")  # uncomment to skip
    def test_01_run_genetree_genome_context(self):
        method = '01_run_genetree_genome_context'

        print ("\n\nRUNNING: test_"+method)
        print ("==================================================\n\n")

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        
        debug_workspace_name = 'dylan:1471889662419'
        debug_workspace_id = 16750
        debug_genetree_ref = '16750/78/1'
        
        ret = self.serviceImpl.run_genetree_genome_context (self.ctx,
                     #{'workspace_name': self.wsName,
                     {'workspace_name': debug_workspace_name,
                      'input_genetree_ref': debug_genetree_ref,
                      'genome_disp_name_config': "obj_name",
                      'slice_width': 20.0,
                      'max_rows': 100,
                      'prevalence_color_threshold': 25.0
                     }
                    )
        print ("DONE")
