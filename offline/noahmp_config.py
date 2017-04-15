# -*- coding: utf8 -*-
import sys
import os.path
import datetime
import f90nml


NML_FIELDS = ['static_parameter_file',
              'initialization_file',
              'restart_file',

              'input_directory',
              'input_frequency',
              'output_directory',
              'output_frequency',
              'restart_directory',
              'restart_frequency',

              'start_year',
              'start_month',
              'start_day',
              'start_hour',
              'start_minute',
              'start_second',

              'end_year',
              'end_month',
              'end_day',
              'end_hour',
              'end_minute',
              'end_second',

              'interval_seconds',

              'opt_veg',
              'opt_run',
              'opt_btr',
              'opt_rad',
              'opt_tub',
              'opt_can',
              'opt_inf',
              'opt_snf',
              'opt_tbot']



class Config(object):
    def __init__(self, cfgfile):
        self.indir = '.'
        self.infreq = None
        self.outdir = '.'
        self.outfreq = None
        self.resdir = '.'
        self.resfreq = None

        self.constfile = 'domain.nc'
        self.initfile = 'init.nc'

        self.datetimebeg = None
        self.datetimeend = None
        self.timestep = 0

        self.parse_cfg(cfgfile)

    def parse_cfg(self, cfgfile):
        if not os.path.isfile(cfgfile):
            print('ERR: Unable to find configuration file', cfgfile)
            sys.exit(1)
        nml = f90nml.read(cfgfile)
        cfg = nml['NOAHMP_OFFLINE']
        # Mandatory fields
        for var in NML_FIELDS:
            if var not in cfg:
                print('ERR: Unable to find {:s} in configuration file {:s}'.format(var, cfgfile))
                sys.exit(1)
        # Initialization
        self.constfile = cfg['static_parameter_file']
        self.initfile = cfg['initialization_file']
        self.resfile = cfg['restart_file']
        # Input & Output
        self.indir = cfg['input_directory']
        self.infreq = cfg['input_frequency']
        self.outdir = cfg['output_directory']
        self.outfreq = cfg['output_frequency']
        self.resdir = cfg['restart_directory']
        self.resfreq = cfg['restart_frequency']
        # Physics
        self.opt_veg = cfg['opt_veg']
        self.opt_run = cfg['opt_run']
        self.opt_btr = cfg['opt_btr']
        self.opt_rad = cfg['opt_rad']
        self.opt_inf = cfg['opt_inf']
        self.opt_snf = cfg['opt_snf']
        self.opt_tub = cfg['opt_tub']
        self.opt_can = cfg['opt_can']
        self.opt_tbot = cfg['opt_tbot']
        # Model Temporal settings
        self.timestep = datetime.timedelta(seconds=cfg['interval_seconds'])
        self.begdatetime = datetime.datetime(cfg['start_year'], cfg['start_month'], cfg['start_day'],
                                             cfg['start_hour'], cfg['start_minute'], cfg['start_second'])
        self.enddatetime = datetime.datetime(cfg['end_year'], cfg['end_month'], cfg['end_day'],
                                             cfg['end_hour'], cfg['end_minute'], cfg['end_second'])
        pass
