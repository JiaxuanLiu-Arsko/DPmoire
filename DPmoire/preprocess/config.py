import os, yaml, re
from typing import Any
import numpy as np
from ._keys import *
from _pyio import TextIOWrapper
import warnings
import copy
class Config(object):
    config_dict = None

    def __getitem__(self, key:str):
        return self.config_dict[key]
    
    def __setitem__(self, key:str, val):
        self.config_dict[key] = copy.deepcopy(val)
        return key
    
    def items(self):
        return self.config_dict.items()
    
    def keys(self):
        return self.config_dict.keys()
    
    def _as_dict(self):
        return self.config_dict
    
    def __repr__(self):
        return str(dict(self))

    def update(self, config_dict:dict):
        for key in config_dict.keys():
            if key in CONFIG_KEYS:
                self.__setitem__(key, config_dict[key])
            else:
                warnings.warn("Key not found")
        self.check_config()

    def __init__(self, config_dict:dict):
        self.config_dict = copy.deepcopy(DEFAULTS)
        self.update(config_dict)

    @classmethod
    def from_yaml(cls, config_file:str|TextIOWrapper):
        if isinstance(config_file, str):
            with open(config_file, "r") as f:
                config_dict = yaml.load(f, Loader=yaml.FullLoader)
        elif isinstance(config_file, TextIOWrapper):
            config_dict = yaml.load(f, Loader=yaml.FullLoader)
        else:
            raise Exception("config_file must be a string or TextIOWrapper object!")
        c = Config(config_dict=config_dict)
        return c
    
    @classmethod
    def from_dict(cls, config_dict:dict):
        c = Config(config_dict=config_dict)
        return c

    def check_config(self):
        if self.config_dict is None:
            raise Exception("Global config_dict is None!")
        for key, value in self.config_dict.items():
            if value is None:
                raise Exception(f"{key} is not set!")
            if "dir" in key:
                if not os.path.exists(value):
                    raise Exception(f"{key} doesn't exist!")
                else:
                    self.config_dict[key] = os.path.abspath(value)
        if not os.path.exists(f"{self.config_dict['script_dir']}/{self.config_dict['DFT_script']}"):
            raise Exception(f"{self.config_dict['DFT_script']} not in {self.config_dict['script_dir']}")
        if not os.path.exists(f"{self.config_dict['script_dir']}/{self.config_dict['learn_script']}"):
            raise Exception(f"{self.config_dict['learn_script']} not in {self.config_dict['script_dir']}")
        if not os.path.exists(f"{self.config_dict['input_dir']}/top_layer.poscar"):
            raise Exception(f"top_layer.poscar not in {self.config_dict['input_dir']}")
        if not os.path.exists(f"{self.config_dict['input_dir']}/bot_layer.poscar"):
            raise Exception(f"bot_layer.poscar not in {self.config_dict['input_dir']}")
        if not os.path.exists(f"{self.config_dict['input_dir']}/rlx_INCAR"):
            raise Exception(f"min_INCAR.poscar not in {self.config_dict['input_dir']}")
        if not os.path.exists(f"{self.config_dict['input_dir']}/MD_INCAR"):
            raise Exception(f"MD_INCAR.poscar not in {self.config_dict['input_dir']}")
