
CONFIG_KEYS = ["n_sectors", "work_dir", "input_dir", "POTCAR_dir", "n_nodes", \
         "script_dir", "DFT_script", "username", "learn_script", \
        "do_relaxation", "OUTCAR_collect_freq", "VASP_ML", "init_mlff",\
        "twist_val", "max_val_n", "min_val_n", "d", "sc", "sym_reduce"]
DEFAULTS = {
    "n_sectors":9,
    "work_dir":None,
    "input_dir":None,
    "POTCAR_dir":None,
    "n_nodes":1,
    "script_dir":None,
    "DFT_script":None,
    "username":None,
    "learn_script":None,
    "do_relaxation":True,
    "OUTCAR_collect_freq":7,
    "d":None,
    "VASP_ML":True, 
    "sc":2, 
    "init_mlff":True, 
    "twist_val":False, 
    "max_val_n":0,
    "min_val_n":0, 
    "sym_reduce":False
}