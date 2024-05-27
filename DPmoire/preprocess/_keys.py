
CONFIG_KEYS = ["n_sectors", "work_dir", "input_dir", "POTCAR_dir", "n_nodes", \
         "script_dir", "DFT_script", "username", "learn_script", "test_dir", \
        "do_relaxation", "OUTCAR_collect_freq", "VASP_ML", "init_mlff", 
        "d", "sc"]
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
    "test_dir":None,
    "do_relaxation":True,
    "OUTCAR_collect_freq":7,
    "d":None,
    "VASP_ML":True, 
    "sc":2, 
    "init_mlff":True
}