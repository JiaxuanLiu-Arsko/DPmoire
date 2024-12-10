import argparse, os, yaml, time

from ..preprocess import Config, EnvironmentHandler
from ..data import Dataset
from ..train import NequIPTrainer
from ..dft import RelaxationHandler, MDHandler, ValidationHandler
from copy import deepcopy
import numpy as np
def main(args=None):

    parser = argparse.ArgumentParser(
        description="A module to automatically train a Machine Learning forcefield for a Twisted Bilayer Material"
    )
    parser.add_argument(
        "config", help="YAML file configuring the module"
    )

    parser.add_argument(
        "--mode", default="all", type=str, help="mode: all, run, train, or build_val."
    )

    args = parser.parse_args(args=args)

    if not os.path.isfile(args.config):
        print("YAML not found!")
    config = Config.from_yaml(args.config)
    start_time = time.time()
    dataset = Dataset()
    val_dataset = Dataset()
    job_list = []
    env_handler = EnvironmentHandler(config=config)
    if config["sym_reduce"]:
        stackings = env_handler.find_sym_reduced_stackings()
    else:
        stackings = np.array([[i, j] for i in range(config["n_sectors"]) for j in range(config["n_sectors"])])
    md_handler = MDHandler(config=config, stackings=stackings)
    if args.mode == "build_val":
        if not config["twist_val"]:
            raise ValueError("twist_val in config.yaml is not True!")
        work_dir = config["work_dir"]
        angles = env_handler.gen_val_environment( config["min_val_n"], config['max_val_n'], f"{work_dir}/validation")
        val_handler = ValidationHandler(config=config, angles=angles, val_dir=f"{work_dir}/validation")
        val_handler.run_calculation()
        _, job_list = val_handler.get_running_jobs()
        val_handler.wait_until_finished()
        val_dataset = val_handler.postprocess()
        val_dataset.save_extxyz(f"{work_dir}/valid.extxyz")
    if args.mode == "all" or args.mode == "run":
        work_dir = config["work_dir"]
        if config["twist_val"]:
            angles = env_handler.gen_val_environment( config["min_val_n"], config['max_val_n'], f"{work_dir}/validation")
            val_handler = ValidationHandler(config=config, angles=angles, val_dir=f"{work_dir}/validation")
            val_handler.run_calculation()
            _, job_list = val_handler.get_running_jobs()
        if config["init_mlff"]:
            env_handler.gen_init_environment(f"{work_dir}/init_mlff", 0)
            os.system(f"cp {config['script_dir']}/{config['DFT_script']} {work_dir}/init_mlff")
            md_handler = MDHandler(config=config, existing_job=job_list, stackings=stackings)
            md_handler.submit_job(work_dir=f"{work_dir}/init_mlff")
            md_handler.wait_until_finished()
            env_handler.gen_init_environment(f"{work_dir}/init_mlff", 1)
            md_handler.submit_job(work_dir=f"{work_dir}/init_mlff")
            job_list = md_handler.wait_until_finished()
            os.system(f"cp {work_dir}/init_mlff/ML_ABN {config['input_dir']}/ML_AB ")
            os.system(f"cp {work_dir}/init_mlff/ML_FFN {config['input_dir']}/ML_FF ")
        if config["do_relaxation"]:
            env_handler.gen_POSCAR(config["work_dir"], stackings=stackings)
            env_handler.gen_environment('rlx_INCAR', config["work_dir"], stackings=stackings)
            rlx_handler = RelaxationHandler(config=config, existing_job=job_list, stackings=stackings)
            rlx_handler.run_calculation()
            job_list = rlx_handler.wait_until_finished()
            rlx_dataset = rlx_handler.postprocess()
            dataset.load_dataset_class(rlx_dataset)
            print(f"Relaxation finished. {time.time() - start_time}s consumed.")
            start_time = time.time()
        else:
            if os.path.exists(f"{work_dir}/rlx_data.extxyz"):
                dataset.load_dataset_extxyz(f"{work_dir}/rlx_data.extxyz")

        md_handler = MDHandler(config=config, existing_job=job_list, stackings=stackings)
        env_handler.gen_environment('MD_INCAR', config["work_dir"], stackings=stackings)
        md_handler.run_calculation()
        md_dataset = md_handler.postprocess()
        dataset.load_dataset_class(md_dataset)
        if config["twist_val"]:
            val_handler.wait_until_finished()
            val_dataset = val_handler.postprocess()
            val_dataset.save_extxyz(f"{work_dir}/valid.extxyz")
        print(f"MD finished. {time.time() - start_time}s consumed.")
        start_time = time.time()

    if args.mode == "all" or args.mode == "train":
        if args.mode == "train":
            if os.path.exists(f"{config['work_dir']}/rlx_data.extxyz"):
                dataset.load_dataset_extxyz(f"{config['work_dir']}/rlx_data.extxyz")
                print("Relaxation data loaded.")
            else:
                print("Warninng! No rlx_data.extxyz found! No relaxation data loaded!")
            if os.path.exists(f"{config['work_dir']}/MD_data.extxyz"):
                dataset.load_dataset_extxyz(f"{config['work_dir']}/MD_data.extxyz")
                print("MD data loaded.")
            else:
                md_dataset = md_handler.make_dataset()
                dataset.load_dataset_class(md_dataset)
                print("Warninng! No MD_data.extxyz found, load MD data from ML_ABN files!")
            if os.path.exists(f"{config['work_dir']}/valid.extxyz"):
                val_dataset.load_dataset_extxyz(f"{config['work_dir']}/valid.extxyz")
            else:
                print("Warning! No valid.extxyz found!")
        trainer = NequIPTrainer(config=config, dataset=dataset, val_dataset=val_dataset)
        trainer.preprocess(RCUT=env_handler.find_RCUT())
        trainer.submit_training()
        trainer.postprocess()
        print(f"Training finished. Time cost = {time.time()-start_time} secs.")

if __name__ == "__main__":
    main()