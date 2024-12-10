import os, time
import copy
import re
import numpy as np



def check_job_list_status(job_list:list=None):
    if len(job_list)==0:
        return []
    job_ids_str = ",".join(map(str, job_list))
    command = f'sacct -j {job_ids_str} --format JobID,State'
    with os.popen(command) as process: 
        response = process.readlines()[2:]
    results = {}
    for status in response:
        words = status.split()
        if words[0] in job_list:
            if len(re.findall(r'(PENDING)|(RUNNING)', words[1]))>0:
                results[words[0]] = 0  # running or pending
            elif len(re.findall('COMPLETED', words[1])):
                results[words[0]] = 1  # normally completed
            elif len(re.findall('FAILED', words[1])):
                results[words[0]] = 2  # job failed
            else:
                results[words[0]] = 3  # other status
    if len(results)<len(job_list):
        print("-----------------------------")
        print(job_list)
        print(results)
        print(command)
        print("-----------------------------")
    return results


def check_job_status(job_id:str):
    with os.popen(f"sacct -j {job_id} --format State") as process:
        status = process.readlines()[2]
    if len(re.findall(r'(PENDING)|(RUNNING)', status))>0:
        return 0 # running or pending
    elif len(re.findall('COMPLETED', status)):
        return 1 # normally completed
    elif len(re.findall('FAILED', status)):
        return 2 # job failed
    else:
        return 3 # other status

class DFTHandler:
    existing_job = None
    job_list = None
    n_nodes = None
    script_name = None
    job_work_dir = None
    auto_resub = None
    def __init__(self, script_name:str, n_nodes:int, existing_job:list=None, auto_resub:bool = False):
        self.job_list = []
        self.existing_job = copy.deepcopy(existing_job)
        self.n_nodes = n_nodes
        self.script_name = script_name
        self.job_work_dir = {}
        self.auto_resub = auto_resub

    def run_calculation(self):
        pass

    def postprocess(self):
        pass
        
    
    def get_running_jobs(self, job_list:list=None):
        if job_list is None:
            job_list = self.job_list 
        results = check_job_list_status(job_list)
        jobs_left = 0
        running_jobs = []
        for i, job_id in enumerate(job_list):
            if results[job_id] == 0:
                jobs_left += 1
                running_jobs.append(job_id)
        return jobs_left, running_jobs

    def submit_job(self, work_dir:str):
        os.chdir(work_dir)
        jobs_left, _ = self.get_running_jobs(job_list=self.job_list)
        existing_jobs_left, _ = self.get_running_jobs(job_list=self.existing_job)
        jobs_left += existing_jobs_left
        while jobs_left>=self.n_nodes:
            time.sleep(30)
            existing_jobs_left, _ = self.get_running_jobs(job_list=self.existing_job)
            jobs_left, _ = self.get_running_jobs(job_list=self.job_list)
            jobs_left += existing_jobs_left
        with os.popen(f"sbatch {self.script_name}") as process:
            job_id = process.readlines()[0].split()[-1]
        self.job_list.append(job_id)
        self.job_work_dir[job_id] = work_dir
        time.sleep(5)

    def wait_until_finished(self):
        while len(self.job_list)>0:
            finished_jobs = []
            failed_jobs = []
            results = check_job_list_status(job_list=self.job_list)
            for idx, job_id in enumerate(self.job_list):
                if results[job_id] == 0:
                    continue
                elif results[job_id] == 1:
                    finished_jobs.append(job_id)
                else:
                    failed_jobs.append(job_id)
            
            self.job_list = [v for v in self.job_list if v not in finished_jobs]
            self.job_list = [v for v in self.job_list if v not in failed_jobs]
            if self.auto_resub:
                for idx, job_id in enumerate(failed_jobs):
                    self.submit_job(self.job_work_dir[job_id])
                    print(f"Slurm job {job_id} failed/canceled. DPmoire will resubmit a job in {self.job_work_dir[job_id]}.")
                    if job_id in self.job_work_dir:
                        del self.job_work_dir[job_id]
            else:
                for idx, job_id in enumerate(failed_jobs):
                    self.submit_job(self.job_work_dir[job_id])
                    print(f"Slurm job {job_id} in {self.job_work_dir[job_id]} failed/canceled. `auto_resub` tag was set to False, there will be no resubmission.")
            time.sleep(30)
        _, existing_jobs = self.get_running_jobs(job_list=self.existing_job)
        return existing_jobs