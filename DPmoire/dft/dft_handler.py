import os, time
import copy
class DFTHandler:
    existing_job = None
    job_list = None
    n_nodes = None
    script_name = None
    def __init__(self, script_name:str, n_nodes:int, existing_job:list=None):
        self.job_list = []
        self.existing_job = copy.deepcopy(existing_job)
        self.n_nodes = n_nodes
        self.script_name = script_name

    def run_calculation(self):
        pass

    def postprocess(self):
        pass

    def check_job(self, job_id:int):
        if len(os.popen(f"squeue | grep {job_id}").readlines())>0:
            return True
        else:
            return False
        
    def get_running_jobs(self, additional_str=""):
        squeue_return = os.popen(f"squeue {additional_str}").readlines()
        job_ids_return = []
        for lines in squeue_return[1:]:
            job_ids_return.append(int(lines.split()[0]))
        return job_ids_return

    def check_job_list(self, job_list:list=None):
        if job_list is None:
            job_list = self.job_list 
        jobs_left = 0
        finished_list = []
        running_jobs = self.get_running_jobs()
        for i, job_id in enumerate(job_list):
            if job_id in running_jobs:
                jobs_left += 1
            else:
                finished_list.append(i)
        finished_list.reverse()
        for idx in finished_list:
            job_list.pop(idx)
        return jobs_left, job_list

    def submit_job(self, work_dir:str):
        time.sleep(5)
        os.chdir(work_dir)
        jobs_left, _ = self.check_job_list(job_list=self.job_list)
        while jobs_left>=self.n_nodes:
            time.sleep(30)
            existing_jobs_left, _ = self.check_job_list(job_list=self.existing_job)
            jobs_left, _ = self.check_job_list(job_list=self.job_list)
            jobs_left += existing_jobs_left
        job_id = int(os.popen(f"sbatch {self.script_name}").readlines()[0].split()[-1])
        self.job_list.append(job_id)

    def wait_until_finished(self):
        jobs_left, _ = self.check_job_list(job_list=self.job_list)
        while jobs_left>0:
            time.sleep(30)
            jobs_left, _ = self.check_job_list(job_list=self.job_list)
        _, existing_jobs = self.check_job_list(job_list=self.existing_job)
        return existing_jobs