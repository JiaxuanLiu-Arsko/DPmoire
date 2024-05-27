import os, time

class DFTHandler:
    job_list = None
    n_nodes = None
    script_name = None
    def __init__(self, script_name:str, n_nodes:int, existing_job:list=None):
        if existing_job is None:
            self.job_list = []
        else:
            self.job_list = existing_job
        pass
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
    
    def check_job_list(self, job_list:list):
        jobs_left = 0
        finished_list = []
        for i, job_id in enumerate(job_list):
            if self.check_job(job_id=job_id):
                jobs_left += 1
            else:
                finished_list.append(i)
        finished_list.reverse()
        for idx in finished_list:
            job_list.pop(idx)
        return jobs_left

    def submit_job(self, work_dir:str):
        time.sleep(5)
        os.chdir(work_dir)
        jobs_left = self.check_job_list(job_list=self.job_list)
        while jobs_left>self.n_nodes:
            time.sleep(30)
            jobs_left = self.check_job_list(job_list=self.job_list)
        job_id = int(os.popen(f"sbatch {self.script_name}").readlines()[0].split()[-1])
        self.job_list.append(job_id)

    def wait_until_finished(self):
        jobs_left = self.check_job_list(job_list=self.job_list)
        while jobs_left>0:
            time.sleep(30)
            jobs_left = self.check_job_list(job_list=self.job_list)