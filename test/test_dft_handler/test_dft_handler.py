from DPmoire.dft.dft_handler import DFTHandler, check_job_list_status
import os
import time
job_existing = [strs.split()[0] for strs in os.popen(f"squeue").readlines()[1:3]]
job_existing.append("100")
dft_handler = DFTHandler("DFT.sh", 81, job_existing)
print(job_existing)
print(check_job_list_status(job_existing))
print(check_job_list_status(["100"]))
dft_handler.submit_job(".")
print(dft_handler.job_list)
time.sleep(5)
print(dft_handler.get_running_jobs())
dft_handler.wait_until_finished()
dft_handler.get_running_jobs([job_existing[0]])