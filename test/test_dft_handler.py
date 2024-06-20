from DPmoire.dft.dft_handler import DFTHandler
import os
job_existing = [int(strs.split()[0]) for strs in os.popen(f"squeue").readlines()[1:]]
dft_handler = DFTHandler("DFT.sh", 81, [job_existing[0]])
print(job_existing[1])
print(dft_handler.check_job_list([job_existing[1], 1749999]))
dft_handler.submit_job(".")
dft_handler.wait_until_finished()
dft_handler.check_job_list(job_existing[0])