import os
windows = ['04enhancer_peak_distribution01','04enhancer_peak_distribution02', '04enhancer_peak_distribution03', '04enhancer_peak_distribution04', '04enhancer_peak_distribution05']
for i in windows:
    os.system("nohup python -u "+i+".py >> "+i+".log 2>&1 &")