# list_samples.py example
import os
import sys

directory = sys.argv[1]
samples_file = sys.argv[2]

samples = set()
for filename in os.listdir(directory):
    if filename.endswith(".fq.gz") or filename.endswith(".fastq.gz"):
        sample_name = "_".join(filename.split("_")[:-1])
        samples.add(sample_name)

with open(samples_file, 'w') as f:
    for sample in sorted(samples):
        f.write(sample + '\n')