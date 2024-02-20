import os
import sys

directory = sys.argv[1]
samples = set()

for filename in os.listdir(directory):
    if filename.endswith(".fq.gz") or filename.endswith(".fastq.gz"):
        sample_name = "_".join(filename.split("_")[:-1])
        samples.add(sample_name)

for sample in sorted(samples):
    print(sample)

