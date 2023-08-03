
## This bash script is being updated by me. 

# The following two folders are necessary to be updated. 


# /home/madhu/datasets/scripts
# original_folder="/mnt/labshare/share/public_data/DNAseq/forkseq/"

original_folder="/home/madhu/datasets/scripts"

# The folder where the bash script is saved
# basecalled_folder="/mnt/labshare/share/public_data/DNAseq/forkseq/guppy_v5.0.11/"

basecalled_folder="/home/madhu/datasets/download/GSM3528751/data/extracted_data/fast5/extracted_DATA"

Dataset_name=$1
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/labshare/share/tools/guppy/download2023Jan17/v6.1.3/opt/ont/guppy/lib/:/mnt/labshare/share/tools/guppy/v5.0.11/opt/ont/guppy/lib; /mnt/labshare/share/tools/guppy/v5.0.11/opt/ont/guppy/bin/guppy_basecaller -i ${original_folder}/${Dataset_name} -s ${basecalled_folder}/${Dataset_name} -c dna_r9.4.1_450bps_hac.cfg \
                 --recursive  --num_callers 2 --gpu_runners_per_device 1 --cpu_threads_per_caller 4  --device cuda:1 --fast5_out -q 4000
