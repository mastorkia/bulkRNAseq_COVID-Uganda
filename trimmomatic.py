import os
import sys
import argparse
import logging

'''Usage:trimmomatic_adaptor.py -f <path containing Raw files > -p <fasta containing adaptor sequences> -c <config file with software information> -a <single or pair end SE/PE>'''
#python trimmomatic_adaptor.py -f fastq.list -p /share/data/RNA_Seq/Matt_Cummings/Working/Raw_fastq/Htseq_adaptors.fa -c config_mapping.yml -a 'SE/PE'
#important adaptor fasta needs to be in the same location as raw files

def msg():
    return ''' This script runs adaptor filtering using Trimmomatic software.
               Exiting.... '''

def parse_arguments():
    '''Adding the command line arguments to be provided while running the script'''
    parser = argparse.ArgumentParser(description='Process command line arguments.', usage=msg())
    parser.add_argument('-f', '--file_list', metavar='FILE', required=True, help='fastq file list')
    parser.add_argument('-c', '--config', metavar='FILE', required=True, help='config file')
    parser.add_argument('-p', '--primers', metavar='FILE', required=True, help='adaptor fasta file')
    parser.add_argument('-a', '--aim', metavar='SE/PE', required=True, help='SE/PE')
    args = parser.parse_args()
    return args

cmd_line_args = parse_arguments()
file_list = cmd_line_args.file_list
print("User provided file list: " + file_list)
primers = cmd_line_args.primers
print("User provided adaptor fasta file: " + primers)
aim = cmd_line_args.aim
print("The aim of the current analysis is: " + aim)


def make_config_hash():

    print("Reading the config file")

    config_dict = {}
    config_file = cmd_line_args.config
    if os.path.isfile(config_file):
        print("Received config file: " + config_file)
    else:
        sys.exit()

    config = open( config_file, 'r')

    for line in config:
        this_line = line.strip()
        this_key = this_line.split(":")[0]
        this_val = this_line.split(":")[1]
        config_dict[this_key] = this_val

    return config_dict


def make_sample_dict():

    sample_dict = {}

    #Iterate through the file list
    if os.path.isfile(file_list):
        print("Received fastq file list: " + file_list)

    fastq_list = open(file_list, 'r')

    for fastq in fastq_list:
        this_fastq = fastq.strip()
        if os.path.isfile(this_fastq):
            print("Found: " + this_fastq)
        else:
            print("this fastq does not exist. Please check the file")

        dir_name = os.path.dirname(this_fastq)
        fastq_file = os.path.basename(this_fastq)
        fastq_name = fastq_file.split(".")[0]

        sample_dict[fastq_name] = dir_name

    return sample_dict


def trimmomatic_adaptor():

    my_config_dict = make_config_hash()
    print("Printing config hash")
    print(my_config_dict)

    sample_dictionary = make_sample_dict()
    print("Printing sample dictionary")
    print(sample_dictionary)

    for fastq_name, dir_name in sample_dictionary.items():

        singularity_image = my_config_dict['singularity_image']
        print("Singularity image: " + singularity_image)
        jar_file = my_config_dict['jar_file']
        print("Trimmomatic version: " + jar_file)
        mem = my_config_dict['local_mem']
        num_cores = my_config_dict['local_cores']
        print("Memory Assigned: " + mem)
        fastq_name2 = fastq_name.split("_")[0]
        job_name = fastq_name2 + "_trimmomatic"
        print(job_name)
        primers_path = os.path.basename(primers)
        output_name = "../Trimmed_reads/output_" + fastq_name2 + "/"


        if aim == 'SE':
            cmd_mkdir = "echo mkdir " + output_name + \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=2g" + " -e " + dir_name + "/" + " -o " +  dir_name + "/" +  " -j y -cwd"

            cmd = "echo singularity exec -B /share/data/ " + singularity_image + " java -jar /opt/trimmomatic/" + jar_file + " SE -threads "+ num_cores + " -trimlog " + fastq_name2 + ".log " + \
            fastq_name + ".fq.gz " + output_name + fastq_name2 + "_filtered.fq.gz " + "ILLUMINACLIP:" + primers_path + ":2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36" + \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=" + mem + " -pe smp " +  num_cores + " -e " + dir_name + "/" + " -o " +  dir_name + "/" +  " -j y -cwd"


        if aim == 'PE':
            cmd_mkdir = "echo mkdir " + output_name + \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=2g" + " -e " + dir_name + "/" + " -o " +  dir_name + "/" +  " -j y -cwd"

            cmd = "echo singularity exec -B /share/data/ " + singularity_image + " java -jar /opt/trimmomatic/" + jar_file + " PE -threads "+ num_cores + " -trimlog ../Trimmed_reads/" + fastq_name2 + ".log " + \
            fastq_name2 + "_R1_001.fastq.gz " + fastq_name2 +"_R2_001.fastq.gz " + "-baseout " + output_name + fastq_name2 + "_filtered.fq.gz "  + "ILLUMINACLIP:" + primers_path + ":2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36" + \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=" + mem + " -pe smp " +  num_cores + " -e " + dir_name + "/" + " -o " +  dir_name + "/" +  " -j y -cwd"


    #generate a SGE command
        print(cmd_mkdir)
        os.chdir(dir_name)
        current_dir = os.getcwd()
        print("\nWorking dir: " + current_dir)
        return_val = os.system(cmd_mkdir)
        
        print(cmd)
        os.chdir(dir_name)
        current_dir = os.getcwd()
        print("\nWorking dir: " + current_dir)
        return_val = os.system(cmd)
        

        if return_val == 0:
            print("trimmomatic job successfully submitted for sample:  " + fastq_name)

        else:
            print("Return value: " + str(return_val))
            print("trimmomatic job cannot be successfully submitted for sample: " + fastq_name)

trimmomatic_adaptor()
