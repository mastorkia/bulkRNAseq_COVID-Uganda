import os
import sys
import argparse
import logging

'''Usage:STAR_mapping.py -f <path containing samples > -g <whole path containing genome> -c <config file with software information> -a <single or pair end SE/PE>'''
#python STAR_mapping.py -c config_mapping.yml -a PE -f sample.list -g /share/data/Illumina_Databases/human/star/genocode_GRCh38-v41/

#(needs to be the output of trimmomatic)
#ls -d /share/data/RNA_Seq/Matt_Cummings/Working/Trimmed_reads/output_* > sample.list
#cat sample.lst
#/share/data/RNA_Seq/ETT/star_Aug29_2022/test_raw_data/output_MP22-PRE-Control-Null
#/share/data/RNA_Seq/ETT/star_Aug29_2022/test_raw_data/output_MP30-PRE-Case-Null
#/share/data/RNA_Seq/ETT/star_Aug29_2022/test_raw_data/output_MP45-PRE-Case-Null

def msg():
    return ''' This script runs star mapping program to generate mapping alignemnts.
               Exiting.... '''

def parse_arguments():
    '''Adding the command line arguments to be provided while running the script'''
    parser = argparse.ArgumentParser(description='Process command line arguments.', usage=msg())
    parser.add_argument('-f', '--file_list', metavar='FILE', required=True, help='sample file list')
    parser.add_argument('-c', '--config', metavar='FILE', required=True, help='config file')
    parser.add_argument('-g', '--genome_path', metavar='DIR', required=True, help='genome folder whole path')
    parser.add_argument('-a', '--aim', metavar='SE/PE', required=True, help='SE/PE')
    args = parser.parse_args()
    return args

cmd_line_args = parse_arguments()
file_list = cmd_line_args.file_list
print("User provided file : " + file_list)
genome_path = cmd_line_args.genome_path
print("User provided genome path : " + genome_path)
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
        print("Received sample file list: " + file_list)

    sample_list = open(file_list, 'r')

    for sample in sample_list:
        this_sample = sample.strip()
        if os.path.isdir(this_sample):
            print("Found: " + this_sample)
        else:
            print("this directory does not exist. Please check the file")

        dir_name = os.path.abspath(this_sample)
        sample_file = os.listdir(dir_name)[0]
        sample_dict[sample_file] = dir_name

    return sample_dict


def STAR_mapping():

    my_config_dict = make_config_hash()
    print("Printing config hash")
    print(my_config_dict)

    all_path_genome = os.path.abspath(genome_path)
    genome_dir = os.path.dirname(all_path_genome)
    
    genome_annot = all_path_genome.split(".")[0]

    sample_dictionary = make_sample_dict()
    print("Printing sample dictionary")
    print(sample_dictionary)

    for sample_file, dir_name in sample_dictionary.items():

        STAR_container = my_config_dict['star_container_loc']
        print("Software location: " + STAR_container)
        mem = my_config_dict['local_mem']
        num_cores = my_config_dict['local_cores']
        print("Memory Assigned: " + mem)
        MultimappingMax = my_config_dict['max_multiMapping']
        SAM_multiMax = my_config_dict['SAM_multiMax']
        sample_file2 = sample_file.split(".")[0]
        sample_file3 = sample_file2.split("_1P")[0]
        sample_file4 = sample_file3.split("_")[0]
        job_name = sample_file3 + "_STAR"
        print(job_name)

        if aim == 'SE':
            cmd = "echo singularity exec -B /share/data " + STAR_container + " STAR --runThreadN " + num_cores + " --genomeDir " + genome_annot + \
            " --readFilesIn " + sample_file2 + ".fq.gz " +" --readFilesCommand zcat --sjdbGTFfile " + genome_annot + "/gencode.v41.annotation.gtf --limitBAMsortRAM 1076403423 "  +\
            " --outFileNamePrefix aligned_" + sample_file2 + " --outSAMmultNmax " + SAM_multiMax + " --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax " + MultimappingMax +  \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=" + mem + " -pe smp " +  num_cores + " -e " + dir_name + "/" + " -o " +  dir_name + "/" +  " -j y -cwd"


        if aim == 'PE':

            cmd = "echo singularity exec -B /share/data " + STAR_container + " STAR --runThreadN " + num_cores + " --genomeDir " + genome_annot + \
            " --readFilesIn " + sample_file3 + "_1P.fq.gz " + sample_file3 + "_2P.fq.gz" +" --readFilesCommand zcat --sjdbGTFfile " + genome_annot + "/gencode.v41.annotation.gtf "  +\
            " --outFileNamePrefix aligned_" + sample_file4 + " --outSAMmultNmax " + SAM_multiMax + " --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax " + MultimappingMax +  \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=" + mem + " -pe smp " +  num_cores + " -e " + dir_name + "/" + " -o " +  dir_name + "/" +  " -j y -cwd"


    #generate a SGE command
        print(cmd)

        os.chdir(dir_name)
        current_dir = os.getcwd()
        print("\nWorking dir: " + current_dir)
        return_val = os.system(cmd)

        if return_val == 0:
            print("STAR_sorted job successfully submitted for sample:  " + sample_file3)

        else:
            print("Return value: " + str(return_val))
            print("STAR_sorted job cannot be successfully submitted for sample: " + sample_file3)

STAR_mapping()
