import os
import re
import Calculate_probe_coordinates as Calc_probe_coord
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import subprocess
from subprocess import Popen


# Paths
# Obtain the path where the script is being run
path = '/mnt/fastdata/mbp15ja/capturec-pilot/capture_c_perl'

# Probes coordinates
path_probes = '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw2'

# Restriction site fragments
path_fragments = os.path.join(path_probes, "digest.dir")

# Flashed reads
flash_dir = '/fastdata/mbp15ja/capturec-pilot/capture_c_fish_work_desktop/flashed.dir/'




# Creates a script with the commands specified to be run. Each script can run in parallel.
def write_parallel_processing_script(script_file, command):
    
    with IOTools.openFile(script_file, "w") as outf:

        out_command = ""
        out_command += "#! /bin/bash" + "\n"
        out_command += "# -V" + "\n"
        out_command += "source ~/.bashrc" + "\n"
        out_command += command + "\n"
    
        outf.write(out_command)
    
    
    
    return script_file


if __name__ == "__main__":
    
    # Obtain the path where the script is being run
    path = '/mnt/fastdata/mbp15ja/capturec-pilot/capture_c_perl'
#     
#     # Generate the coordinates as specified in the Perl Capture C Manual 
#     
#     # Using pregenerated digest fragments from Capture C Sudlab run
#     # To be implemented: From pipeline_caputrec: splitDigest, mergeDigest, digest2fragments
#     outfile_probe_fragments = os.path.join(path, "probe_fragments.bed")
#     outfile_lookup_out = os.path.join(path, "lookup_out")
#     
#     infile_probes = os.path.join(path_probes, "probes.bed.gz")
#     infile_fragments = os.path.join(path_fragments, "fragments.bed.gz")
#     
#     Calc_probe_coord.getProbeFragments(infile_probes,
#                                infile_fragments, 
#                                outfile_probe_fragments,
#                                outfile_lookup_out)
#     
#     
#     
#     
#     # Deduplicate restriction fragments
    # each DpnII fragment only ONCE in the oligo coordinate file
    outfile_probe_fragments_dedup = os.path.join(path, "dedup_probe_fragments.bed")
    
#     statement = ' '.join(["sort", "--field-separator=$'\t'", "--key=1,1", "--key=2,2", "--key=3,3", "--key=5,5", "--key=6,6", ("--unique < "+outfile_probe_fragments), (" > "+outfile_probe_fragments_dedup)])
#     
#     process = subprocess.Popen(statement, shell=True, stdout=subprocess.PIPE)
#     process.wait()
#     
#     if(process.returncode != 0 ):
#         print "Error in the deduplication of files"    
#     
#     
    # Output oligo file with 9 column line format
    outfile_oligo_coord_dedup = os.path.join(path, "oligo_coord_dedup.txt")
    
    (probe_collisions, exclus_collisions) = Calc_probe_coord.formatProbeFragments(outfile_probe_fragments_dedup, outfile_oligo_coord_dedup)
    
    print(probe_collisions)
    
    print(exclus_collisions)
    
#     # Using pregenerated flashed reads: flashReads() method (before deduplicating) from Capture_C_Fish
#     # To be implementec here
#     
#     
#     # Get all files in directory
#     flashed_files = os.listdir(flash_dir)
#     
#     # Execution statement
#     combine_statement = ""
#     
#     # Alignment
#     alignment_statement = ""
#     
#     # Alignment script register
#     script_register = []
#     
#     
#     
#     for flashed_file in flashed_files:
#         match = re.match("(.+).extendedFrags.fastq.gz", flashed_file)
#         if match:
#             id_file = match.group(1)
#             
#             # Create file with it's full path to collect flashed file
#             infile = os.path.join(flash_dir, id_file)
#             
#             # Create the output file in the script directory
#             outfile = os.path.join(path, id_file)
#             
#             combine_statement = ' '.join(["zcat", (infile+".notCombined.fastq.gz"), (infile+".extendedFrags.fastq.gz >"), (outfile+".fastq;  ")])
#             
#             alignment_statement = (' '.join(["bowtie -p 1 -m 2 --best --strata --sam --chunkmb 256 --sam /shared/sudlab1/General/mirror/genomes/bowtie/hg19", (outfile+".fastq"), (outfile+".sam"), "2>", (outfile+".log") ]))
#             
#             parallel_processing_statement = write_parallel_processing_script((outfile+".sh"), (combine_statement+alignment_statement))
#             
#             script_register.append(parallel_processing_statement)
#             
# 
#     
#     
#     # Parallel execution
#     processes = []    
#     
#     # Execution command
#     exec_command = ""
#     
#     for script in script_register:
#         
#         exec_command = "qsub -cwd -q openmp.q -pe openmp 1 -l rmem=4G -l mem=4G "+ script
#         
#         print(exec_command)
#          
#         process = Popen(exec_command, shell=True)
#           
#         processes.append(process)
#         
#     # wait for completion
#     for p in processes: p.wait()
    
    
        
        
           
        