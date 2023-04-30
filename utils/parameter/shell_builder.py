import os


def general_shell_builder(slurm, execution_path, log_path, module_list, file_name):
    # if exist, delete the old shell script
    if os.path.exists(execution_path):
        os.remove(execution_path)

    # create execution path
    if not os.path.exists(os.path.dirname(execution_path)):
        os.makedirs(os.path.dirname(execution_path))

    # build general format of shell script for the path
    shell_file = open(execution_path, 'w')
    shell_file.write("#!/bin/bash\n")
    shell_file.write("#SBATCH --partition=" + slurm['partition'] + "\n")
    shell_file.write("#SBATCH --time=" + slurm['time'] + "\n")
    shell_file.write("#SBATCH --mem=" + slurm['mem'] + "\n")
    shell_file.write("#SBATCH --cpus-per-task=" + str(slurm['cpus-per-task']) + "\n")
    shell_file.write("#SBATCH -A " + slurm['A'] + "\n")
    shell_file.write("#SBATCH --job-name=" + file_name + "\n")
    shell_file.write("#SBATCH --output=\"" + log_path + file_name + ".out\"\n\n")

    # build slurm log output path if not exist
    dir_builder(log_path)

    # load module
    for module in module_list:
        shell_file.write("module load " + module + "\n")
    shell_file.write("\n")

    shell_file.close()


def path_builder(dir_path, file_name, suffix):
    return dir_path + file_name + suffix


def slurm_script():
    return


def shell_runner(execution_path, ctrl):
    ctrl.logger.write_log(ctrl, "Start executing shell script: " + execution_path)
    signal = os.system(execution_path)
    if signal != 0:
        ctrl.logger.write_log(ctrl, "Error: shell script execution failed.")
    return

def dir_builder(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return

def paired_end_processor(ctrl, paired_end):
    if paired_end == 'y':
        return '-1 ' + ctrl.parameters['rnaseq_source']['fastq1'] + ' -2 ' + ctrl.parameters['rnaseq_source']['fastq2']
    else:
        return '-U ' + ctrl.parameters['rnaseq_source']['fastq1']