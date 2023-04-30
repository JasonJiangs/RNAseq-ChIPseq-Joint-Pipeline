from utils.parameter.shell_builder import general_shell_builder, shell_runner, dir_builder


class Joint():
    def __init__(self, tools):
        self.qc_path = None
        self.tools = tools


def fastqc(ctrl, condition):
    chipseq_ctrl = ctrl.chipseq_controller
    rnaseq_ctrl = ctrl.rnaseq_controller
    joint_ctrl = ctrl.joint_controller
    # assign qc_pth if qc_path is not exist
    if joint_ctrl.qc_path is None:
        joint_ctrl.qc_path = ctrl.base_path + '/result/joint/qc/'

    if condition == 'before':
        file_format = 'fastq'
        file_name = 'fastqc_before'
        script_path = ctrl.base_path + '/result/shell_script/fastqc_before.sh'
    elif condition == 'after':
        file_format = 'bam'
        script_path = ctrl.base_path + '/result/shell_script/fastqc_after.sh'
        file_name = 'fastqc_after'
    else:
        ctrl.logger.write_log(ctrl, 'Error: qc condition is not correct.')
        exit(-1)

    module_list = ['fastqc', 'gcc']
    log_path = ctrl.base_path + '/result/shell_log/'

    general_shell_builder(ctrl.slurm, script_path, log_path, module_list, file_name)

    chipseq_file_str = ' '
    if chipseq_ctrl.source_path['paired-end'] == 'y':
        nw_list = []
        for i in range(len(chipseq_ctrl.source_path['files'])):
            nw_list.append(chipseq_ctrl.source_path['files'][i]+'_1')
            nw_list.append(chipseq_ctrl.source_path['files'][i]+'_2')
        chipseq_file_str = chipseq_file_str.join(nw_list)
    else:
        chipseq_file_str = chipseq_file_str.join(chipseq_ctrl.source_path['files'])
    shell_file = open(script_path, 'a')
    shell_file.write('for i in ' + chipseq_file_str + '\n')
    shell_file.write('do\n')
    shell_file.write('fastqc -f ' + file_format + ' -o ' + joint_ctrl.qc_path + condition + '/ ' +
                     chipseq_ctrl.source_path['dir_path'] + '\"$i\"' + chipseq_ctrl.source_path['suffix'] + '\n')
    shell_file.write('done\n')
    dir_builder(joint_ctrl.qc_path + condition + '/')

    rnaseq_file_str = ' '
    if rnaseq_ctrl.source_path['paired-end'] == 'y':
        nw_list = []
        for i in range(len(rnaseq_ctrl.source_path['files'])):
            nw_list.append(rnaseq_ctrl.source_path['files'][i]+'_1')
            nw_list.append(rnaseq_ctrl.source_path['files'][i]+'_2')
        rnaseq_file_str = rnaseq_file_str.join(nw_list)
    else:
        rnaseq_file_str = rnaseq_file_str.join(rnaseq_ctrl.source_path['files'])
    shell_file.write('for i in ' + rnaseq_file_str + '\n')
    shell_file.write('do\n')
    shell_file.write('fastqc -f ' + file_format + ' -o ' + joint_ctrl.qc_path + condition + '/ ' +
                        rnaseq_ctrl.source_path['dir_path'] + '\"$i\"' + rnaseq_ctrl.source_path['suffix'] + '\n')
    shell_file.write('done\n')
    dir_builder(joint_ctrl.qc_path + condition + '/')

    shell_file.close()

    linux_command = 'sbatch ' + script_path
    shell_runner(linux_command, ctrl)


def samtools(ctrl):

    pass

