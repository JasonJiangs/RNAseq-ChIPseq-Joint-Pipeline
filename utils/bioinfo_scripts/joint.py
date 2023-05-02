from utils.parameter.shell_builder import *
from utils.parameter.html_parser import *


class Joint():
    def __init__(self, tools, bash_path, logger):
        self.qc_path = bash_path + '/result/joint/qc/'
        self.tools = tools
        dir_builder(self.qc_path)
        logger.parameter_log('------------------ Joint Load Start ------------------')
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('bash_path: ' + bash_path)
        logger.parameter_log('------------------ Joint Load Finish ------------------')


    def fastqc(self, ctrl, condition, script_only):
        chipseq_ctrl = ctrl.chipseq_controller
        rnaseq_ctrl = ctrl.rnaseq_controller
        joint_ctrl = ctrl.joint_controller

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

        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['fastqc', 'gcc'], file_name)

        chipseq_file_str = loop_concatanator(chipseq_ctrl.config['paired-end'], chipseq_ctrl.config['files'])
        shell_file = open(script_path, 'a')
        shell_file.write('for i in ' + chipseq_file_str + '\n')
        shell_file.write('do\n')
        shell_file.write('fastqc -f ' + file_format + ' -o ' + joint_ctrl.qc_path + condition + '/ ' +
                         chipseq_ctrl.config['dir_path'] + '\"$i\"' + chipseq_ctrl.config['suffix'] + '\n')
        shell_file.write('done\n')
        dir_builder(joint_ctrl.qc_path + condition + '/')

        rnaseq_file_str = loop_concatanator(rnaseq_ctrl.config['paired-end'], rnaseq_ctrl.config['files'])
        shell_file.write('for i in ' + rnaseq_file_str + '\n')
        shell_file.write('do\n')
        shell_file.write('fastqc -f ' + file_format + ' -o ' + joint_ctrl.qc_path + condition + '/ ' +
                         rnaseq_ctrl.config['dir_path'] + '\"$i\"' + rnaseq_ctrl.config['suffix'] + '\n')
        shell_file.write('done\n')
        dir_builder(joint_ctrl.qc_path + condition + '/')

        shell_file.close()

        ctrl.logger.write_log('fastqc script build finished: ' + script_path)

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)

            # summarize html
            ctrl.logger.write_log('Finish FastQC and start summarizing html files for ' + condition)
            dir_path = joint_ctrl.qc_path + condition + '/'
            # loop files in the dir
            for file in os.listdir(dir_path):
                if file.endswith('.html'):
                    summarize_html(ctrl=ctrl,
                                   logger=ctrl.logger,
                                   sample_name=file.split('.')[0],
                                   sample_path=dir_path+file)


    def samtools(self, ctrl, script_only):
        pass
