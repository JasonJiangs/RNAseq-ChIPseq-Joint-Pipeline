from utils.parameter.shell_builder import *

class Chipseq():
    def __init__(self, config, tools, base_path, logger):
        self.bowtie2_path = base_path + '/result/chipseq/bowtie2/'
        self.macs2_path = base_path + '/result/chipseq/macs2/'
        self.config = config
        self.tools = tools
        dir_builder(self.bowtie2_path)
        dir_builder(self.macs2_path)
        logger.parameter_log('------------------ Chipseq Load Start ------------------')
        logger.parameter_log('config: ' + str(config))
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('bowtie2_path: ' + str(self.bowtie2_path))
        logger.parameter_log('macs2_path: ' + str(self.macs2_path))
        logger.parameter_log('------------------ Chipseq Load Finish ------------------')

    def bowtie2(self, ctrl, script_only):
        paired_end = ctrl.parameter['config_dict']['datasource']['chip-seq']['paired-end']
        script_path = self.bowtie2_path + 'bowtie2.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['gcc', 'bowtie2'], 'bowtie2_result')

        if paired_end == 'y':
            pass
        else:
            pass

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)


    def macs2(self, ctrl, script_only):
        script_path = self.macs2_path + 'macs2.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['gcc', 'macs2'], 'macs2_result')


        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)


