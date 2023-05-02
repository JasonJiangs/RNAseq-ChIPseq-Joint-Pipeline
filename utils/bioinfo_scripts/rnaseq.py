from utils.parameter.shell_builder import *

class RNASeq():
    def __init__(self, config, tools, base_path, logger):
        self.hisat2_path = base_path + '/result/rnaseq/hisat2/'
        self.stringtie_path = base_path + '/result/rnaseq/stringtie/'
        self.deseq2_path = base_path + '/result/rnaseq/deseq2/'
        self.htseq_path = base_path + '/result/rnaseq/htseq/'
        self.config = config
        self.tools = tools
        dir_builder(self.hisat2_path)
        dir_builder(self.stringtie_path)
        dir_builder(self.deseq2_path)
        dir_builder(self.htseq_path)
        logger.parameter_log('------------------ RNASeq Load Start ------------------')
        logger.parameter_log('config: ' + str(config))
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('hisat2_path: ' + str(self.hisat2_path))
        logger.parameter_log('stringtie_path: ' + str(self.stringtie_path))
        logger.parameter_log('deseq2_path: ' + str(self.deseq2_path))
        logger.parameter_log('htseq_path: ' + str(self.htseq_path))
        logger.parameter_log('------------------ RNASeq Load Finish ------------------')


    def hisat2(self, ctrl, script_only):
        paired_end = ctrl.parameter['config_dict']['datasource']['rna-seq']['paired-end']
        script_path = self.hisat2_path + 'hisat2.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path,
                              ['gcc', 'hisat2'], 'hisat2_result' )

        if paired_end == 'y':
            pass
        else:
            pass

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)
        else:
            ctrl.logger.write_log('Only script is generated: ' + script_path)


    def stringtie(self, ctrl, script_only):
        script_path = self.stringtie_path + 'stringtie.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path,
                              ['gcc', 'stringtie'], 'stringtie_result' )

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)

    def deseq2(self, ctrl, script_only):
        pass

    def htseq(self, ctrl, script_only):
        pass
