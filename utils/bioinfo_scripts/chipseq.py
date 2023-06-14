from utils.parameter.shell_builder import *


class Chipseq():
    def __init__(self, config, tools, base_path, logger):
        self.bowtie2_path = base_path + '/result/chipseq/bowtie2/'
        self.macs2_path = base_path + '/result/chipseq/macs2/'
        self.config = config
        self.tools = tools
        self.logger = logger
        dir_builder(self.bowtie2_path)
        dir_builder(self.macs2_path)
        logger.parameter_log('------------------ Chipseq Load Start ------------------')
        logger.parameter_log('config: ' + str(config))
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('bowtie2_path: ' + str(self.bowtie2_path))
        logger.parameter_log('macs2_path: ' + str(self.macs2_path))
        logger.parameter_log('------------------ Chipseq Load Finish ------------------')

    def bowtie2(self, ctrl, script_only):
        self.logger.write_log('Start bowtie2 script generation.')
        paired_end = ctrl.parameters['config_dict']['datasource']['chip-seq']['paired-end']
        script_path = ctrl.base_path + '/result/shell_script/bowtie2.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['gcc', 'bowtie2'], 'bowtie2_result')

        module_loader(ctrl.total_script_file, ['bowtie2'])
        total_shell_file = open(ctrl.total_script_file, 'a')

        if paired_end == 'y':
            # TODO: bowtie2 paired-end
            self.logger.write_log('bowtie2 paired-end is not implemented yet.')
        else:
            chipseq_list = loop_concatanator('n', self.config['files'])
            shell_file = open(script_path, 'a')
            shell_file.write('for i in ' + chipseq_list + '\n')
            shell_file.write('do\n')
            shell_file.write('bowtie2 -p ' + str(self.tools['bowtie2']['-p']) +
                             ' -x ' + ctrl.parameters['config_dict']['datasource']['mapping-index']['bowtie2'] +
                             ' -U ' + self.config['dir_path'] + "\"$i\".fastq.gz" +
                             ' -S ' + self.bowtie2_path + '\"$i\".sam\n')
            shell_file.write('done\n')
            shell_file.close()

            total_shell_file.write('for i in ' + chipseq_list + '\n')
            total_shell_file.write('do\n')
            total_shell_file.write('bowtie2 -p ' + str(self.tools['bowtie2']['-p']) +
                                   ' -x ' + ctrl.parameters['config_dict']['datasource']['mapping-index']['bowtie2'] +
                                   ' -U ' + self.config['dir_path'] + "\"$i\".fastq.gz" +
                                   ' -S ' + self.bowtie2_path + '\"$i\".sam\n')
            total_shell_file.write('done\n')
            total_shell_file.write('\n')
            total_shell_file.close()
            self.logger.write_log('Finish bowtie2 script generation with single-end mode.')

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)

    def macs2(self, ctrl, script_only):
        joint_ctrl = ctrl.joint_controller
        script_path = ctrl.base_path + '/result/shell_script/macs2.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['gcc/9.2.0', 'macs2'], 'macs2')
        shell_file = open(script_path, 'a')
        lb = lambda x, n: n if x == 'y' else ''
        shell_file.write('macs2 callpeak -t ' + joint_ctrl.samtools_path + '*.sort.bam' +
                            ' -c ' + joint_ctrl.samtools_path + '*.sort.bam' +
                            ' -f ' + self.tools['macs2']['-f'] + ' -g ' + self.tools['macs2']['-g'] +
                            ' -n Replicate01' +
                            ' --outdir ' + self.macs2_path + ' ' + lb(self.tools['macs2']['-bdg'], '-bdg'))
        shell_file.write('\n')
        shell_file.write('macs2 callpeak -t ' + joint_ctrl.samtools_path + '*.sort.bam' +
                            ' -c ' + joint_ctrl.samtools_path + '*.sort.bam' +
                            ' -f ' + self.tools['macs2']['-f'] + ' -g ' + self.tools['macs2']['-g'] +
                            ' -n Replicate02' +
                            ' --outdir ' + self.macs2_path + ' ' + lb(self.tools['macs2']['-bdg'], '-bdg'))
        shell_file.write('\n')
        shell_file.close()
        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)
