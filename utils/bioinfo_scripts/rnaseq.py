from utils.parameter.shell_builder import *


class RNASeq():
    def __init__(self, config, tools, base_path, logger):
        self.hisat2_path = base_path + '/result/rnaseq/hisat2/'
        self.stringtie_path = base_path + '/result/rnaseq/stringtie/'
        self.deseq2_path = base_path + '/result/rnaseq/deseq2/'
        self.htseq_path = base_path + '/result/rnaseq/htseq/'
        self.prepDE = base_path + '/result/rnaseq/prepDE/'
        self.config = config
        self.tools = tools
        self.logger = logger
        dir_builder(self.hisat2_path)
        dir_builder(self.stringtie_path)
        dir_builder(self.deseq2_path)
        dir_builder(self.htseq_path)
        dir_builder(self.prepDE)
        logger.parameter_log('------------------ RNASeq Load Start ------------------')
        logger.parameter_log('config: ' + str(config))
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('hisat2_path: ' + str(self.hisat2_path))
        logger.parameter_log('stringtie_path: ' + str(self.stringtie_path))
        logger.parameter_log('deseq2_path: ' + str(self.deseq2_path))
        logger.parameter_log('htseq_path: ' + str(self.htseq_path))
        logger.parameter_log('prepDE_path: ' + str(self.prepDE))
        logger.parameter_log('------------------ RNASeq Load Finish ------------------')

    def hisat2(self, ctrl, script_only):
        self.logger.write_log('Start hisat2 script generation.')
        paired_end = ctrl.parameters['config_dict']['datasource']['rna-seq']['paired-end']
        script_path = ctrl.base_path + '/result/shell_script/hisat2.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path,
                              ['gcc', 'hisat2'], 'hisat2_result')

        module_loader(ctrl.total_script_file, ['hisat2'])
        total_shell_file = open(ctrl.total_script_file, 'a')

        if paired_end == 'y':
            rnaseq_file_list = loop_concatanator('n', self.config['files'])
            shell_file = open(script_path, 'a')
            shell_file.write('for i in ' + rnaseq_file_list + '\n')
            shell_file.write('do\n')
            shell_file.write('hisat2 -t -p ' + str(self.tools['hisat2']['-p']) + ' -x ' +
                             ctrl.mapping_index_list['hisat2'] + ' -1 ' + self.config['dir_path'] +
                             '\"$i\"_1.fastq.gz -2 ' + self.config['dir_path'] + '\"$i\"_2.fastq.gz -S '
                             + self.hisat2_path + '\"$i\"' + '.sam\n')
            shell_file.write('done\n')
            shell_file.close()

            total_shell_file.write('for i in ' + rnaseq_file_list + '\n')
            total_shell_file.write('do\n')
            total_shell_file.write('hisat2 -t -p ' + str(self.tools['hisat2']['-p']) + ' -x ' +
                                   ctrl.mapping_index_list['hisat2'] + ' -1 ' + self.config['dir_path'] +
                                   '\"$i\"_1.fastq.gz -2 ' + self.config['dir_path'] + '\"$i\"_2.fastq.gz -S '
                                   + self.hisat2_path + '\"$i\"' + '.sam\n')
            total_shell_file.write('done\n')
            total_shell_file.write('\n')
            total_shell_file.close()
            self.logger.write_log('Finish hisat2 script generation with paired-end mode.')
        else:
            # TODO: single end
            self.logger.write_log('Finish hisat2 script generation with single-end mode.')

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
                              ['gcc', 'stringtie'], 'stringtie_result')

        module_loader(ctrl.total_script_file, ['stringtie'])
        total_shell_file = open(ctrl.total_script_file, 'a')

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)
        else:
            ctrl.logger.write_log('Only script is generated: ' + script_path)

    def prepDE(self, ctrl, script_only):
        prepDE_path = os.getcwd() + '/tools/prepDE.py'
        module_loader(ctrl.total_script_file, ['gcc/7.1.0', 'openmpi/3.1.4', 'python/2.7.16'])
        total_shell_file = open(ctrl.total_script_file, 'a')

        # delete prepDE_list.txt if exist
        prepDE_list_path = self.prepDE + 'prepDE_list.txt'
        if os.path.isfile(prepDE_list_path):
            os.remove(prepDE_list_path)

        # make prepDE_list.txt
        prepDE_list_file = open(prepDE_list_path, 'w')
        # get all files in stringtie_path
        stringtie_file_list = os.listdir(self.stringtie_path)
        for i in range(len(stringtie_file_list)):
            file = stringtie_file_list[i]
            if file.endswith('.gtf'):
                file_name = file.split('.')[0]
                if i == len(stringtie_file_list) - 1:
                    prepDE_list_file.write(file_name + ' ' + self.stringtie_path + file)
                else:
                    prepDE_list_file.write(file_name + ' ' + self.stringtie_path + file + '\n')
        prepDE_list_file.close()

        total_shell_file.write('python ' + prepDE_path + ' -i ' + prepDE_list_path)
        total_shell_file.close()

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + ctrl.total_script_file)
            linux_command = 'sbatch ' + ctrl.total_script_file
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + ctrl.total_script_file)
        else:
            ctrl.logger.write_log('Only script is generated: ' + ctrl.total_script_file)

    def deseq2(self, ctrl, script_only):
        pass

    def htseq(self, ctrl, script_only):
        pass
