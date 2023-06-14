from utils.parameter.shell_builder import *


class RNASeq():
    def __init__(self, config, tools, base_path, logger):
        self.hisat2_path = base_path + '/result/rnaseq/hisat2/'
        self.stringtie_path = base_path + '/result/rnaseq/stringtie/'
        self.deseq2_path = base_path + '/result/rnaseq/deseq2/'
        self.htseq_path = base_path + '/result/rnaseq/htseq/'
        self.prepDE = base_path + '/result/rnaseq/prepDE/'
        self.getTPM = base_path + '/result/rnaseq/getTPM/'
        self.getFPKM = base_path + '/result/rnaseq/getFPKM/'
        self.config = config
        self.tools = tools
        self.logger = logger
        dir_builder(self.hisat2_path)
        dir_builder(self.stringtie_path)
        dir_builder(self.deseq2_path)
        dir_builder(self.htseq_path)
        dir_builder(self.prepDE)
        dir_builder(self.getTPM)
        dir_builder(self.getFPKM)
        logger.parameter_log('------------------ RNASeq Load Start ------------------')
        logger.parameter_log('config: ' + str(config))
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('hisat2_path: ' + str(self.hisat2_path))
        logger.parameter_log('stringtie_path: ' + str(self.stringtie_path))
        logger.parameter_log('deseq2_path: ' + str(self.deseq2_path))
        logger.parameter_log('htseq_path: ' + str(self.htseq_path))
        logger.parameter_log('prepDE_path: ' + str(self.prepDE))
        logger.parameter_log('getTPM_path: ' + str(self.getTPM))
        logger.parameter_log('getFPKM_path: ' + str(self.getFPKM))
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
        joint_ctrl = ctrl.joint_controller
        script_path = ctrl.base_path + '/result/shell_script/stringtie.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path,
                              ['gcc/9.2.0', 'stringtie'], 'stringtie')

        module_loader(ctrl.total_script_file, ['stringtie'])
        total_shell_file = open(ctrl.total_script_file, 'a')
        rnaseq_file_list = loop_concatanator('n', self.config['files'])
        shell_file = open(script_path, 'a')
        shell_file.write('for i in ' + rnaseq_file_list + '\n')
        shell_file.write('do\n')
        lb = lambda x, n: n if x == 'y' else ''
        shell_file.write('stringtie ' + lb(self.tools['stringtie']['-e'], '-e ') +
                            lb(self.tools['stringtie']['-B'], '-B ') + '-p $SLURM_CPUS_PER_TASK -G ' +
                            ctrl.annotation_source['stringtie'] + ' -o ' + self.stringtie_path + '\"$i\"' + '.gtf ' +
                            joint_ctrl.samtools_path + '\"$i\"' + '.sort.bam\n')
        shell_file.write('done\n')
        shell_file.close()

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            ctrl.logger.write_log('Finish executing shell script: ' + script_path)
        else:
            ctrl.logger.write_log('Only script is generated: ' + script_path)

    def prepDEpy(self, ctrl, script_only):
        list_path = ctrl.base_path + '/result/shell_script/list_writer.txt'
        script_path = ctrl.base_path + '/result/shell_script/prepDE.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path,
                                ['gcc/7.1.0', 'openmpi/3.1.4'], 'prepDE')
        shell_file = open(script_path, 'a')
        shell_file.write('module load python/2.7.16\n')
        current_pwd = os.getcwd()
        shell_file.write('python ' + current_pwd + '/tools/prepDE.py'
                         ' -i ' + list_path +
                         ' -g ' + self.prepDE + 'gene_count_matrix.csv' +
                         ' -t ' + self.prepDE + 'transcript_count_matrix.csv' +
                         ' -l 50')

    def getTPMpy(self, ctrl, script_only):
        list_path = ctrl.base_path + '/result/shell_script/list_writer.txt'
        script_path = ctrl.base_path + '/result/shell_script/getTPM.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path,
                                ['gcc/7.1.0', 'openmpi/3.1.4'], 'getTPM')
        shell_file = open(script_path, 'a')
        shell_file.write('module load python/2.7.16\n')
        current_pwd = os.getcwd()
        shell_file.write('python ' + current_pwd + '/tools/getTPM.py'
                         ' -i ' + list_path +
                         ' -g ' + self.prepDE + 'gene_TPM_matrix.csv' +
                         ' -t ' + self.prepDE + 'transcript_TPM_matrix.csv' +
                         ' -l 50')

    def getFPKMpy(self, ctrl, script_only):
        list_path = ctrl.base_path + '/result/shell_script/list_writer.txt'
        script_path = ctrl.base_path + '/result/shell_script/getFPKM.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path,
                                ['gcc/7.1.0', 'openmpi/3.1.4'], 'getFPKM')
        shell_file = open(script_path, 'a')
        shell_file.write('module load python/2.7.16\n')
        current_pwd = os.getcwd()
        shell_file.write('python ' + current_pwd + '/tools/getFPKM.py'
                         ' -i ' + list_path +
                         ' -g ' + self.prepDE + 'gene_FPKM_matrix.csv' +
                         ' -t ' + self.prepDE + 'transcript_FPKM_matrix.csv' +
                         ' -l 50')

    def list_writer(self, ctrl):
        script_path = ctrl.base_path + '/result/shell_script/list_writer.txt'
        # if exist, delete the old shell script
        if os.path.exists(script_path):
            os.remove(script_path)
        # create execution path
        if not os.path.exists(os.path.dirname(script_path)):
            os.makedirs(os.path.dirname(script_path))
        list_file = open(script_path, 'w')
        file_list = self.config['files']
        for item in file_list:
            list_file.write(item + ' ' + self.stringtie_path + item + '.gtf\n')
        list_file.close()


    def deseq2(self, ctrl, script_only):
        pass

    def htseq(self, ctrl, script_only):
        pass
