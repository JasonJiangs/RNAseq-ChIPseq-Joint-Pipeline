from utils.parameter.shell_builder import *
from utils.parameter.html_parser import *


class Joint():
    def __init__(self, tools, bash_path, logger):
        self.qc_path = bash_path + '/result/joint/qc/'
        self.samtools_path = bash_path + '/result/joint/samtools/'
        self.fastp_path = bash_path + '/result/joint/fastp/'
        self.tools = tools
        self.logger = logger
        dir_builder(self.qc_path)
        dir_builder(self.samtools_path)
        dir_builder(self.fastp_path)
        logger.parameter_log('------------------ Joint Load Start ------------------')
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('bash_path: ' + bash_path)
        logger.parameter_log('qc_path: ' + str(self.qc_path))
        logger.parameter_log('samtools_path: ' + str(self.samtools_path))
        logger.parameter_log('fastp_path: ' + str(self.fastp_path))
        logger.parameter_log('------------------ Joint Load Finish ------------------')

    def fastp(self, ctrl, script_only):
        self.logger.write_log('Start to write fastp.sh script.')
        chipseq_ctrl = ctrl.chipseq_controller
        rnaseq_ctrl = ctrl.rnaseq_controller

        if rnaseq_ctrl.config['paired-end']:
            file_name = 'fastp'
            script_path = ctrl.base_path + '/result/shell_script/fastp.sh'
        else:
            file_name = 'fastp'
            script_path = ctrl.base_path + '/result/shell_script/fastp.sh'

        # RNA-seq file
        if rnaseq_ctrl.config['paired-end'] == 'y':
            first_rnaseq_path = rnaseq_ctrl.config['dir_path_original'] + rnaseq_ctrl.config['files'][0] + '_1.fastq.gz'
        else:
            first_rnaseq_path = rnaseq_ctrl.config['dir_path_original'] + rnaseq_ctrl.config['files'][0] + '.fastq.gz'
        # ChIP-seq file
        if chipseq_ctrl.config['paired-end'] == 'y':
            first_chipseq_path = chipseq_ctrl.config['dir_path_original'] + chipseq_ctrl.config['files'][0] + '_1.fastq.gz'
        else:
            first_chipseq_path = chipseq_ctrl.config['dir_path_original'] + chipseq_ctrl.config['files'][0] + '.fastq.gz'
        anaconda_general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['fastp'], file_name)
        shell_file = open(script_path, 'a')

        shell_file.write('# check the read length for RNA-seq\n')
        shell_file.write('zcat ' + first_rnaseq_path + ' | awk \'{if(NR==2) {print length($1); exit}}\'\n')
        shell_file.write('# check the read length for ChIP-seq\n')
        shell_file.write('zcat ' + first_chipseq_path + ' | awk \'{if(NR==2) {print length($1); exit}}\'\n')
        shell_file.write('\n')

        if rnaseq_ctrl.config['paired-end'] == 'y':
            shell_file.write('for i in ' + loop_concatanator('n', rnaseq_ctrl.config['files']) + '\n')
            shell_file.write('do\n')
            shell_file.write('fastp -i ' + rnaseq_ctrl.config['dir_path_original'] + '\"$i\"' + '_1.fastq.gz \\\n')
            shell_file.write('-I ' + rnaseq_ctrl.config['dir_path_original'] + '\"$i\"' + '_2.fastq.gz \\\n')
            shell_file.write('-o ' + rnaseq_ctrl.config['dir_path'] + '\"$i\"' + '_1.fastq.gz \\\n')
            shell_file.write('-O ' + rnaseq_ctrl.config['dir_path'] + '\"$i\"' + '_2.fastq.gz \\\n')
            shell_file.write('--cut_by_quality3 --cut_front \\\n')
            shell_file.write('--cut_mean_quality ' + str(self.tools['fastp']['cut_mean_quality']) +
                             ' --cut_window_size ' + str(self.tools['fastp']['cut_window_size']) + ' \\\n')
            shell_file.write('done\n')
        else:
            pass

        shell_file.write('\n')
        if chipseq_ctrl.config['paired-end'] == 'y':
            pass
        else:
            shell_file.write('for i in ' + loop_concatanator('n', chipseq_ctrl.config['files']) + '\n')
            shell_file.write('do\n')
            shell_file.write('fastp -i ' + chipseq_ctrl.config['dir_path_original'] + '\"$i\"' + '.fastq.gz \\\n')
            shell_file.write('-o ' + chipseq_ctrl.config['dir_path'] + '\"$i\"' + '.fastq.gz \\\n')
            shell_file.write('--cut_by_quality3 --cut_front \\\n')
            shell_file.write('--cut_mean_quality ' + str(self.tools['fastp']['cut_mean_quality']) +
                             ' --cut_window_size ' + str(self.tools['fastp']['cut_window_size']) + ' \\\n')
            shell_file.write('done\n')

        shell_file.close()
        self.logger.write_log('Finish to write fastp.sh script.')



    def fastqc(self, ctrl, condition, script_only):
        chipseq_ctrl = ctrl.chipseq_controller
        rnaseq_ctrl = ctrl.rnaseq_controller

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

        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['gcc', 'fastqc'], file_name)

        module_loader(ctrl.total_script_file, ['gcc', 'fastqc'])
        total_shell_file = open(ctrl.total_script_file, 'a')

        chipseq_file_str = loop_concatanator(chipseq_ctrl.config['paired-end'], chipseq_ctrl.config['files'])
        shell_file = open(script_path, 'a')
        shell_file.write('for i in ' + chipseq_file_str + '\n')
        shell_file.write('do\n')
        shell_file.write('fastqc -f ' + file_format + ' -o ' + self.qc_path + condition + '/ ' +
                         chipseq_ctrl.config['dir_path'] + '\"$i\"' + chipseq_ctrl.config['suffix'] + '\n')
        shell_file.write('done\n')

        total_shell_file.write('for i in ' + chipseq_file_str + '\n')
        total_shell_file.write('do\n')
        total_shell_file.write('fastqc -f ' + file_format + ' -o ' + self.qc_path + condition + '/ ' +
                               chipseq_ctrl.config['dir_path'] + '\"$i\"' + chipseq_ctrl.config['suffix'] + '\n')
        total_shell_file.write('done\n')
        total_shell_file.write('\n')
        dir_builder(self.qc_path + condition + '/')

        rnaseq_file_str = loop_concatanator(rnaseq_ctrl.config['paired-end'], rnaseq_ctrl.config['files'])
        shell_file.write('for i in ' + rnaseq_file_str + '\n')
        shell_file.write('do\n')
        shell_file.write('fastqc -f ' + file_format + ' -o ' + self.qc_path + condition + '/ ' +
                         rnaseq_ctrl.config['dir_path'] + '\"$i\"' + rnaseq_ctrl.config['suffix'] + '\n')
        shell_file.write('done\n')

        total_shell_file.write('for i in ' + rnaseq_file_str + '\n')
        total_shell_file.write('do\n')
        total_shell_file.write('fastqc -f ' + file_format + ' -o ' + self.qc_path + condition + '/ ' +
                               rnaseq_ctrl.config['dir_path'] + '\"$i\"' + rnaseq_ctrl.config['suffix'] + '\n')
        total_shell_file.write('done\n')
        total_shell_file.write('\n')
        dir_builder(self.qc_path + condition + '/')

        shell_file.close()
        total_shell_file.close()

        ctrl.logger.write_log('fastqc script build finished: ' + script_path)

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)

            # summarize html
            ctrl.logger.write_log('Finish FastQC and start summarizing html files for ' + condition)
            dir_path = self.qc_path + condition + '/'
            # loop files in the dir
            for file in os.listdir(dir_path):
                if file.endswith('.html'):
                    summarize_html(ctrl=ctrl,
                                   logger=ctrl.logger,
                                   sample_name=file.split('.')[0],
                                   sample_path=dir_path + file)

    def multiqc(self, ctrl, script_only):
        pass

    def samtools(self, ctrl, script_only):
        self.logger.write_log('Start building samtools script')
        chipseq_ctrl = ctrl.chipseq_controller
        rnaseq_ctrl = ctrl.rnaseq_controller

        script_path = ctrl.base_path + '/result/shell_script/samtools.sh'
        general_shell_builder(ctrl.slurm, script_path, ctrl.slurm_log_path, ['gcc', 'samtools'], 'samtools')

        module_loader(ctrl.total_script_file, ['gcc', 'samtools'])
        total_shell_file = open(ctrl.total_script_file, 'a')

        chipseq_file_str = loop_concatanator(chipseq_ctrl.config['paired-end'], chipseq_ctrl.config['files'])
        rnaseq_file_str = loop_concatanator(rnaseq_ctrl.config['paired-end'], rnaseq_ctrl.config['files'])
        shell_file = open(script_path, 'a')
        shell_file.write('for i in ' + chipseq_file_str + '\n')
        shell_file.write('do\n')
        shell_file.write('samtools view -@ $SLURM_CPUS_PER_TASK -b -S ' +
                         chipseq_ctrl.config['dir_path'] + '\"$i\".sam' + ' > ' +
                         self.samtools_path + '\"$i\"' + '.bam\n')
        shell_file.write('samtools sort -@ $SLURM_CPUS_PER_TASK ' + self.samtools_path + '\"$i\".bam' +
                            ' -o ' + self.samtools_path + '\"$i\".sort.bam\n')
        shell_file.write('samtools index ' + self.samtools_path + '\"$i\".sorted.bam\n')
        shell_file.write('done\n')

        shell_file.write('\n')

        shell_file.write('for i in ' + rnaseq_file_str + '\n')
        shell_file.write('do\n')
        shell_file.write('samtools view -@ $SLURM_CPUS_PER_TASK -b -S ' +
                            rnaseq_ctrl.config['dir_path'] + '\"$i\".sam' + ' > ' +
                            self.samtools_path + '\"$i\"' + '.bam\n')
        shell_file.write('samtools sort -@ $SLURM_CPUS_PER_TASK ' + self.samtools_path + '\"$i\".bam' +
                            ' -o ' + self.samtools_path + '\"$i\".sort.bam\n')
        shell_file.write('samtools index ' + self.samtools_path + '\"$i\".sorted.bam\n')
        shell_file.write('done\n')

        ctrl.logger.write_log('samtools script build finished: ' + script_path)

        if script_only == 'n':
            ctrl.logger.write_log('Start executing shell script: ' + script_path)
            linux_command = 'sbatch ' + script_path
            shell_runner(linux_command, ctrl)
            self.logger.write_log('samtools execution finished' + script_path)
        else:
            self.logger.write_log('samtools script only mode')
