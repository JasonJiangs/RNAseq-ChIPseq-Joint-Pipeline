from utils.parameter.Parser import source_file_parser
from utils.parameter.shell_builder import dir_builder, total_log_general_shell_builder
from utils.bioinfo_scripts.joint import Joint
from utils.bioinfo_scripts.chipseq import Chipseq
from utils.bioinfo_scripts.rnaseq import RNASeq


class Controller():
    def __init__(self, parameters, logger):
        self.logger = logger
        self.logger.parameter_log('==================== Controller Load Start ====================')
        self.logger.parameter_log('Parameters:')
        self.parameters = parameters
        self.model = parameters['model']
        self.phase = parameters['phase']
        self.rnaseq_source = None
        self.chipseq_source = None
        self.mapping_index_source = None
        self.annotation_source = None
        self.base_path = parameters['config_dict']['resultdestination']
        self.slurm_log_path = self.base_path + '/result/shell_log/'
        self.rnaseq_controller = RNASeq(parameters['config_dict']['datasource']['rna-seq'],
                                        parameters['config_dict']['rnaseq'], self.base_path, self.logger)
        self.chipseq_controller = Chipseq(parameters['config_dict']['datasource']['chip-seq'],
                                          parameters['config_dict']['chipseq'], self.base_path, self.logger)
        self.joint_controller = Joint(parameters['config_dict']['joint'], self.base_path, self.logger)
        self.slurm = parameters['config_dict']['slurm']
        self.total_log = self.base_path + '/result/total_log.txt'
        self.shell_script_path = self.base_path + '/result/shell_script/'
        source_file_parser(self)
        dir_builder(self.base_path)
        dir_builder(self.slurm_log_path)
        dir_builder(self.total_log)
        dir_builder(self.shell_script_path)
        total_log_general_shell_builder(self.total_log, self.slurm)
        self.logger.parameter_log('model: ' + self.model)
        self.logger.parameter_log('phase: ' + self.phase)
        self.logger.parameter_log('rnaseq_source: ' + str(self.rnaseq_source))
        self.logger.parameter_log('chipseq_source: ' + str(self.chipseq_source))
        self.logger.parameter_log('mapping_index_source: ' + str(self.mapping_index_source))
        self.logger.parameter_log('annotation_source: ' + str(self.annotation_source))
        self.logger.parameter_log('slurm: ' + str(self.slurm))
        self.logger.parameter_log('base_path: ' + str(self.base_path))
        self.logger.parameter_log('==================== Controller Load Finish ====================')

    def execute(self):
        # if script only mode
        script_only = self.parameters['config_dict']['script-only']

        if self.phase == '1':
            self.logger.write_log("Start to run phase 1.")
            # phase 1
            phase1_execution(self, script_only)
            self.logger.write_log("End to run phase 1.")

        elif self.phase == '12':
            self.logger.write_log("Start to run phase 12.")
            phase1_execution(self, script_only)
            phase2_execution(self, script_only)
            self.logger.write_log("End to run phase 12.")

        elif self.phase == '123':
            self.logger.write_log("Start to run phase 123.")
            phase1_execution(self, script_only)
            phase2_execution(self, script_only)
            phase3_execution(self, script_only)
            self.logger.write_log("End to run phase 123")

        elif self.phase == '2':
            self.logger.write_log("Start to run phase 2.")
            phase2_execution(self, script_only)
            self.logger.write_log("End to run phase 2.")

        elif self.phase == '3':
            self.logger.write_log("Start to run phase 3.")
            phase3_execution(self, script_only)
            self.logger.write_log("End to run phase 3.")

        else:
            self.logger.write_log("Invalid phase number.")
            exit(-1)


def phase1_execution(ctrl, script_only):
    # quality control before read alignment
    ctrl.joint_controller.fastqc(ctrl, 'before', script_only)
    # read alignment
    ctrl.chipseq_controller.bowtie2(ctrl, script_only)
    ctrl.rnaseq_controller.hisat2(ctrl, script_only)
    # quality control after read alignment
    ctrl.joint_controller.fastqc(ctrl, 'after', script_only)


def phase2_execution(ctrl, script_only):
    pass


def phase3_execution(ctrl, script_only):
    pass


def model_execution(ctrl, script_only):
    pass
