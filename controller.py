from utils.parameter.Parser import source_file_parser
from utils.bioinfo_scripts.joint import Joint, fastqc, samtools
from utils.bioinfo_scripts.chipseq import Chipseq, bowtie2
from utils.bioinfo_scripts.rnaseq import RNASeq, hisat2
from utils.clustering.cluster_model import Cluster
from utils.regression.regression_model import Regression

class Controller():
    def __init__(self, parameters, logger):
        self.parameters = parameters
        self.model = parameters['model']
        self.phase = parameters['phase']
        self.logger = logger
        self.rnaseq_source = None
        self.chipseq_source = None
        self.mapping_index_source = None
        self.annotation_source = None
        self.rnaseq_controller = RNASeq(parameters['config_dict']['datasource']['rna-seq'], parameters['config_dict']['rnaseq'])
        self.chipseq_controller = Chipseq(parameters['config_dict']['datasource']['chip-seq'], parameters['config_dict']['chipseq'])
        self.joint_controller = Joint(parameters['config_dict']['joint'])
        self.slurm = parameters['config_dict']['slurm']
        self.base_path = parameters['config_dict']['resultdestination']

    def execute(self):
        # parse source file paths
        source_file_parser(self)

        if self.phase == '1':
            self.logger.write_log(self, "Start to run phase 1.")
            # phase 1
            phase1_execution(self)
            self.logger.write_log(self, "End to run phase 1.")

        elif self.phase == '12':
            self.logger.write_log(self, "Start to run phase 12.")
            phase1_execution(self)
            phase2_execution(self)
            self.logger.write_log(self, "End to run phase 12.")

        elif self.phase == '123':
            self.logger.write_log(self, "Start to run phase 123.")
            phase1_execution(self)
            phase2_execution(self)
            phase3_execution(self)
            self.logger.write_log(self, "End to run phase 123")

        elif self.phase == '2':
            self.logger.write_log(self, "Start to run phase 2.")
            phase2_execution(self)
            self.logger.write_log(self, "End to run phase 2.")

        elif self.phase == '3':
            self.logger.write_log(self, "Start to run phase 3.")
            phase3_execution(self)
            self.logger.write_log(self, "End to run phase 3.")

        else:
            self.logger.write_log(self, "Invalid phase number.")
            exit(1)


def phase1_execution(ctrl):
    # quality control before read alignment
    fastqc(ctrl, 'before')
    # read alignment
    bowtie2(ctrl)
    hisat2(ctrl)
    # quality control after read alignment
    fastqc(ctrl, 'after')

def phase2_execution(ctrl):
    samtools(ctrl)
    pass

def phase3_execution(ctrl):
    pass

def model_execution(ctrl):
    pass