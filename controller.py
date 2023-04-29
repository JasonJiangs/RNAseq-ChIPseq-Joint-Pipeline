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
        self.quality_control = parameters['quality_control']
        self.logger = logger
        self.rnaseq_source = None
        self.chipseq_source = None
        self.mapping_index_source = None
        self.annotation_source = None
        self.rnaseq_controller = RNASeq()
        self.chipseq_controller = Chipseq()

    def execute(self):
        # parse source file paths
        source_file_parser(self)

        if self.quality_control == 'o':
            self.logger.write_log(self, "Start quality control.")

            # quality control before read alignment
            fastqc(self)
            # read alignment
            bowtie2(self)
            hisat2(self)
            # quality control after read alignment
            fastqc(self)

            self.logger.write_log(self, "Finish quality control.")
        elif self.quality_control == 'a':
            self.logger.write_log(self, "Start workflows after quality control.")

            self.logger.write_log(self, "End workflows after quality control.")
        elif self.quality_control == 'oa':
            self.logger.write_log(self, "Start the whole workflow.")
            fastqc(self)

            self.logger.write_log(self, "Finish the whole workflow.")
        else:
            self.logger.write_log(self, "Exit.")
