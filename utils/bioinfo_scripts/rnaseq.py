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
        pass

    def stringtie(self, ctrl, script_only):
        pass

    def deseq2(self, ctrl, script_only):
        pass

    def htseq(self, ctrl, script_only):
        pass
