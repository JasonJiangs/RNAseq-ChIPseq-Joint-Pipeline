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

        pass

    def macs2(self, ctrl, script_only):
        pass



