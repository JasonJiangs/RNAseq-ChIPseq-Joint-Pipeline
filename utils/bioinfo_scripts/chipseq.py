class Chipseq():
    def __init__(self, config, tools, logger):
        self.config = config
        self.tools = tools

        logger.parameter_log('------------------ Chipseq Load Start ------------------')
        logger.parameter_log('config: ' + str(config))
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('------------------ Chipseq Load Finish ------------------')


def bowtie2(ctrl, script_only):
    pass
