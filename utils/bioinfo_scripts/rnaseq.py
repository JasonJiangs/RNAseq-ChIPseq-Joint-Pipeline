class RNASeq():
    def __init__(self, config, tools, logger):
        self.config = config
        self.tools = tools
        logger.parameter_log('------------------ RNASeq Load Start ------------------')
        logger.parameter_log('config: ' + str(config))
        logger.parameter_log('tools: ' + str(tools))
        logger.parameter_log('------------------ RNASeq Load Finish ------------------')


def hisat2(ctrl, script_only):
    pass
