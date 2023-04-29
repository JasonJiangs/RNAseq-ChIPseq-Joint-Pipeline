import argparse
from utils.parameter.Parser import load_parameters
from controller import Controller
from utils.logger.logger import Logger

def main():
    # intake command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest='config', help="Config file path")
    parser.add_argument("-qc", dest='quality_control', help="quality control running choice")
    parser.add_argument("-m", dest='model', help="Model running choice")
    args = parser.parse_args()
    parameters = load_parameters(args.config, args.quality_control, args.model)

    # create logger
    logger = Logger()

    # create controller
    ctrl = Controller(parameters=parameters, logger=logger)

    logger.write_log("Finish loading parameters, start to execute.", parameters['config_dict']['resultdestination']['log'])

    ctrl.execute()

    logger.write_log("Finish executing.", parameters['config_dict']['resultdestination']['log'])


if __name__ == '__main__':
    main()