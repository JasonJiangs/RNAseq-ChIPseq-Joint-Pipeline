import argparse
from utils.parameter.parser import load_parameters
from controller import Controller
from utils.logger.logger import Logger


# debug mode under linux
# import pdb
# pdb.set_trace()

def main():
    # intake command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest='config', help="Config file path")
    parser.add_argument("-p", dest='phase', help="choose which phase to run")
    parser.add_argument("-m", dest='model', help="Model running choice")
    args = parser.parse_args()
    parameters = load_parameters(args.config, args.phase, args.model)

    # create logger
    logger = Logger(parameters['config_dict']['resultdestination'])

    logger.write_log("Start loading parameters.")

    # create controller
    ctrl = Controller(parameters=parameters, logger=logger)

    logger.write_log("Finish loading parameters, start to execute.")

    ctrl.execute()

    logger.write_log("Finish executing.")


if __name__ == '__main__':
    main()
