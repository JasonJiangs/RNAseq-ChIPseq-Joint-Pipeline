from utils.parameter.shell_builder import dir_builder
import datetime


class Logger():
    def __init__(self, path):
        self.path = path + '/execution_log.txt'
        dir_builder(path)
        logfile = open(self.path, 'w')
        logfile.write('Execution Log initialized ......\n')
        logfile.close()

    def write_log(self, data: str):
        # append write log to file
        log_file = open(self.path, 'a')
        # get time
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_file.write('*' + time + '* ' + data + '\n')
        log_file.close()
        print(data)

    def parameter_log(self, data: str):
        # append write log to file
        log_file = open(self.path, 'a')
        log_file.write(data + '\n')
        log_file.close()
