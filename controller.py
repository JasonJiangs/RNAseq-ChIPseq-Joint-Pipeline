
class Controller():
    def __init__(self, parameters, logger):
        self.parameters = parameters
        self.model = None
        self.view = None
        self.quality_control = None
        self.logger = logger

    def execute(self):
        pass