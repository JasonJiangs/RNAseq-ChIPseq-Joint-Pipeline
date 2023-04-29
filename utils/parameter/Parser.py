from platform import system
import os
import yaml

def load_parameters(config, quality_control, model):
    parameters = load_defaults()
    if config is None and quality_control is None and model is None:
        return parameters

    if config is not None:
        # read config file from config.yaml
        with open(config, "r") as stream:
            try:
                parameters['config_dict'] = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

    if quality_control is not None:
        parameters['quality_control'] = quality_control

    if model is not None:
        parameters['model'] = model

    return parameters


def load_defaults():
    parameters = {
        "config_dict":
            {'datasource':
                 {'rna-seq': '/path/to/rna-seq',
                  'chip-seq': '/path/to/chip-seq',
                  'mapping-index': '/path/to/mapping-index',
                  'annotation': '/path/to/annotation'
                  },
             'resultdestination':
                 {'rna-seq': '/path/to/rna-seq-results',
                  'chip-seq': '/path/to/chip-seq-results',
                  'joint-analysis': '/path/to/joint-analysis-results',
                  'modeling': '/path/to/modeling-results',
                  'log': '/path/to/log'}
             },
        "quality_control": "oa",
        "model": "cr",
        "OS": system(),
        "sys_path": os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    }
    return parameters


