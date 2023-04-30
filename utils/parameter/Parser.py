import os
import yaml

def load_parameters(config, phase, model):
    parameters = load_defaults()
    if config is None and phase is None and model is None:
        return parameters

    if config is not None:
        # read config file from config.yaml
        with open(config, "r") as stream:
            try:
                parameters['config_dict'] = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

    if phase is not None:
        parameters['phase'] = phase

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
             'resultdestination': '/scratch/pfq7pm/test_pipeline/proj_result'
             },
        "phase": "1",
        "model": "cr"
    }
    return parameters


# return four lists storing file paths of file directories
def source_file_parser(ctrl):
    rna_seq_source = {}
    chip_seq_source = {}
    mapping_index_source = {}
    annotation_source = {}

    # parse rna-seq source file paths
    rnaseq_dirpath = ctrl.parameters['config_dict']['datasource']['rna-seq']['dir_path']
    rnaseq_file_list = ctrl.parameters['config_dict']['datasource']['rna-seq']['files']
    rnaseq_suffix = ctrl.parameters['config_dict']['datasource']['rna-seq']['suffix']
    for file in rnaseq_file_list:
        if ctrl.parameters['config_dict']['datasource']['rna-seq']['paired-end'] == 'y':
            rna_seq_source[file] = [os.path.join(rnaseq_dirpath, file + '_1' + rnaseq_suffix),
                                    os.path.join(rnaseq_dirpath, file + '_2' + rnaseq_suffix)]
        else:
            rna_seq_source[file] = os.path.join(rnaseq_dirpath, file + rnaseq_suffix)

    # parse chip-seq source file paths
    chipseq_dirpath = ctrl.parameters['config_dict']['datasource']['chip-seq']['dir_path']
    chipseq_file_list = ctrl.parameters['config_dict']['datasource']['chip-seq']['files']
    chipseq_suffix = ctrl.parameters['config_dict']['datasource']['chip-seq']['suffix']
    for file in chipseq_file_list:
        if ctrl.parameters['config_dict']['datasource']['chip-seq']['paired-end'] == 'y':
            chip_seq_source[file] = [os.path.join(chipseq_dirpath, file + '_1' + chipseq_suffix),
                                     os.path.join(chipseq_dirpath, file + '_2' + chipseq_suffix)]
        else:
            chip_seq_source[file] = os.path.join(chipseq_dirpath, file + chipseq_suffix)

    # parse mapping-index source file paths
    for root, dirs, files in os.walk(ctrl.parameters['config_dict']['datasource']['mapping-index']):
        for file in files:
            # if file include 'bowtie'
            if file.find('bowtie') != -1:
                mapping_index_source['bowtie'] = os.path.join(root, file)
            elif file.find('hisat2') != -1:
                mapping_index_source['hisat2'] = os.path.join(root, file)
            elif file.find('bwa') != -1:
                mapping_index_source['bwa'] = os.path.join(root, file)
            elif file.find('star') != -1:
                mapping_index_source['star'] = os.path.join(root, file)
            elif file.find('salmon') != -1:
                mapping_index_source['salmon'] = os.path.join(root, file)
            else:
                pass

    # parse annotation source file paths
    for root, dirs, files in os.walk(ctrl.parameters['config_dict']['datasource']['annotation']):
        for file in files:
            # if file include 'gtf'
            if file.find('.gtf') != -1:
                annotation_source['gtf'] = os.path.join(root, file)
            elif file.find('.gff') != -1:
                annotation_source['gff'] = os.path.join(root, file)
            elif file.find('.bed') != -1:
                annotation_source['bed'] = os.path.join(root, file)
            elif file.find('.ref') != -1:
                annotation_source['ref'] = os.path.join(root, file)
            else:
                pass

    ctrl.rnaseq_source = rna_seq_source
    ctrl.chipseq_source = chip_seq_source
    ctrl.mapping_index_source = mapping_index_source
    ctrl.annotation_source = annotation_source

