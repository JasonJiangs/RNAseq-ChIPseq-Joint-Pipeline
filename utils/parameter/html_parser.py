from lxml import etree


def summarize_html(ctrl, logger, sample_name, sample_path):
    logger.write_log('=========> FastQC brief for ' + sample_name + ': ')
    logger.parameter_log('File path: ' + sample_path)
    passed, warned, failed, h2_alt_dict = html_conclude(sample_path, logger)
    # stop if the number of fail and warn is exceed the threshold

    # log the result respectively
    for key in h2_alt_dict:
        logger.parameter_log(key + ': ' + h2_alt_dict[key])

    logger.write_log('In conclusion: ' + str(passed) + ' passed, ' +
                     str(warned) + ' warned, ' + str(failed) + ' failed  <=========')

    if failed > ctrl.parameters['config_dict']['qc-toleration']['fail'] or \
            warned > ctrl.parameters['config_dict']['qc-toleration']['warn']:
        logger.write_log('Error: the number of fail or warn is exceed the threshold.')
        exit(-1)


def html_conclude(path, logger):
    html = etree.parse(path, etree.HTMLParser())
    # get the content as a list from every h2 tag
    h2_list = html.xpath('//h2/text()')[1:]
    # get the content of the alt as a list from every img inside a h2 tag
    alt_list = html.xpath('//h2/img/@alt')

    # exist if the length of h2_list and alt_list is not equal
    if len(h2_list) != len(alt_list):
        logger.write_log('Error: the length of h2_list and alt_list is not equal.')
        exit(-1)
    # if the length of h2_list and alt_list is 0, return 0
    if len(h2_list) and len(alt_list) == 0:
        logger.write_log('Error: the length of h2_list and alt_list is 0.')
        exit(-1)

    # combine h2_list and alt_list into a dictionary
    h2_alt_dict = dict(zip(h2_list, alt_list))
    # count the number of passed, warned and failed
    passed = 0
    warned = 0
    failed = 0
    for key in h2_alt_dict:
        if h2_alt_dict[key] == 'pass':
            passed += 1
        elif h2_alt_dict[key] == 'warn':
            warned += 1
        elif h2_alt_dict[key] == 'fail':
            failed += 1
        return passed, warned, failed, h2_alt_dict
