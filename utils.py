import importlib
import argparse
import os
import settings
import sys


def import_module_from_string(class_path):
    """
    Try to import a class and return it
    Takes a python-style path e.g. my_module.foo.Bar if we want class Bar from the my_module/foo.py module
    """
    try:
        parts = class_path.split('.')
        module_path, class_name = '.'.join(parts[:-1]), parts[-1]
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    except (ImportError, AttributeError) as ex:
        msg = "Could not import: %s" % class_path
        raise ImportError(msg)


class CLParser(object):

    class MakeAbsolutePathAction(argparse.Action):
        '''
            this attaches an action (if specified with 'action=MakeAbsolutePathAction')
            relative paths are converted to absolute paths for strictness
        '''
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, os.path.realpath(os.path.abspath(values)))

    class MakeAbsolutePathActionForList(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            s = [os.path.realpath(os.path.abspath(v)) for v in values]
            setattr(namespace, self.dest, s)

    @classmethod
    def setup_args(cls):
        parser = argparse.ArgumentParser(description = 'Creates a figure showing read coverage')
        parser.add_argument('-t',
                            '--transcript',
                            required = True,
                            default = None,
                            help = 'The transcript identifier (e.g. ENST if Ensembl GTF file.)',
                            dest = 'selected_transcript')

        parser.add_argument('-a',
                            '--annotation',
                            required = True,
                            default = None,
                            action = cls.MakeAbsolutePathAction,
                            help = 'Path to an annotation file (e.g. GTF) which defines the transcript features',
                            dest = 'annotation_filepath')

        parser.add_argument('-u',
                            '--upstream_pad',
                            required = False,
                            default = 0,
                            help = 'Specifies how many base pairs to pad the transcript by in the upstream direction.',
                            type = int,
                            dest = 'upstream_padding')

        parser.add_argument('-d',
                            '--downstream_pad',
                            required = False,
                            default = 0,
                            help = 'Specifies how many base pairs to pad the transcript by in the downstream direction.',
                            type = int,
                            dest = 'downstream_padding')

        parser.add_argument('-o',
                            '--output',
                            required = True,
                            action = cls.MakeAbsolutePathAction,
                            help = 'Path to the output file.  The output format will be implied by the file extension',
                            dest = 'output_path')

        parser.add_argument('coverage_files',
                            nargs = '+',
                            action = cls.MakeAbsolutePathActionForList,
                            help = 'Any number of paths to input files, separated by space.  These can '
                                   'be any type of files that contain coverage information and are understood by this program.  '
                                   'Note that the coverage overlay order matches the order here.')

        """
        parser.add_argument('-no_gene',
                            required = False,
                            default = None,
                            action = 'store_true',
                            dest = 'include_gene_model')
        """

        return parser

    @classmethod
    def parse_cl_args(cls):
        parser = cls.setup_args()
        return vars(parser.parse_args())


class UnsupportedFiletype(Exception):
    pass


class MultipleFiletypeException(Exception):
    pass


def locate_provider(file_ext, available_classes):
    '''
    Returns the class for a particular file extension.
    '''
    for c in available_classes:
        print 'try class %s' % c
        clz = import_module_from_string(c)
        if file_ext == clz.ACCEPTED_FILE_EXTENSION:
            return clz
    raise UnsupportedFiletype("Could not find an appropriate class for reading your file with extension: '%s' " % file_ext)


def get_data_source_provider(f, source_type, many=False):

    if not many:
        file_ext = f.split('.')[-1].lower()
    else:
        file_ext_set = set([x.split('.')[-1].lower() for x in f])
        if len(file_ext_set) != 1:
            raise MultipleFiletypeException('Multiple filetypes is not supported at the moment.')
        file_ext = list(file_ext_set)[0]

    if source_type == settings.ANNOTATION_SRC:
        available_classes = settings.GENE_DETAIL_DATA_SOURCES
    elif source_type == settings.COVERAGE_SRC:
        available_classes = settings.COVERAGE_DATA_SOURCES
    else:
        #TODO: custom exception and/or message
        raise Exception('??')
    clz = locate_provider(file_ext, available_classes)
    try:
        return clz(f, many=many)
    except Exception as ex:
        print 'Error: could not create a data source.'
        sys.exit(1)


def plot_gene_model(**kwargs):
    # see if they explicitly specified whether the gene model should be in the image
    # If not, fall back to the default from the settings module.
    try:
        return kwargs['include_gene_model']
    except KeyError:
        return settings.INCLUDE_GENE_MODEL


def configure_data_provider(**kwargs):
    try:

        # instantiate the class that will supply data for the plot.  We will then
        # inject the appropriate elements into that.
        data_provider_clz = import_module_from_string(settings.DATA_PROVIDER)
        data_provider = data_provider_clz()

        # annotation data- has the region to look at
        annotation_filepath = kwargs['annotation_filepath']
        annotation_source = get_data_source_provider(annotation_filepath, settings.ANNOTATION_SRC)
        data_provider.add_data_source(settings.ANNOTATION_SRC, annotation_source)

        # files containing coverage info
        input_files = kwargs['coverage_files']
        coverage_data_source = get_data_source_provider(input_files, settings.COVERAGE_SRC, many=True)
        data_provider.add_data_source(settings.COVERAGE_SRC, coverage_data_source)

        return data_provider

    except Exception as ex:
        #TODO: catch specific exceptions, handle appropriately.
        raise ex


def get_artists(**kwargs):
    cvg_artist_clz = import_module_from_string(settings.COVERAGE_PROFILE_ARTIST)

    gene_detail_artist_clz = None
    try:
        draw_gene = kwargs['include_gene_model']
    except KeyError:
        draw_gene = settings.INCLUDE_GENE_MODEL
    if draw_gene:
        gene_detail_artist_clz = import_module_from_string(settings.GENE_DETAIL_ARTIST)

    additional_artists = [import_module_from_string(x) for x in settings.OTHER_ARTISTS]

    return cvg_artist_clz, gene_detail_artist_clz, additional_artists

