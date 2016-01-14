import matplotlib.pyplot as plt
import settings


class CoveragePlotManager(object):

    def __init__(self):

        # artists that are dependent on data being supplied to them:
        self.coverage_artist = None
        self.gene_detail_artist = None

        # we can have multiple artists that are independent of data
        # for example, there could be artists for formatting axes ticks, etc.
        self.generic_artists = []

        # Initialize:
        self.data_provider = None
        self.figure, self.ax = plt.subplots()

    def set_data_provider(self, data_provider):
        self.data_provider = data_provider

    def add_artists(self, cvg_artist_clz, gene_detail_artist_clz, *other_artists):
        '''
        artist_clz should be a class (e.g. not an instance)
        Here, we instantiate it and give it 'access' to the axis on which to draw
        '''
        self.coverage_artist = cvg_artist_clz(self.ax)

        if gene_detail_artist_clz:
            self.gene_detail_artist = gene_detail_artist_clz(self.ax)

        self.generic_artists.extend([clz(self.ax) for clz in other_artists])

    def draw(self, **kwargs):
        '''
        The main draw method
        '''

        # get the transcript information
        annotation_data_source = self.data_provider.get_data_source(settings.ANNOTATION_SRC)
        window_interval = annotation_data_source.get_window_interval(kwargs['transcript_id'],
                                                                     kwargs['upstream_padding'],
                                                                     kwargs['downstream_padding'])

        # get the coverage data


        # get the transcript features so the artist knows where to draw
        if self.gene_detail_artist:
            transcript_features = annotation_data_source.get_genomic_feature_set()

        

        self.coverage_artist.draw()
        self.gene_detail_artist.draw()
        for a in self.generic_artists:
            a.draw()

    def save_figure(self, output_filepath):
        self.figure.savefig(output_filepath)