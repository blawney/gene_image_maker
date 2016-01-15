import matplotlib
print matplotlib.get_backend()
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
        self.fig = plt.figure(figsize=settings.FIGURE_SIZE)
        self.ax = self.fig.add_subplot(111)

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
        print 'ann data src: %s' % annotation_data_source
        try:
            print kwargs
            window_interval = annotation_data_source.get_window_interval(kwargs['selected_transcript'],
                                                                     kwargs['upstream_padding'],
                                                                     kwargs['downstream_padding'])
        except Exception as ex:
            print 'EXXXXX'
            print ex.message
        print 'window: %s' % window_interval

        # get the coverage data as a Pandas DataFrame
        # each row is the coverage for a single sample
        coverage_data_source = self.data_provider.get_data_source(settings.COVERAGE_SRC)
        print coverage_data_source
        cvg_data_df = coverage_data_source.get_coverage_data(window_interval)
        self.coverage_artist.draw(cvg_data_df)


        # get the transcript features so the artist knows where to draw
        # Note that we do NOT need to check whether the plot has specified that a gene model should
        # be plotted.  That fact is already accounted for- if gene_detail_artist is not None, then we're good.
        if self.gene_detail_artist:
            transcript_features = annotation_data_source.get_genomic_features()
            self.gene_detail_artist.draw(transcript_features)

        """
        import numpy as np
        import pandas as pd
        y1 = 100*np.ones(1000) + np.random.randint(-5, 5, 1000)
        y2 = 50*np.ones(1000) + np.random.randint(-5, 5, 1000)
        yvals = np.vstack((y1,y2))
        cvg_data_df = pd.DataFrame(yvals)
        cvg_data_df.index = ['SampleA', 'SampleB']
        """
        #self.gene_detail_artist.draw()
        for a in self.generic_artists:
            a.draw()

    def save_figure(self, output_filepath):
        self.fig.savefig(output_filepath, bbox_inches='tight')