# Here we specify the concrete classes used.

# class for exposing the data to the plotting methods
# should derive from data_provider.BaseDataProvider
# Unless special extras are needed this probably doesn't need to change
DATA_PROVIDER = 'data_provider.DefaultDataProvider'

# Specifies how to draw the coverage plot
COVERAGE_PROFILE_ARTIST = 'artist_implementations.CoverageProfileArtist'

# Specifies how to draw the gene image/model
GENE_DETAIL_ARTIST = 'artist_implementations.GeneDetailArtist'

# These are implementations of the PlotFeatureArtist "interface" that do not depend on data
# (e.g. other styling, such as ticks, fonts, etc)
OTHER_ARTISTS = ['artist_implementations.FormattingArtist']

# Provides the coverage data (e.g. a count matrix)
COVERAGE_DATA_SOURCES = ['data_sources.BAMCoverageDataSource',
                         'data_sources.BEDCoverageDataSource']

# Provides the data on the gene (e.g. the exon boundaries, UTRs, etc.)
# so the gene model can be drawn
GENE_DETAIL_DATA_SOURCES = ['data_sources.GTFGeneDetailDataSource',
                            'data_sources.BEDGeneDetailDataSource']

# boolean to specify whether the gene model is plotted with the coverage profile
INCLUDE_GENE_MODEL = True

# CONSTs for consistent reference (doesn't matter what the actual string values are as long as they are unique)
ANNOTATION_SRC = 'ANN'
COVERAGE_SRC = 'CVG'

FIGURE_SIZE = (20,20)

