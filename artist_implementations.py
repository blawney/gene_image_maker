class PlotFeatureArtist(object):

    def __init__(self, ax):
        self.ax = ax

    def draw(self):
        raise NotImplementedError('Implement this method in the derived class!')


class GeneDetailArtist(PlotFeatureArtist):
    def draw(self):
        print 'in gene detail draw'


class CoverageProfileArtist(PlotFeatureArtist):
    def draw(self):
        ax = self.ax
        print 'in coverage profile draw'
        import numpy as np
        ax.scatter(np.random.random(10), np.random.random(10))
        ax.set_xlabel('Xaxis')


class FormattingArtist(PlotFeatureArtist):
    def draw(self):
        print 'in other draw'
