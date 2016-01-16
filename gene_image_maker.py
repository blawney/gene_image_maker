import sys
import utils
from artist import CoveragePlotManager


def plot(**kwargs):

    try:
        manager = CoveragePlotManager()

        # inject the proper data sources and pass that provider to the plot manager
        data_provider = utils.configure_data_provider(**kwargs)
        manager.set_data_provider(data_provider)
        print 'add artists'
        # Add artists to the manager
        cvg_artist_clz, gene_detail_artist_clz, additional_artists = utils.get_artists(**kwargs)
        print 'pt1'
        manager.add_artists(cvg_artist_clz, gene_detail_artist_clz, *additional_artists)

        manager.draw(**kwargs)
        manager.save_figure(kwargs['output_path'])

    except Exception as ex:
        print ex.message
        sys.exit(1)


def main():
    cl_args = utils.CLParser.parse_cl_args()
    plot(**cl_args)


if __name__ == "__main__":
    main()

