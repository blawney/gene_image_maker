import numpy as np
import matplotlib as mpl

class PlotFeatureArtist(object):

    def __init__(self, ax):
        self.ax = ax

    def draw(self):
        raise NotImplementedError('Implement this method in the derived class!')


class GeneDetailArtist(PlotFeatureArtist):

    # some constants for basic customization
    INCLUDED_FEATURES = ['utr', 'exon'] # this also establishes a precendence for plotting.
    FEATURE_SIZES = {'exon': 0.2, 'utr': 0.1, 'intron':0.05}
    GENE_MODEL_COLOR = 'gray'
    GENE_MODEL_OFFSET = 0.05 # percentage of total plot height that gene model is below horizontal axis
    DIRECTION_ARROW_STALK_THICKNESS = 0.025
    DIRECTION_ARROW_HEAD_HEIGHT = 0.1

    @staticmethod
    def _prepare_features(transcript_features):
        transcript_features.sort(key=lambda f: f.iv.start)
        return filter(lambda f: f.type.lower() in GeneDetailArtist.INCLUDED_FEATURES, transcript_features)

    @staticmethod
    def _get_breakpoints(feature_list):
        points = set()
        for f in feature_list:
            points.add(f.iv.start)
            points.add(f.iv.end)
        return sorted(list(points))

    @staticmethod
    def _overlap(feature_range, other_range):
        return (feature_range[1] > other_range[0]) and (feature_range[0] < other_range[1])

    @staticmethod
    def _assign_gene_elements(features, breakpoints):
        range_dict = {}
        for i in range(len(breakpoints)-1):
            range_dict[(breakpoints[i], breakpoints[i+1])] = []
        for f in features:
            found = False
            feature_range = (f.iv.start, f.iv.end)
            for range_key in range_dict.keys():
                if GeneDetailArtist._overlap(feature_range, range_key):
                    range_dict[range_key].append(f.type.lower())
        return range_dict

    def draw(self, transcript_features):
        print 'in gene detail draw'
        filtered_features = GeneDetailArtist._prepare_features(transcript_features)
        breakpoints = GeneDetailArtist._get_breakpoints(filtered_features)
        range_and_feature_map = GeneDetailArtist._assign_gene_elements(filtered_features, breakpoints)

        # get strand from the features
        strand = filtered_features[0].iv.strand

        ax = self.ax
        ph = ax.get_ylim()[1]
        gene_model_top = -GeneDetailArtist.GENE_MODEL_OFFSET * ph
        gene_model_height = ph *np.max(GeneDetailArtist.FEATURE_SIZES.values())
        gene_model_midline = gene_model_top - 0.5 * gene_model_height

        # draw a box across the whole transcript.  This will eventually be covered by exons, etc so this is
        # effectively drawing introns
        h = GeneDetailArtist.FEATURE_SIZES['intron'] * ph
        transcript_patches = []
        transcript_patches.append(
            ax.add_patch(
                mpl.patches.Rectangle(
                    (breakpoints[0], gene_model_midline - 0.5 * h),
                    breakpoints[-1]-breakpoints[0],
                    h,
                    facecolor='gray',
                    linewidth=0
                )
            )
        )

        print range_and_feature_map
        for feature_range, contents in range_and_feature_map.items():
            for f in GeneDetailArtist.INCLUDED_FEATURES:
                if f in contents:
                    print 'plot feature: %s' % f
                    h = GeneDetailArtist.FEATURE_SIZES[f] * ph
                    transcript_patches.append(
                        ax.add_patch(
                            mpl.patches.Rectangle(
                                (feature_range[0], gene_model_midline - 0.5 * h),
                                feature_range[1]-feature_range[0],
                                h,
                                edgecolor = 'black',
                                facecolor='gray',
                                linewidth=0
                            )
                        )
                    )
                    break


        for l in transcript_patches:
            l.set_clip_on(False)

        aspect_ratio = ax.get_aspect()
        arrow_base_x = breakpoints[0]
        arrow_base_y = gene_model_midline
        t = GeneDetailArtist.DIRECTION_ARROW_STALK_THICKNESS * ph
        H = GeneDetailArtist.DIRECTION_ARROW_HEAD_HEIGHT * ph
        ex_h = GeneDetailArtist.FEATURE_SIZES['exon'] * ph
        l = 1.5*t + 0.5*H + 0.5*ex_h
        rect1 = mpl.patches.Rectangle(
            (arrow_base_x, arrow_base_y-l),
            t*aspect_ratio,
            l,
            facecolor='gray',
            linewidth=0
        )
        rect2 = mpl.patches.Rectangle(
            (arrow_base_x, arrow_base_y-l),
            l*aspect_ratio,
            t,
            facecolor='gray',
            linewidth=0,
        )
        rect1.set_clip_on = True

        xs = arrow_base_x + l*aspect_ratio
        ys = arrow_base_y - l + 0.5*t - 0.5*H
        tri_h = 0.5 * np.sqrt(3) * H * aspect_ratio
        pt_array = np.empty((4,2))
        pt_array[0] = [xs, ys]
        pt_array[3] = [xs, ys]
        pt_array[1] = [xs + tri_h, ys + 0.5*H]
        pt_array[2] = [xs, ys + H]
        tri = mpl.patches.Polygon(pt_array,facecolor='gray', edgecolor='gray')

        arrow_group = mpl.collections.PatchCollection([rect1, rect2, tri], match_original = True)
        arrow_group.set_clip_on(False)

        # this uses a bit of linear algebra to flip the arrow across the vertical axis so we don't have to create
        # two drawing cases with different orientations, etc.
        if strand == '-':
            mirror_transform = mpl.transforms.Affine2D().from_values(-1,0,0,1,0,0).translate(2*arrow_base_x + breakpoints[-1] - breakpoints[0],0) + ax.transData
            arrow_group.set_transform(mirror_transform)
        ax.add_collection(arrow_group)



class CoverageProfileArtist(PlotFeatureArtist):
    def draw(self, df):
        ax = self.ax
        print 'in coverage profile draw'
        colors=['#97A861', "#9D6188"]
        #df = 5*df
        # plot the coverage profile:
        for k in range(df.shape[0]):
            fill = ax.fill_between(df.columns, df.iloc[k,:], 0)
            fill.set_facecolor(colors[k])
            fill.set_edgecolor(colors[k])
            fill.set_alpha(0.5)
            fill.set_linewidth(0.5)

        # for aspect ratio, limits, etc.
        figure_w_to_h = 5
        max_cvg = df.max().max()
        aspect_ratio = (df.columns.max() - df.columns.min())/float(max_cvg * figure_w_to_h)
        print aspect_ratio
        print type(aspect_ratio)
        #aspect_ratio = 1.6
        #print aspect_ratio
        ax.set_aspect(aspect_ratio)


class FormattingArtist(PlotFeatureArtist):

    DEFAULT_FONT = {'family':'serif', 'size':14, 'weight':'light'}
    MARKER_INC = 100

    def draw(self):
        print 'in other draw'
        ax = self.ax
        #ax.set_yticks([max_depth-1]) #ticks need to be less than the actual data max for them to show up
        ax.set_xticks([])
        ax.get_yaxis().tick_left()
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_position('zero')
        #ax.spines['left'].set_bounds(0, 150)
        ax.spines['right'].set_visible(False)
        x_axis_font = mpl.font_manager.FontProperties(**self.DEFAULT_FONT)
        y_axis_font = mpl.font_manager.FontProperties(**self.DEFAULT_FONT)
        [t.set_fontproperties(x_axis_font) for t in ax.get_xticklabels()]
        [t.set_fontproperties(x_axis_font) for t in ax.get_yticklabels()]

        ymin, ymax = ax.get_ylim()
        #ymax_new = 1.05 * ymax
        #ax.set_ylim((ymin, ymax_new))

        ax.set_yticks([ymax/2, ymax])

        # temp- REMOVE
        #ymax = -60

        """
        xmin, xmax = ax.get_xlim()
        axis_range_x = xmax - xmin
        marker_fraction = 0.25 # approx how long should the length indicator be, as a fraction of total x-axis
        marker_length = marker_fraction * axis_range_x

        # this makes it some 'usual' number so we don't have a scale marker with 823bp- make it like 800bp
        marker_length = int(marker_length/FormattingArtist.MARKER_INC)*FormattingArtist.MARKER_INC

        marker_left = xmin
        marker_right = xmin + marker_length
        marker_location = -0.0 * ymax # scales with the size of the figure

        length_marker_lines = []
        #length_marker_lines.append( ax.add_line(mpl.lines.Line2D((marker_left, marker_right),(marker_location, marker_location))))

        vert_marker_height = 0.05 * ymax # percent of total plot height
        length_marker_lines.append( ax.add_line(mpl.lines.Line2D((marker_left, marker_left),(marker_location-vert_marker_height, marker_location))))
        length_marker_lines.append( ax.add_line(mpl.lines.Line2D((marker_right, marker_right),(marker_location-vert_marker_height, marker_location))))

        #length_marker_lines.append(ax.add_line(mpl.lines.Line2D((marker_right, marker_right),(marker_location*(1-vert_marker_height), marker_location*(1+vert_marker_height)))))

        midpoint = 0.5 * (marker_left + marker_right)
        #ax.text(midpoint, marker_location, str(int(marker_length)) + ' bp', backgroundcolor=[1,1,1,0], fontproperties=x_axis_font, horizontalalignment='center', verticalalignment='top')


        ax.annotate(
            '', xy=(xmin, -20), xycoords='data',
            xytext=(xmin + marker_length, -20), textcoords='data', annotation_clip = False,
            arrowprops={'arrowstyle': '<->'})

        #ax.text(xmin, 0, 'ENSTXYZ', backgroundcolor='white', fontproperties=x_axis_font, horizontalalignment='left', verticalalignment='top')


        for l in length_marker_lines:
            l.set_color('k')
            l.set_clip_on(False)

        #ax.annotate('', xy=(0.8, -0.1), xycoords='axes fraction', xytext=(1, -0.1),
        #    arrowprops=dict(arrowstyle="<->", color='b'))

        """
