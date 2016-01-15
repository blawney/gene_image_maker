import os
try:
    import HTSeq
    import pysam
    import numpy as np
    import pandas as pd
except ImportError as ex:
    print "Failed when importing dependencies: %s"  % ex.message
    raise ex


class TranscriptNotFound(Exception):
    pass


class BaseDataSource(object):
    pass


class AnnotationSourceMixin(object):

    def get_transcript_interval(self):
        '''
        Returns a HTSeq.GenomicInterval which represents the coordinates of the transcript of interest
        '''
        raise NotImplementedError('Not implemented in derived class!')

    def get_genomic_feature_set(self):
        '''
        returns a set of HTSeq.GenomicFeature's
        '''
        raise NotImplementedError('Not implemented in derived class!')


class CoverageDataSourceMixin(object):

    def get_coverage_data(self, genomic_interval):
        raise NotImplementedError('Not implemented in derived class!')


class FileDataSource(BaseDataSource):

    def __init__(self, filepath_or_filepathlist, many=False):
        self.multiple_files = many

        # put into a list for consistent handling of single or multiple files
        if not many:
            filepathlist = [filepath_or_filepathlist,]
        else:
            filepathlist = filepath_or_filepathlist

        try:
            print 'checking file: %s' % filepathlist
            for f in filepathlist:
                os.stat(f)
            self.filepathlist = filepathlist
        except OSError as ex:
            print 'Could not stat file: %s' % f
            raise ex


class BEDDataSource(FileDataSource):
    ACCEPTED_FILE_EXTENSION = 'bed'


class BAMCoverageDataSource(FileDataSource, CoverageDataSourceMixin):
    ACCEPTED_FILE_EXTENSION = 'bam'

    def get_coverage_data(self, genomic_interval):
        coverage_mtx = np.empty((len(self.filepathlist), genomic_interval.length))
        for i, fp in enumerate(self.filepathlist):
            print fp
            coverage_mtx[i] = self._get_coverage_profile(fp, genomic_interval)
        sample_names = ['.'.join(x.split('.')[:-1]) for x in self.filepathlist]
        df = pd.DataFrame(coverage_mtx)
        df.index = sample_names
        df.columns = np.arange(genomic_interval.start, genomic_interval.end)
        return df

    @staticmethod
    def _get_coverage_profile(bam_file, window):
        """
        Reads through the BAM file and returns a numpy array giving the coverage for the transcript.
        """
        bam = pysam.AlignmentFile(bam_file, 'rb')
        return np.array(bam.count_coverage(window.chrom, window.start, window.end)).sum(axis=0)


class BEDCoverageDataSource(BEDDataSource):
    def get_coverage_data(self, genomic_interval):
        pass


class GTFGeneDetailDataSource(FileDataSource, AnnotationSourceMixin):

    ACCEPTED_FILE_EXTENSION = 'gtf'

    def __init__(self, filepath_or_filepathlist, **kwargs):
        print 'creating GTF data source'
        super(GTFGeneDetailDataSource, self).__init__(filepath_or_filepathlist, **kwargs)

        # recall that the base class init method can handle multiple files, and stores them in a list
        self.gtf = HTSeq.GFF_Reader(self.filepathlist[0])
        self.feature_set = None

    def parse_gtf_features(self, transcript_id):
        self.feature_set = set()
        try:
            for feature in self.gtf:
                try:
                    if feature.attr['transcript_id'] == transcript_id:
                        self.feature_set.add(feature)
                except:
                    # if a line in the GTF file doesn't have a transcript field (e.g. a line denoting the gene)
                    # then the try block above will throw a KeyError on that dictionary lookup.
                    # Silently catch and move on
                    pass
        except Exception as ex:
            print 'Exception occurred while parsing features from GTF file: %s' % ex.message
            raise ex

    def get_window_interval(self, transcript_id, upstream_padding, downstream_padding):
        print 'search for %s' % transcript_id
        if not self.feature_set:
            self.parse_gtf_features(transcript_id)
        transcript_found = False
        for f in self.feature_set:
            if f.type == "transcript":
                if f.iv.strand == '-':
                    new_end = f.iv.end + upstream_padding
                    new_start = f.iv.start - downstream_padding
                else:
                    new_start = f.iv.start - upstream_padding
                    new_end = f.iv.end + downstream_padding
                return HTSeq.GenomicInterval(f.iv.chrom, new_start, new_end, f.iv.strand)

        raise TranscriptNotFound("Could not locate a GenomicFeature of type "
                                 "'transcript', so could not create a GenomicInterval.")

    def get_genomic_features(self):
        if not self.feature_set:
            self.parse_gtf_features(transcript_id)
        return list(self.feature_set)


class BEDGeneDetailDataSource(BEDDataSource, AnnotationSourceMixin):

    def get_transcript_interval(self):
        pass

    def get_genomic_feature_set(self):
        pass