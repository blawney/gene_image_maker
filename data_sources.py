import os
try:
    import HTSeq
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
        for i,fp in enumerate(self.filepathlist):
            coverage_mtx[i] = self._get_coverage_profile(fp, genomic_interval)
        return coverage_mtx

    @staticmethod
    def _get_coverage_profile(bam_file, window):
        """
        Reads through the BAM file and returns a numpy array giving the coverage for the transcript.
        """
        bam = HTSeq.BAM_Reader(bam_file)
        coverage = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
        try:
            for alnmt in bam[window]:
                if alnmt.aligned:
                    coverage[alnmt.iv]+=1
        except ValueError:
            print """
            Exception when reading the BAM file.
            This is common for two situations:
                1: There is no .bai file for your BAM file
                2: It is possible that your BAM file and GTF file do not have the same chromosome names.
             Check the chromosome names in the sam or index file and those in the GTF for agreement (and fix as necessary).
                 """
        # We now have coverage, which is a generator for tuples.
        # Each tuple has a GenomicInterval and an integer for the read-depth.  To eventually plot, we need to make this into a numpy array
        cvg_list = []
        it = coverage.steps() #an iterator
        try:
            step = it.next() #get the first object from the iterator so we can enter the while loop
            while step:
                if step[0].start<=window.end and step[0].end>=window.start:  #if step overlaps with window
                    if step[0].start<=window.start and step[0].end>=window.end: #if the step is longer than the window (unlikely, but possible)
                        cvg_list.extend(window.length*[step[1]])
                    elif step[0].start<=window.start:
                        cvg_list.extend(abs(step[0].end-window.start)*[step[1]])
                    elif step[0].end>=window.end:
                        cvg_list.extend(abs(window.end-step[0].start)*[step[1]])
                    else:
                        cvg_list.extend(step[0].length*[step[1]])
                step = it.next()
        except StopIteration: #when the generator is done, it throws an exception, which we catch and move on
            pass
        if len(cvg_list) == 0:
            print "Could not find any coverage data for the genomic region %s in BAM file %s " % (window, bam_file)
        return np.array(cvg_list)


class BEDCoverageDataSource(BEDDataSource):
    def get_coverage_data(self, genomic_interval):
        pass


class GTFGeneDetailDataSource(FileDataSource, AnnotationSourceMixin):

    ACCEPTED_FILE_EXTENSION = 'gtf'

    def __init__(self, filepath_or_filepathlist, **kwargs):
        super(GTFGeneDetailDataSource, self).__init__(filepath_or_filepathlist, **kwargs)

        # recall that the base class init method can handle multiple files, and stores them in a list
        self.gtf = HTSeq.GFF_Reader(self.filepathlist[0])
        self.feature_set = None

    def parse_gtf_features(self, transcript_id):
        self.feature_set = set()
        try:
            for feature in self.gtf:
                if feature.attr['transcript_id'] == transcript_id:
                    self.feature_set.add(feature)
        except Exception as ex:
            print 'Exception occurred while parsing features from GTF file: %s' % ex.message
            raise ex

    def get_window_interval(self, transcript_id, upstream_padding, downstream_padding):
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

    def get_genomic_feature_set(self):
        if not self.feature_set:
            self.parse_gtf_features(transcript_id)
        return self.feature_set


class BEDGeneDetailDataSource(BEDDataSource, AnnotationSourceMixin):

    def get_transcript_interval(self):
        pass

    def get_genomic_feature_set(self):
        pass