import exceptions


class BaseDataProvider(object):

    def __init__(self, *args, **kwargs):
        self.data_sources = {}

    def get_data_source(self, key):
        raise exceptions.NotImplementedError

    def add_data_source(self, key, data_source):
        '''
        Returns something.
        '''
        raise exceptions.NotImplementedError


class DefaultDataProvider(BaseDataProvider):

    def __init__(self, *args, **kwargs):
        super(DefaultDataProvider, self).__init__(*args, **kwargs)


    def get_data_source(self, key):
        if key in self.data_sources.keys():
            print 'check key %s' % key
            data_src = self.data_sources[key]
            return data_src
        else:
            raise DataSourceNotAvailableException('The data source %s was not set on the data provider' % key)


    def add_data_source(self, key, data_source):
        print 'add data src: %s, %s' % (key, data_source)
        self.data_sources[key] = data_source