from pip._vendor.six.moves.urllib import parse as urllib_parse


class Index(object):
    def __init__(self, url):
        self.url = url
        self.netloc = urllib_parse.urlsplit(url).netloc
        self.simple_url = self.url_to_path('simple')
        self.pypi_url = self.url_to_path('pypi')

    def url_to_path(self, path):
        return urllib_parse.urljoin(self.url, path)


PyPI = Index('https://pypi.org/')
