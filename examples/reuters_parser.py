# Author: Eustache Diemert <eustache@diemert.fr>
# License: BSD 3 clause

import os.path
import fnmatch
import re
import sgmllib


class ReutersTopics():
    """Utility class to read official topic names from the relevant file."""
    TOPICS_FILENAME = 'all-topics-strings.lc.txt'
    def __init__(self, topics_path):
        self.topics_ = open(topics_path).read().split('\n')
        self.topics_ = dict([ (self.topics_[i],i) for i in range(len(self.topics_)) ])

    def topic_ids(self):
        return self.topics_.values()



class ReutersParser(sgmllib.SGMLParser):
    """Utility class to parse a SGML file and yield documents one at a time."""
    def __init__(self, verbose=0):
        sgmllib.SGMLParser.__init__(self, verbose)
        self._reset()

    def _reset(self):
        self.in_title = 0
        self.in_body = 0
        self.in_topics = 0
        self.in_topic_d = 0
        self.title = ""
        self.body = ""
        self.topics = []
        self.topic_d = ""

    def parse(self, fd):
        self.docs = []
        for chunk in fd:
            self.feed(chunk)
            for doc in self.docs:
                yield doc
            self.docs = []
        self.close()

    def handle_data(self, data):
        if self.in_body:
            self.body += data
        elif self.in_title:
            self.title += data
        elif self.in_topic_d:
            self.topic_d += data

    def start_reuters(self, attributes):
        pass

    def end_reuters(self):
        self.body = re.sub(r'\s+', r' ', self.body)
        self.docs.append({'title':self.title,
                          'body':self.body,
                          'topics':self.topics})
        self._reset()

    def start_title(self, attributes):
        self.in_title = 1

    def end_title(self):
        self.in_title = 0

    def start_body(self, attributes):
        self.in_body = 1

    def end_body(self):
        self.in_body = 0

    def start_topics(self, attributes):
        self.in_topics = 1

    def end_topics(self):
        self.in_topics = 0

    def start_d(self, attributes):
        self.in_topic_d = 1

    def end_d(self):
        self.in_topic_d = 0
        self.topics.append(self.topic_d)
        self.topic_d = ""



class ReutersStreamReader():
    """Reads documents form the directory where the Reuters dataset has been 
    uncompressed. Documents are represented as dictionaries with 'body' (str), 
    'title' (str), 'topics' (list(str)) keys. 
    """
    def __init__(self, data_path):
        self.data_path = data_path
        self.topics = ReutersTopics(os.path.join(data_path,
                                                 ReutersTopics.TOPICS_FILENAME))
        self.classes = self.topics.topic_ids()

    def iterdocs(self):
        for root, _dirnames, filenames in os.walk(self.data_path):
            for filename in fnmatch.filter(filenames, '*.sgm'):
                path = os.path.join(root, filename)
                parser = ReutersParser()
                for doc in parser.parse(open(path)):
                    yield doc



###############################################################################

if __name__ == '__main__':
    """Test streamer by printing the first document"""
    import sys
    path = sys.argv[1] # path to *.sgm files 
    data_streamer = ReutersStreamReader(path, 10)
    for doc in data_streamer.iterdocs():
        print doc
        break
