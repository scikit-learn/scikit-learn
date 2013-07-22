from pyearth._record import ForwardPassRecord, ForwardPassIteration, PruningPassRecord, PruningPassIteration
from pyearth._util import gcv
from nose.tools import assert_true, assert_equal, assert_list_equal

class TestForwardPassRecord(object):

    def __init__(self):
        #Create a record
        self.num_samples = 1000
        self.num_variables = 10
        self.penalty = 3.0
        self.sst = 100.0
        self.record = ForwardPassRecord(self.num_samples, self.num_variables, self.penalty, self.sst)
        self.record.append(ForwardPassIteration(0, 3, 3, 63.0, 3))
        self.record.append(ForwardPassIteration(0, 3, 14, 34.0, 5))
        self.record.append(ForwardPassIteration(3, 6, 12, 18.0, 7))
        self.mses = [self.sst,63.0,34.0,18.0]
        self.sizes = [1,3,5,7]

    def test_statistics(self):
        mses = [self.record.mse(i) for i in range(len(self.record))]
        mses_ = [self.mses[i] for i in range(len(self.record))]
        gcvs = [self.record.gcv(i) for i in range(len(self.record))]
        gcvs_ = [gcv(self.mses[i], self.sizes[i], self.num_samples, self.penalty) for i in range(len(self.record))]
        rsqs = [self.record.rsq(i) for i in range(len(self.record))]
        rsqs_ = [1 - (self.mses[i] / self.sst) for i in range(len(self.record))]
        grsqs = [self.record.grsq(i) for i in range(len(self.record))]
        grsqs_ = [1 - (self.record.gcv(i) / gcv(self.sst, 1, self.num_samples, self.penalty)) for i in range(len(self.record))]
        assert_list_equal(mses,mses_)
        assert_list_equal(gcvs,gcvs_)
        assert_list_equal(rsqs,rsqs_)
        assert_list_equal(grsqs,grsqs_)
        
        
class TestPruningPassRecord(object):

    def __init__(self):
        #Create a record
        self.num_samples = 1000
        self.num_variables = 10
        self.penalty = 3.0
        self.sst = 100.0
        self.record = PruningPassRecord(self.num_samples, self.num_variables, self.penalty, self.sst, 7, 18.0)
        self.record.append(PruningPassIteration(2, 6, 25.0))
        self.record.append(PruningPassIteration(1, 5, 34.0))
        self.record.append(PruningPassIteration(3, 4, 87.0))
        self.mses = [18.0,25.0,34.0,87.0]
        self.sizes = [7,6,5,4]

    def test_statistics(self):
        mses = [self.record.mse(i) for i in range(len(self.record))]
        mses_ = [self.mses[i] for i in range(len(self.record))]
        gcvs = [self.record.gcv(i) for i in range(len(self.record))]
        gcvs_ = [gcv(self.mses[i], self.sizes[i], self.num_samples, self.penalty) for i in range(len(self.record))]
        rsqs = [self.record.rsq(i) for i in range(len(self.record))]
        rsqs_ = [1 - (self.mses[i] / self.sst) for i in range(len(self.record))]
        grsqs = [self.record.grsq(i) for i in range(len(self.record))]
        grsqs_ = [1 - (self.record.gcv(i) / gcv(self.sst, 1, self.num_samples, self.penalty)) for i in range(len(self.record))]
        assert_list_equal(mses,mses_)
        assert_list_equal(gcvs,gcvs_)
        assert_list_equal(rsqs,rsqs_)
        assert_list_equal(grsqs,grsqs_)

if __name__ == '__main__':
    import nose
    nose.run(argv=[__file__, '-s', '-v'])
    
    