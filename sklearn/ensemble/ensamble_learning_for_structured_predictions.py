from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
import scipy.spatial.distance as distance
import math
from numpy import array
import heapq

class EnsambleLearningforStructuredPredictions():

    def __init__(self, beta=0.5, L=2):
        self.BETA = beta
        self.EXPERTS = ['svm', 'GradientBoostingClassifier', 'linear-svc']
        self.NUMBER_OF_LAST_Ws = 1
        self.W_T = []
        self.number_of_substructures_L = L
        self.P_EXPERTS = len(self.EXPERTS)
        self.trained_models = []
        self.sub_Xs = []
        self.expert_path = []
        self.expert_weights = []
        self.params = {'beta':beta, 'L':L}


    def fit(self,X,Y):
        self.initW()
        self.init_experts(X,Y)
        for t in range(0,len(X)):
            self.create_W(X[t],Y[t]) # maybe send 'sub_X[t]' instead of X[t]
            # if t % 1000 == 0:
                # print "finished round {}".format(t)
        self.get_expert_path()

    def get_expert_path(self):
        W = self.W_T[-1]
        trained_models_index = 0
        for row in W:
            maxIndx = row.index(max(row))
            self.expert_path.append(self.trained_models[trained_models_index][maxIndx])
            self.expert_weights.append(row[maxIndx])
            trained_models_index += 1

    def get_params(self):
        return self.params


    def set_params(self,**params):
        if not params:
            return self

        for key, value in params.iteritems():
            self.params[key] = value

        return self

    def predict(self, X):
        interval = len(X[0]) / self.number_of_substructures_L

        final_predicted_Y = []
        for row in X:
            current_Y = []
            start_index = 0
            for expert in range(self.number_of_substructures_L):
                sub_x = row[start_index: start_index + interval]
                start_index += interval
                y = self.expert_path[expert].predict(sub_x)
                current_Y.append(y[0])
            best_y_per_currect_row = self.get_optimal_prediction(current_Y)
            final_predicted_Y.append(best_y_per_currect_row)
        return final_predicted_Y

    def predict_proba(self, X):
        interval = len(X[0]) / self.number_of_substructures_L

        final_predicted_Y = []
        for row in X:
            current_Y = []
            start_index = 0
            for expert in range(self.number_of_substructures_L):
                sub_x = row[start_index: start_index + interval]
                start_index += interval
                try:
                    y = self.expert_path[expert].predict_proba(sub_x)
                except: # in case of the expert is not implementing 'predict_proba'
                    continue
                current_Y.append(y[0])
            best_y_per_currect_row = self.get_optimal_proba(current_Y)
            final_predicted_Y.append(best_y_per_currect_row)
        return final_predicted_Y

    def get_optimal_proba(self,Y):
        best_Y = 0

        for y in array(Y).flat:
            if y > best_Y:
                best_Y = y
        return best_Y

    def get_optimal_prediction(self, Y):
        best_y = -1
        W = self.W_T[-1]
        hist = dict((x,Y.count(x)) for x in Y)
        first_max, second_max = self.get_first_and_second_max(hist)

        if second_max == None:
            return first_max

        if(hist[first_max] > hist[second_max]):
            return first_max
        else:
            # check weights
            sum_first = 0
            sum_second = 0
            for i in range(len(Y)):
                curr_y = Y[i]
                if curr_y == first_max:
                    sum_first += self.expert_weights[i]
                elif curr_y == second_max:
                    sum_second += self.expert_weights[i]
            if sum_first >= sum_second:
                return first_max
            else:
                return second_max

    def get_first_and_second_max(self,hist):
        lst = heapq.nlargest(2, hist,key=lambda k : hist[k])
        if len(lst) > 1:
            return lst[0], lst[1]
        else:
            return lst[0], None


    def create_W(self,X,Y):
        W_t_plus_one = []
        for k in range(0,self.number_of_substructures_L):
            sigma_p = self.calc_sigma_p(self.sub_Xs[k],Y,k)
            w_t_row = []
            for j in range(0, self.P_EXPERTS):
                mone = self.calc_mone(k,j,Y)
                result = mone / sigma_p
                w_t_row.append(result)

            W_t_plus_one.append(w_t_row)
        self.W_T.append(W_t_plus_one)


    def calc_mone(self,k,j,Y):
        W_t = self.W_T[-1]
        W_t_value = W_t[k][j]
        predicted_y = self.trained_models[k][j].predict(self.sub_Xs[k][j])

        hamming_dist = distance.hamming(predicted_y[0], Y)
        return W_t_value * math.pow(self.BETA,hamming_dist)


    def initW(self):
        val = float(1) / float(self.P_EXPERTS)
        W = [[val for x in range(0, self.P_EXPERTS)] for x in range(0,self.number_of_substructures_L)]
        self.W_T.append(W)


    def calc_sigma_p(self,X,Y,k):
        sum = 0
        for j in range(self.P_EXPERTS):
            W_t = self.W_T[-1]
            W_t_value = W_t[k][j]
            predicted_y = self.trained_models[k][j].predict(X[j])

            hamming_dist = distance.hamming(predicted_y, Y)
            sum += W_t_value * math.pow(self.BETA,hamming_dist)
            # print "hamming_dist =[%s]   W_t_value = [%s]   math.pow(self.BETA,hamming_dist) = [%s]" % (hamming_dist, W_t_value,math.pow(self.BETA,hamming_dist))
            # print "sum = [%s]" % sum
        return  sum


    def init_experts(self,X,Y):

        interval = len(X[0]) / self.number_of_substructures_L
        start_index = 0
        for substruct in range(self.number_of_substructures_L):
            sub_X = [[x for x in x[start_index:start_index+interval]] for x in X]
            self.sub_Xs.append(sub_X)
            start_index += interval

            model_row = []
            for clf_name in self.EXPERTS:
                clf = self.get_classifier(clf_name)
                clf.fit(sub_X, Y)
                model_row.append(clf)

            self.trained_models.append(model_row)


    def get_classifier(self, name):
        if name == 'svm':
            svm = SVC(probability=True)
            # print "SVM"
            # print(svm)
            return svm
        elif name == 'random_forest':
            rf = RandomForestClassifier()
            # print rf
            return rf
        elif name == 'linear-svc':
            lsvc = LinearSVC()
            # print "LinearSVC"
            # print lsvc
            return lsvc
        elif name == "GradientBoostingClassifier":
            g = GradientBoostingClassifier()
            # print "GradientBoostingClassifier"
            # print GradientBoostingClassifier
            return g