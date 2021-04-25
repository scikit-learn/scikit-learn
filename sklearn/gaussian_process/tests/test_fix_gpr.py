import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
import pytest
import warnings
import math


def test_case_0():
    noise=1e-5
    signal_var = 8.98576054e+05
    length_factor = np.array([5.91326520e+02, 1.32584051e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    #Setup Training Data
    x_train=np.array([[ 0., 0.],
    [ 1.54919334,-0.77459667],
    [-1.54919334, 0.        ],
    [ 0.         ,-1.54919334],
    [ 0.77459667  ,0.77459667],
    [-0.77459667  ,1.54919334]])
    y_train=np.array([[-2.14882017e-10],
    [-4.66975823e+00],
    [ 4.01823986e+00],
    [-1.30303674e+00],
    [-1.35760156e+00],
    [ 3.31215668e+00]])
    #fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[-1.93649167,-1.93649167],
    [ 1.93649167, -1.93649167],
    [-1.93649167,  1.93649167],
    [ 1.93649167,  1.93649167]])
    #Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2=np.sqrt(np.diagonal(cov))
    #Compare
    flag = True
    eps = 1e-5
    for i in range(0,len(std1)):
        if math.fabs(std1[i]-std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_1():
    noise = 3e-5
    signal_var = 8.98576054e+05
    length_factor = np.array([1.18859267e+02, 4.46284042e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[0.71388468, -2.74775681],
                        [0.29702229, 2.7254978],
                        [4.98924586, 3.53565024],
                        [-0.55163687, -2.08462431],
                        [0.52342054, 3.42620408],
                        [-2.31651775, -0.11808817]])
    y_train = np.array([[-2.14882017e-10],
                        [-4.66975823e+00],
                        [4.01823986e+00],
                        [-1.30303674e+00],
                        [-1.35760156e+00],
                        [3.31215668e+00]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[-1.93649167, -1.93649167],
                       [1.93649167, -1.93649167],
                       [-1.93649167, 1.93649167],
                       [1.93649167, 1.93649167]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_2():
    noise = 1e-5
    signal_var = 8.68797886e+05
    length_factor = np.array([2.13041806e+02, 3.35316558e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[2.62922842, -1.54594557],
                        [0.55296607, 0.80854811],
                        [-3.66069924, 0.52653935],
                        [2.41656578, -3.24047799],
                        [0.16370801, -4.96687776],
                        [-1.42230462, 3.61330179]])
    y_train = np.array([[-3.80432314],
                        [-3.1269489],
                        [3.91081338],
                        [2.32182974],
                        [-0.53099622],
                        [-3.81636951]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[0.54811996, -4.87313084],
                       [-3.85052863, -0.56658422],
                       [-4.29348578, 2.77510038],
                       [-3.49051566, 4.36136012]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_3():
    noise = 1e-5
    signal_var = 5.4895927e+05
    length_factor = np.array([1.85666783e+02, 3.25192115e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[0.80201264, -2.65846114],
                        [-1.8659259, -1.90181379],
                        [-2.42783625, -4.82627679],
                        [-4.75016133, 4.07930783],
                        [-4.70858611, 3.17706985],
                        [0.30492441, -0.64503186]])
    y_train = np.array([[-2.96362039],
                        [-4.07664738],
                        [0.05367682],
                        [1.3998386],
                        [-0.79427188],
                        [-2.60490784]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[2.11263668, -4.26520676],
                       [-4.88826261, -3.85332622],
                       [-1.69552295, 1.78326764],
                       [4.46614447, -1.1270797]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_4():
    noise = 1e-5
    signal_var = 5.47311595e+05
    length_factor = np.array([3.0556732e+02,  4.1673091e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[-4.31631744, 2.84307697],
                        [1.37280169, -1.27588091],
                        [0.7600307, -1.84746334],
                        [4.72788714, 1.09226148],
                        [-1.84420944, 1.50509334],
                        [-1.38126041, 4.2586101]])
    y_train = np.array([[-1.37203372],
                        [2.11269714],
                        [-2.20488749],
                        [4.78006842],
                        [-2.41567304],
                        [-1.05789616]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[-2.98131042, -2.83254693],
                       [-4.91164618, 0.33522311],
                       [0.02995569, 3.43234924],
                       [-2.73569117, -2.67387784]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_5():
    noise = 1e-5
    signal_var = 5.13464393e+05
    length_factor = np.array([2.512199e+02, 2.52732124e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[1.70729543, -1.33570758],
                        [-0.35153569, 3.41125631],
                        [0.17698545, -4.6952475],
                        [-0.97741929, 1.37824434],
                        [-2.34283391, 1.84257246],
                        [2.32144872, -0.60873786]])
    y_train = np.array([[1.63078094],
                        [-0.28267535],
                        [3.8005446],
                        [4.86606608],
                        [-4.30706923],
                        [-2.96144082]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[4.61189266, 1.11900243],
                       [2.2858796, -1.85619925],
                       [2.82782816, -0.66851802],
                       [-0.70499443, -1.6715043]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_6():
    noise = 1e-5
    signal_var = 7.00283167e+05
    length_factor = np.array([1.69365166e+02, 3.82141751e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[2.62922842, -1.54594557],
                        [0.55296607, 0.80854811],
                        [-3.66069924, 0.52653935],
                        [2.41656578, -3.24047799],
                        [0.16370801, -4.96687776],
                        [-1.42230462, 3.61330179]])
    y_train = np.array([[-4.57938803],
                        [3.94670756],
                        [-4.88687977],
                        [-4.4543016],
                        [4.51344949],
                        [-1.38553635]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[-1.52972462, -4.6938267],
                       [0.26004387, -3.45781447],
                       [1.5570652, 4.81672381],
                       [1.55777735, -1.89621684]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_7():
    noise = 1e-5
    signal_var = 7.80107872e+05
    length_factor = np.array([1.73880921e+02, 1.90185775e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[-2.00705433, 3.68402189],
                        [-3.17284422, 3.08148533],
                        [3.36232363, 1.33984468],
                        [-0.47052089, 0.38719207],
                        [0.68414026, 4.84590746],
                        [2.59202292, -3.51154102]])
    y_train = np.array([[3.22760786],
                        [-2.2403138],
                        [-3.45937986],
                        [4.05560292],
                        [2.92191944],
                        [-3.30199415]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[-0.19757781, 2.34556484],
                       [-2.9008185, 2.35625496],
                       [3.62826885, -1.56896901],
                       [-2.28506659, 3.6865107]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_8():
    noise = 1e-5
    signal_var = 3.22261793e+05
    length_factor = np.array([3.4518944e+02, 1.58767155e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[2.75140186, 4.5832697],
                        [1.6915996, 3.6972948],
                        [0.2949518, -3.95877262],
                        [-4.45227312, -4.94865286],
                        [2.63342597, 3.59550855],
                        [3.30183712, -1.83339388]])
    y_train = np.array([[0.47113625],
                        [4.38526708],
                        [2.13907661],
                        [-3.71215321],
                        [-2.57122424],
                        [-1.72059563]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[-2.50589975, 1.16609763],
                       [1.95174744, 4.87319831],
                       [-3.54742382, 1.64943713],
                       [4.75849572, -3.31399285]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True


def test_case_9():
    noise = 1e-5
    signal_var = 7.69449168e+05
    length_factor = np.array([4.30117815e+02, 4.57265035e+03])
    kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
    m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
    # Setup Training Data
    x_train = np.array([[1.1836609, -2.8671946],
                        [-0.567545, 1.4199508],
                        [-3.72871151, -3.99840598],
                        [4.5380178, 0.47819713],
                        [1.56294593, 2.15874747],
                        [3.52000758, -4.64742481]])
    y_train = np.array([[-1.75511956],
                        [-2.13036887],
                        [-1.62971815],
                        [4.07889753],
                        [-1.18763414],
                        [4.30117815]])
    # fit. We must do this to register teh GPR as "fitted"
    m.fit(x_train, y_train)
    theta = np.exp(m.kernel_.theta)
    x_test = np.array([[-3.86479357, 3.03253336],
                       [-2.82328255, 3.03209198],
                       [0.00277192, 2.47264377],
                       [-3.45211314, -0.21028548]])
    # Predict the std_dev in two ways
    pred1, std1 = m.predict(x_test, return_std=True)
    pred2, cov = m.predict(x_test, return_cov=True)
    std2 = np.sqrt(np.diagonal(cov))
    # Compare
    flag = True
    eps = 1e-5
    for i in range(0, len(std1)):
        if math.fabs(std1[i] - std2[i]) >= eps:
            flag = False
    assert flag == True