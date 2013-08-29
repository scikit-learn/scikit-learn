import numpy as np
import collections
from scipy import stats


def create_random_codebook(n_classes, code_size, random_state, estimator=None):
    code_book = random_state.random_sample((n_classes, code_size))
    code_book[code_book > 0.5] = 1

    if hasattr(estimator, "decision_function"):
        code_book[code_book != 1] = -1
    else:
        code_book[code_book != 1] = 0

    return code_book


def create_decoc_codebook(n_classes, X, y, random_state):
    code_book = np.ones(shape=(n_classes, n_classes-1), dtype=np.int)

    class_from = 0
    class_to = n_classes-1

    sigma = calculate_sigma(X)
    mutualInformationMatrix = calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = calculate_class_frequencies(y)

    parse_decoc(np.unique(y), code_book, class_from,
                class_to, random_state,
                mutualInformationMatrix, classFrequencies)

    return code_book


def calculate_class_frequencies(y):
    classes = np.unique(y)
    n_classes = len(classes)
    classFrequencies = np.zeros((n_classes,), dtype=np.int)

    for clsIdx in range(n_classes):
        cls = classes[clsIdx]
        y_selected = y[y == cls]
        classFrequencies[clsIdx] = len(y_selected)

    return classFrequencies


def calculate_mutual_information_matrix(X, y, sigma):
    classes = np.unique(y)
    n_classes = len(classes)
    mutualInformationMatrix = np.zeros((n_classes, n_classes))

    for clsLeftIdx in range(n_classes):
        for clsRightIdx in range(clsLeftIdx + 1):
            class_left = classes[clsLeftIdx]
            class_right = classes[clsRightIdx]

            X_left = X[y == class_left]
            X_right = X[y == class_right]

            mutInfAt = calculate_mutual_information_atom(X_left, X_right, sigma)
            mutualInformationMatrix[clsLeftIdx, clsRightIdx] = mutInfAt

    for clsLeftIdx in range(n_classes):
        for clsRightIdx in range(clsLeftIdx+1, n_classes):
            mt_atom = mutualInformationMatrix[clsRightIdx, clsLeftIdx]
            mutualInformationMatrix[clsLeftIdx, clsRightIdx] = mt_atom

    return mutualInformationMatrix


def calculate_mutual_information_atom(X_left, X_right, sigma):
    mutInfAt = sum([calculate_parzen_estimate(x_l, x_r, sigma)
                    for x_l in X_left
                    for x_r in X_right])
    return mutInfAt


def calculate_parzen_estimate(x, y, sigma):
    if x.shape != y.shape:
        raise ValueError('X and Y vectors should have the same length!')

    estimate = np.sqrt(np.matrix(x)*np.matrix(y).T)
    parzen_estimate = stats.norm.pdf(estimate, sigma)
    return parzen_estimate


def parse_decoc(current_classes, code_book, class_from,
                class_to, random_state,
                mutualInformationMatrix, classFrequencies):

    if class_from >= class_to:
        return

    y_left_unique, y_right_unique = sffs(current_classes, random_state,
                                         mutualInformationMatrix,
                                         classFrequencies)

    code_book[:, class_from] = 0
    code_book[y_left_unique, class_from] = 1
    code_book[y_right_unique, class_from] = -1

    left_count = y_left_unique.size

    parse_decoc(y_left_unique, code_book,
                class_from + 1, class_from + left_count,
                random_state, mutualInformationMatrix, classFrequencies)

    parse_decoc(y_right_unique, code_book,
                class_from + left_count, class_to,
                random_state, mutualInformationMatrix, classFrequencies)


def sffs(current_classes, random_state,
         mutualInformationMatrix, classFrequencies):

    y_left, y_right = random_split(current_classes, random_state)
    current_qmi = calculate_qm_information(y_left,
                                           y_right,
                                           mutualInformationMatrix,
                                           classFrequencies)

    bp_params = collections.namedtuple('BinaryPartitionParams',
                                       ['y_left', 'y_right',
                                        'qmi', ], verbose=False)

    bp_params.y_left = y_left
    bp_params.y_right = y_right
    bp_params.qmi = current_qmi

    qmi_improves = True
    while qmi_improves:
        add_class_to_binary_partition(bp_params, random_state,
                                      mutualInformationMatrix,
                                      classFrequencies)

        remove_class_to_binary_partition(bp_params, random_state,
                                         mutualInformationMatrix,
                                         classFrequencies)

        qmi_improves = (abs(current_qmi - bp_params.qmi) > 0.001)
        current_qmi = bp_params.qmi

    return bp_params.y_left, bp_params.y_right


def add_class_to_binary_partition(bp_params, random_state,
                                  mutualInformationMatrix,
                                  classFrequencies):

    labels_right = bp_params.y_right

    if len(labels_right) > 0:
        random_state.shuffle(labels_right)
        labels_right = np.delete(labels_right, -1)

    for label in labels_right:
        y_pot_left = np.concatenate((bp_params.y_left,
                                     np.array([label])), axis=0)
        y_pot_right = bp_params.y_right[bp_params.y_right != label]

        qmi = calculate_qm_information(y_pot_left,
                                       y_pot_right,
                                       mutualInformationMatrix,
                                       classFrequencies)

        if qmi > bp_params.qmi:
            bp_params.qmi = qmi
            bp_params.y_left = y_pot_left
            bp_params.y_right = y_pot_right


def remove_class_to_binary_partition(bp_params, random_state,
                                     mutualInformationMatrix,
                                     classFrequencies):

    labels_left = bp_params.y_left

    if len(labels_left) > 0:
        random_state.shuffle(labels_left)
        labels_left = np.delete(labels_left, -1)

    for label in labels_left:
        y_pot_left = bp_params.y_left[bp_params.y_left != label]
        y_pot_right = np.concatenate((bp_params.y_right,
                                      np.array([label])), axis=0)

        qmi = calculate_qm_information(y_pot_left,
                                       y_pot_right,
                                       mutualInformationMatrix,
                                       classFrequencies)

        if qmi > bp_params.qmi:
            bp_params.qmi = qmi
            bp_params.y_left = y_pot_left
            bp_params.y_right = y_pot_right


def random_split(current_classes, random_state):
    left_classes = np.array([current_classes[0]])
    right_classes = np.array([current_classes[-1]])

    current_classes = np.delete(current_classes, 0)
    current_classes = np.delete(current_classes, -1)

    if len(current_classes) > 0:
        random_assingments = np.array([random_state.randint(0, 1)
                                       for p in range(len(current_classes))])

        left_classes_random = current_classes[random_assingments == 1]
        right_classes_random = current_classes[random_assingments == 0]

        left_classes = np.append(left_classes,
                                 left_classes_random)

        right_classes = np.append(right_classes,
                                  right_classes_random)

    return left_classes, right_classes


def calculate_qm_information(y_left, y_right,
                             mutualInformationMatrix,
                             classFrequencies):

    vin = calculate_vin(y_left, y_right,
                        mutualInformationMatrix, classFrequencies)

    vall = calculate_vall(y_left, y_right,
                          mutualInformationMatrix, classFrequencies)

    vbtw = calculate_vbtw(y_left, y_right,
                          mutualInformationMatrix, classFrequencies)

    quadraticMI = vin + vall - 2*vbtw

    return quadraticMI


def count_frequencies(y, classFrequencies):
    return sum(classFrequencies[y])

def calculate_vin(y_left, y_right, mutualInformationMatrix, classFrequencies):

    vin_left = sum([mutualInformationMatrix[a, b]
                    for a in y_left
                    for b in y_left])

    vin_right = sum([mutualInformationMatrix[a, b]
                     for a in y_right
                     for b in y_right])

    J_left = count_frequencies(y_left, classFrequencies)
    J_right = count_frequencies(y_right, classFrequencies)
    total = J_left + J_right

    vin = vin_left + vin_right
    vin /= np.square(total)
    return vin


def calculate_vall(y_left, y_right,
                   mutualInformationMatrix, classFrequencies):

    y = np.concatenate((y_left, y_right), axis=1)

    vall_estimates = sum([mutualInformationMatrix[a, b]
                          for a in y
                          for b in y])

    J_left = count_frequencies(y_left, classFrequencies)
    J_right = count_frequencies(y_right, classFrequencies)
    total = J_left + J_right

    vall_classes = float(np.square(J_left) + np.square(J_right))
    vall_classes /= float(np.square(total))

    vall_all = vall_estimates*vall_classes
    vall_all /= float(np.square(total))

    return vall_all


def calculate_vbtw(y_left, y_right,
                   mutualInformationMatrix, classFrequencies):

    y = np.concatenate((y_left, y_right), axis=1)

    vbtw_left = sum([mutualInformationMatrix[a, b]
                     for a in y
                     for b in y_left])

    vbtw_right = sum([mutualInformationMatrix[a, b]
                      for a in y
                      for b in y_right])

    J_left = count_frequencies(y_left, classFrequencies)
    J_right = count_frequencies(y_right, classFrequencies)
    total = J_left + J_right

    vbtw = (J_left/float(total)) * vbtw_left
    vbtw += (J_right/float(total)) * vbtw_right

    vbtw /= float(np.square(total))
    return vbtw


def calculate_sigma(X):
    """
    Finds the sigma uses for calculating Quadratic Mutual Information

    Parameters
    ----------
    X : nump array
        Array containing input daya

    Returns
    --------
    sigma : float
        uses for calculating Quadratic Mutual Information
    """

    sigma = max([np.square(X-x).sum(axis=1).max() for x in X])
    sigma = float(sigma)*0.5
    return sigma
