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
    code_book = np.ones(shape=(n_classes, n_classes-1), dtype=np.float)

    class_from = 0
    class_to = n_classes-1

    parse_decoc(X, y, code_book, class_from,  class_to, random_state)

    return code_book


def parse_decoc(X, y, code_book, class_from,  class_to, random_state):
    if class_from >= class_to:
        return

    bp_params = sffs(X, y, random_state)
    y_left_unique = np.unique(bp_params.y_left)
    y_right_unique = np.unique(bp_params.y_right)

    code_book[:, class_from] = 0
    code_book[y_left_unique, class_from] = 1
    code_book[y_right_unique, class_from] = -1

    left_count = y_left_unique.size

    parse_decoc(bp_params.X_left, bp_params.y_left, code_book,
                class_from + 1, class_from + left_count, random_state)

    parse_decoc(bp_params.X_right, bp_params.y_right, code_book,
                class_from + left_count, class_to, random_state)


def sffs(X, y, random_state):
    sigma = calculate_sigma(X)

    X_left, y_left, X_right, y_right = random_split(X, y, random_state)
    current_qmi = find_quadratic_mutual_information(X_left, X_right, sigma)

    bp_params = collections.namedtuple('BinaryPartitionParams',
                                       ['X_left', 'y_left',
                                        'X_right', 'y_right',
                                        'qmi', 'sigma'], verbose=False)
    bp_params.X_left = X_left
    bp_params.y_left = y_left
    bp_params.X_right = X_right
    bp_params.y_right = y_right
    bp_params.sigma = sigma
    bp_params.qmi = current_qmi

    qmi_improves = True
    while qmi_improves:
        add_class_to_binary_partition(bp_params, random_state)
        remove_class_to_binary_partition(bp_params, random_state)

        qmi_improves = (abs(current_qmi - bp_params.qmi) > 0.001)
        current_qmi = bp_params.qmi

    return bp_params


def add_class_to_binary_partition(bp_params, random_state):

    unique_labels_right = np.unique(bp_params.y_right)

    if len(unique_labels_right) > 0:
        random_state.shuffle(unique_labels_right)
        unique_labels_right = np.delete(unique_labels_right, -1)

    for label in unique_labels_right:
        X_label_samples = bp_params.X_right[bp_params.y_right == label, :]
        y_label_samples = bp_params.y_right[bp_params.y_right == label, :]

        X_pot_left = np.concatenate((bp_params.X_left,
                                     X_label_samples), axis=0)
        X_pot_right = bp_params.X_right[bp_params.y_right != label, :]

        y_pot_left = np.concatenate((bp_params.y_left,
                                     y_label_samples), axis=0)
        y_pot_right = bp_params.y_right[bp_params.y_right != label]

        qmi = find_quadratic_mutual_information(X_pot_left,
                                                X_pot_right,
                                                bp_params.sigma)

        if qmi > bp_params.qmi:
            bp_params.qmi = qmi
            bp_params.X_left = X_pot_left
            bp_params.X_right = X_pot_right
            bp_params.y_left = y_pot_left
            bp_params.y_right = y_pot_right


def remove_class_to_binary_partition(bp_params, random_state):

    unique_labels_left = np.unique(bp_params.y_left)

    if len(unique_labels_left) > 0:
        random_state.shuffle(unique_labels_left)
        unique_labels_left = np.delete(unique_labels_left, -1)

    for label in unique_labels_left:
        X_label_samples = bp_params.X_left[bp_params.y_left == label, :]
        y_label_samples = bp_params.y_left[bp_params.y_left == label, :]

        X_pot_left = bp_params.X_left[bp_params.y_left != label, :]
        X_pot_right = np.concatenate((bp_params.X_right,
                                      X_label_samples), axis=0)

        y_pot_left = bp_params.y_left[bp_params.y_left != label]
        y_pot_right = np.concatenate((bp_params.y_right,
                                      y_label_samples), axis=0)

        qmi = find_quadratic_mutual_information(X_pot_left,
                                                X_pot_right,
                                                bp_params.sigma)

        if qmi > bp_params.qmi:
            bp_params.qmi = qmi
            bp_params.X_left = X_pot_left
            bp_params.X_right = X_pot_right
            bp_params.y_left = y_pot_left
            bp_params.y_right = y_pot_right


def random_split(X, y, random_state):
    unique_labels = np.unique(y)

    left_classes = np.array(unique_labels[0])
    right_classes = np.array(unique_labels[-1])

    unique_labels = np.delete(unique_labels, 0)
    unique_labels = np.delete(unique_labels, -1)

    if len(unique_labels) > 0:
        random_assingments = np.array([random_state.randint(0, 1)
                                       for p in range(len(unique_labels))])

        left_classes_random = unique_labels[random_assingments == 1]
        right_classes_random = unique_labels[random_assingments == 0]

        left_classes = np.append(left_classes,
                                 left_classes_random)

        right_classes = np.append(right_classes,
                                  right_classes_random)

    left_samples = np.array([el in left_classes for el in y])
    right_samples = np.array([el in right_classes for el in y])

    X_left = X[left_samples]
    y_left = y[left_samples]

    X_right = X[right_samples]
    y_right = y[right_samples]

    return X_left, y_left, X_right, y_right


def find_quadratic_mutual_information(X_left, X_right, sigma):
    X = np.concatenate((X_left, X_right), axis=0)
    y_left = np.ones(shape=(1, X_left.shape[0]))
    y_right = np.ones(shape=(1, X_right.shape[0]))*(-1)

    y = np.concatenate((y_left, y_right), axis=1)
    qmi = calculate_quadratic_mutal_information(X, y, sigma)
    return qmi


def calculate_quadratic_mutal_information(X, y, sigma):
    vin = calculate_vin(X, y, sigma)
    vall = calculate_vall(X, y, sigma)
    vbtw = calculate_vbtw(X, y, sigma)

    quadraticMI = vin + vall - 2*vbtw

    return quadraticMI


def calculate_vin(X, y, sigma):
    unique_labels = np.unique(y)
    idx = y.reshape(X.shape[0])
    vin = sum([calculate_parzen_estimate(x_pl, x_pk, sigma)
               for c in unique_labels
               for x_pl in X[idx == c, :]
               for x_pk in X[idx == c, :]])

    vin /= np.square(2)
    return vin


def calculate_vall(X, y, sigma):
    unique_labels = np.unique(y)
    N = float(len(y))

    vall_estimates = sum([calculate_parzen_estimate(x_l, x_k, sigma)
                          for x_l in X
                          for x_k in X])

    vall_classes = sum([np.square(len(y[y == c])/N)
                        for c in unique_labels])
    vall_all = vall_estimates*vall_classes
    vall_all /= np.square(N)

    return vall_all


def calculate_vbtw(X, y, sigma):
    unique_labels = np.unique(y)
    N = float(len(y))

    idx = y.reshape(X.shape[0])
    vbtw = sum([len(y[y == c])/N*calculate_parzen_estimate(x_l, x_pk, sigma)
                for c in unique_labels
                for x_l in X
                for x_pk in X[idx == c, :]])

    vbtw /= np.square(N)
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


def calculate_parzen_estimate(x, y, sigma):
    if x.shape != y.shape:
        raise ValueError('X and Y vectors should have the same length!')

    estimate = np.sqrt(np.matrix(x)*np.matrix(y).T)
    return stats.norm.pdf(estimate, sigma)
