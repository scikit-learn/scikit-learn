import json
import warnings
import numpy as np
import numbers

try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen

from scipy.io.arff import loadarff

_SEARCH_NAME = "https://openml.org/api/v1/json/data/list/data_name/{}/limit/1"

jsons = "https://openml.org/api/v1/json/data/list/data_name/{}"
data_dl = "https://www.openml.org/data/download/{}"


def _get_data_info_by_name(name):
    url_path = urlopen(_SEARCH_NAME.format(name))                                                                               
    json_data = json.load(url_path)
    return json_data['data']['dataset'][0]


def fetch_openml(name_or_id=None, version=1, json_loc=jsons,
                 data_loc=data_dl):
    """Fetch dataset from openml by name or dataset id.

    Parameters
    ----------

    Returns
    -------
    """
    if isinstance(name_or_id, numbers.Integral):
        if version != 1:
            raise ValueError(
                "Dataset id={} and version={} passed, but you can only "
                "specify a numeric id or a version, not both.".format(
                    name_or_id, version))
        data_id = name_or_id
    elif isinstance(name_or_id, str):
        name = name_or_id

    else:
        raise TypeError(
            "Invalid name_or_id {}, should be string or integer.".format(
                name_or_id))

    json_dl = urlretrieve(json_loc.format(name))[0]
    # get the json file
    with open(json_dl, 'r') as tmp:
        json_data = json.load(tmp)['data']['dataset']
    vers = [(idx, val) for idx, item in enumerate(json_data)
            for key, val in item.items() if key == "version"]
    # tell user there are more versions if they dont specify number
    if len(vers) > 1 and name_vers is None:
        msg = ("dataset: {} has versions {}, "
               "default is {}").format(name,
                                       [i[1] for i in vers],
                                       min([i[1] for i in vers]))
        warnings.warn(msg)
    # check if the version specified (if it is) is in the ones gotten
    use = 1 if name_vers is None else name_vers
    for v in vers:
        if v[1] == use:
            to_get = json_data[v[0]]['file_id']
    # download data
    data_tmp = urlretrieve(data_loc.format(to_get))[0]
    # load the data
    data = loadarff(data_tmp)
    data_fmt = np.zeros((data[0].shape[0], len(data[0][0])), dtype=object)
    # scipy returns a tuple so try to put it in the right format
    for idx, row in enumerate(data[0]):
        data_fmt[idx, :] = [val for val in row]
    return data_fmt
