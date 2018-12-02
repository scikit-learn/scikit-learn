__all__ = ("getDataset",)
from io import StringIO
import pandas
from .GitHubDatasets import *
from ..utils import Bunch


defaultRepoName = "vincentarelbundock/Rdatasets"

def getRDatasetsIndex(repoName=None, refspec="master"):
    if not repoName:
        repoName = defaultRepoName
    import csv
    res = {}
    repoRawPath = GH_RAW_BASE+repoName
    index = requests.get(repoRawPath + "/" + refspec + "/datasets.csv").text
    reader = csv.reader(index.splitlines())
    headers = next(reader)
    
    hasDicMapping = {}
    hasPrefix = "has_"
    
    for i in range(len(headers)):
        headers[i] = headers[i].lower()
        if headers[i].startswith(hasPrefix):
            hasDicMapping[headers[i][len(hasPrefix):]] = headers[i]
    
    for r in reader:
        row = dict(zip(headers, r))
        pkg = row['package']
        del(row['package'])
        id = row['item']
        
        hasDic = {}
        for hasTypeName, hasTypeColumnName in hasDicMapping.items():
            hasDic[hasTypeName] = bool(row[hasTypeColumnName])
            del(row[hasTypeColumnName])
        
        row["has"] = hasDic
        
        row["shape"] = (row["rows"], row["cols"])
        del(row["rows"])
        del(row["cols"])
        
        if pkg not in res:
            res[pkg] = {}
        
        helpPath = row["doc"].split("/")
        helpPath.insert(-1, "rst")
        helpPath[-1] = ".".join(helpPath[-1].split(".")[:-1]+["rst"])
        row["doc"] = "/".join(helpPath)
        res[pkg][id] = row
        
    return res

index = None


def getDataset(category, name):
    """Returns a dataset from """ + defaultRepoName
    
    global index
    if not index:
        index = getRDatasetsIndex()
    
    metadata = index[category][name]
    csv = requests.get(metadata["csv"]).text
    desc = requests.get(metadata["doc"]).text
    
    with StringIO(csv) as csvStr:
        data = pandas.read_csv(csvStr)
    
    additionalInfo = parseInfoFromDocs(desc)
    
    return Bunch(X=data, y=None, target=None, feature_names=data.columns, DESCR=metadata["title"], details=desc, url=metadata["csv"], **additionalInfo)

def parseInfoFromDocs(desc):
    import docutils.parsers.rst
    import docutils.utils
    import docutils.frontend
    s = docutils.frontend.OptionParser(components=(docutils.parsers.rst.Parser,)).get_default_values()
    p = docutils.parsers.rst.Parser()
    d = docutils.utils.new_document(None, s)
    p.parse(desc, d)
    parsedDocs = d.asdom()
    sections = {}
    
    for s in parsedDocs.getElementsByTagName("section"):
        tit = s.getElementsByTagName("title")
        nm = tit[0].childNodes[0].data
        sections[nm] = s
        
    parsedParams = getMetaInfoFromDocs(dd)
    return {"columnsDescriptions": parseColumnsInfo(sections["Format"])}

def parseColumnsInfo(formatSect):
    parsedParams = {}
    ft = formatSect.getElementsByTagName("table")[0]
    for row in ft.getElementsByTagName("row"):
        if len(row.childNodes) == 2:
            resNodes = []
            for s in row.childNodes:
                ps = s.getElementsByTagName("paragraph")
                if len(ps) == 1:
                    p = ps[0]
                    ts = p.childNodes
                    assert(len(ts) == 1)
                    resNodes.append(ts[0].data)
            if len(resNodes) == 2:
                k, v = resNodes
                if k[-1] == ":":
                    k = k[:-1]
                parsedParams[k] = v
    
    return parsedParams

