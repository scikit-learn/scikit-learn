import os
import pickle
import typing as T

def run(args: T.List[str]) -> int:
    with open(args[0], "wb") as f:
        pickle.dump(dict(os.environ), f)
    return 0
