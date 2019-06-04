from category_encoders.hashing import HashingEncoder
from nhashing import NHashingEncoder
import time
import pandas as pd
from sklearn.datasets import load_boston
bunch = load_boston()
X = pd.DataFrame(bunch.data, columns=bunch.feature_names)
DL = []
for i in range(1000): DL.append(X)
DF = pd.concat(DL, ignore_index=True).reset_index(drop=True)
he = HashingEncoder(cols=['CHAS', 'RAD'])
start = time.time()
he.fit_transform(DF)
he_time = time.time() - start
nhe = NHashingEncoder(cols=['CHAS', 'RAD'], max_process=4)
start = time.time()
nhe.fit_transform(DF)
nhe_time = time.time() - start
print("HashingEncoder Time : ", he_time, "NHashingEncoder Time : ", nhe_time)