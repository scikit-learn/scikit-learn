import sklearn
import multiprocessing as mp
import random
import time
import sklearn.datasets as datasets

def fetch_data(i):
    print(f"Worker {i} started", flush=True)

    time.sleep(random.random())

    data = datasets.fetch_20newsgroups(
        subset="train",
        shuffle=True,
        random_state=42
    )

    print(f"Worker {i} finished", flush=True)

    return (i, len(data.data))

if __name__ == "__main__":
    with mp.Pool(processes=8) as pool:
        results = pool.map(fetch_data, range(8))

    print(results)