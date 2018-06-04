from __future__ import print_function, division, absolute_import

import pytest
from toolz import partial

from dask import compute, get
from dask.utils import filetexts
from dask.bytes import compression
from dask.bag.text import read_text

compute = partial(compute, get=get)


files = {'.test.accounts.1.json':  ('{"amount": 100, "name": "Alice"}\n'
                                    '{"amount": 200, "name": "Bob"}\n'
                                    '{"amount": 300, "name": "Charlie"}\n'
                                    '{"amount": 400, "name": "Dennis"}\n'),
         '.test.accounts.2.json':  ('{"amount": 500, "name": "Alice"}\n'
                                    '{"amount": 600, "name": "Bob"}\n'
                                    '{"amount": 700, "name": "Charlie"}\n'
                                    '{"amount": 800, "name": "Dennis"}\n')}


expected = ''.join([files[v] for v in sorted(files)])

fmt_bs = ([(fmt, None) for fmt in compression.files] +
          [(fmt, 10) for fmt in compression.seekable_files] +
          [(fmt, None) for fmt in compression.seekable_files])
encodings = ['ascii', 'utf-8'] # + ['utf-16', 'utf-16-le', 'utf-16-be']
fmt_bs_enc = [(fmt, bs, encoding) for fmt, bs in fmt_bs
              for encoding in encodings]


@pytest.mark.parametrize('fmt,bs,encoding', fmt_bs_enc)
def test_read_text(fmt, bs, encoding):
    compress = compression.compress[fmt]
    files2 = dict((k, compress(v.encode(encoding))) for k, v in files.items())
    with filetexts(files2, mode='b'):
        b = read_text('.test.accounts.*.json', compression=fmt, blocksize=bs,
                      encoding=encoding)
        L, = compute(b)
        assert ''.join(L) == expected

        b = read_text(sorted(files), compression=fmt, blocksize=bs,
                      encoding=encoding)
        L, = compute(b)
        assert ''.join(L) == expected

        blocks = read_text('.test.accounts.*.json', compression=fmt, blocksize=bs,
                           encoding=encoding, collection=False)
        L = compute(*blocks)
        assert ''.join(line for block in L for line in block) == expected


def test_errors():
    with filetexts({'.test.foo': b'Jos\xe9\nAlice'}, mode='b'):
        with pytest.raises(UnicodeDecodeError):
            read_text('.test.foo', encoding='ascii').compute()

        result = read_text('.test.foo', encoding='ascii', errors='ignore')
        result = result.compute(get=get)
        assert result == ['Jos\n', 'Alice']
