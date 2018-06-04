import tempfile, os

from traitlets.config.loader import Config
import nose.tools as nt

ip = get_ipython()
ip.magic('load_ext storemagic')

def test_store_restore():
    ip.user_ns['foo'] = 78
    ip.magic('alias bar echo "hello"')
    tmpd = tempfile.mkdtemp()
    ip.magic('cd ' + tmpd)
    ip.magic('store foo')
    ip.magic('store bar')
    
    # Check storing
    nt.assert_equal(ip.db['autorestore/foo'], 78)
    nt.assert_in('bar', ip.db['stored_aliases'])
    
    # Remove those items
    ip.user_ns.pop('foo', None)
    ip.alias_manager.undefine_alias('bar')
    ip.magic('cd -')
    ip.user_ns['_dh'][:] = []
    
    # Check restoring
    ip.magic('store -r')
    nt.assert_equal(ip.user_ns['foo'], 78)
    assert ip.alias_manager.is_alias('bar')
    nt.assert_in(os.path.realpath(tmpd), ip.user_ns['_dh'])
    
    os.rmdir(tmpd)

def test_autorestore():
    ip.user_ns['foo'] = 95
    ip.magic('store foo')
    del ip.user_ns['foo']
    c = Config()
    c.StoreMagics.autorestore = False
    orig_config = ip.config
    try:
        ip.config = c
        ip.extension_manager.reload_extension('storemagic')
        nt.assert_not_in('foo', ip.user_ns)
        c.StoreMagics.autorestore = True
        ip.extension_manager.reload_extension('storemagic')
        nt.assert_equal(ip.user_ns['foo'], 95)
    finally:
        ip.config = orig_config
