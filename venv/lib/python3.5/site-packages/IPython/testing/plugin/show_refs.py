"""Simple script to show reference holding behavior.

This is used by a companion test case.
"""

import gc

class C(object):
   def __del__(self):
      pass
      #print 'deleting object...'  # dbg

if __name__ == '__main__':
   c = C()

   c_refs = gc.get_referrers(c)
   ref_ids = list(map(id,c_refs))

   print('c referrers:',list(map(type,c_refs)))
