"""
Small example showing recursive logging in an object hierarchy.
"""

from time import sleep
import itertools

from monologue.progress_log import HasLog

FIRST_NAMES = itertools.cycle(['Jane', 'Joe', 'Jack'])

class Employee(HasLog):

    def __init__(self, name='Joe Average', verbose=False):
        self.name = name
        self.verbose = verbose

    def work(self, chore_msg):
        log = self._get_logger()
        sleep(.2)
        log.progress('%s says "Done my chores %s"',
                msg_vars=(self.name, chore_msg))


class Boss(HasLog):

    def __init__(self, n_employees=3, verbose=False):
        self.verbose = verbose
        self.n_employees = n_employees

    def yell(self):
        log = self._get_logger()
        log.progress('Get to work!!')
        employes = [Employee(name='%s Average' % n,
                            verbose=log.clone())
                    for _, n in zip(range(self.n_employees),
                                 FIRST_NAMES)]

        for employe in employes:
            employe.work('code')


if __name__ == '__main__':
    boss = Boss(verbose=10000)
    boss.yell()

    from sklearn.progress_log import setup_logger
    import logging
    setup_logger('__main__', level=logging.DEBUG, display_name=True,
              timestamp=True)
    boss.yell()
