"""
tests for logger.py prior to coding

I'll fix docstring formatting, I promise.

We are testing the global Logger level, but this does not harm the
possibility of adding a log handler to a logger, with a specific log level.
However, most progress_* functions do bypass the traditional logging mechanism.

Do you want level indicators in messages?

Creation of loggers
===================
>>> from scikits.learn.utils import get_logger, PROGRESS
... # PROGRESS is a custom level (DEBUG + 5)
>>> from logging import DEBUG, INFO, WARNING
>>> verbose_logger = get_logger("scikits.learn.utils.testlogger_verbose",
... verbosity_offset=-10)
>>> standard_logger = get_logger("scikits.learn.utils.testlogger_standard")
>>> laconic_logger = get_logger("scikits.learn.utils.testlogger_laconic",
... verbosity_offset=+10)

drop in replacement for print: print_msg
========================================

Always print if verbosity not specified
-------------------------------------
>>> verbose_logger.print_msg("Message must be displayed")
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.print_msg("Message must be displayed")
[scikits.learn.utils.testlogger_standard] Message must be displayed
>>> laconic_logger.print_msg("Message must be displayed")
[scikits.learn.utils.testlogger_laconic] Message must be displayed

Print according to verbosity
-------------------------------------
- verbosity=False => DEBUG
- verbosity=True => always print
- verbosity=int => use int as verb level
Maybe this makes the case of “verbosity in (0, 1)” counter intuitive?

False
>>> verbose_logger.print_msg("Message must be displayed", verbosity=False)
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.print_msg("Message must'nt be displayed", verbosity=False)
>>> laconic_logger.print_msg("Message must'nt be displayed", verbosity=False)

True
>>> verbose_logger.print_msg("Message must be displayed", verbosity=True)
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.print_msg("Message must be displayed", verbosity=True)
[scikits.learn.utils.testlogger_standard] Message must be displayed
>>> laconic_logger.print_msg("Message must be displayed", verbosity=True)
[scikits.learn.utils.testlogger_laconic] Message must be displayed

0 (same for 1)
>>> verbose_logger.print_msg("Message mustn' be displayed", verbosity=0)
>>> standard_logger.print_msg("Message must'nt be displayed", verbosity=0)
>>> laconic_logger.print_msg("Message must'nt be displayed", verbosity=0)

DEBUG (10)
>>> verbose_logger.print_msg("Message must be displayed", verbosity=DEBUG)
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.print_msg("Message must'nt be displayed", verbosity=DEBUG)
>>> laconic_logger.print_msg("Message must'nt be displayed", verbosity=DEBUG)

INFO (20)
>>> verbose_logger.print_msg("Message must be displayed", verbosity=INFO)
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.print_msg("Message must be displayed", verbosity=INFO)
[scikits.learn.utils.testlogger_standard] Message must be displayed
>>> laconic_logger.print_msg("Message must'nt be displayed", verbosity=INFO)

WARN (30)
>>> verbose_logger.print_msg("Message must be displayed", verbosity=WARNING)
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.print_msg("Message must be displayed", verbosity=WARNING)
[scikits.learn.utils.testlogger_standard] Message must be displayed
>>> laconic_logger.print_msg("Message must be displayed", verbosity=WARNING)
[scikits.learn.utils.testlogger_laconic] Message must be displayed

Play with verbosity offset
--------------------------
For absolute verbosity, use setLevel (see below).
>>> laconic_logger.print_msg("Message must'nt be displayed", verbosity=False)
>>> laconic_logger.set_offset(-10) # relative to general level which is INFO
>>> laconic_logger.print_msg("Message must be displayed", verbosity=False)
[scikits.learn.utils.testlogger_laconic] Message must be displayed
>>> laconic_logger.add_to_offset(20)
... #add_to_offset will get us back to the initial +10 value
>>> laconic_logger.print_msg("Message must'nt be displayed", verbosity=False)

traditional logging functions still work
========================================
>>> verbose_logger.debug("Message must be displayed")
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.debug("Message mustn't be displayed")
>>> laconic_logger.debug("Message mustn't be displayed")

>>> verbose_logger.info("Message must be displayed")
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.info("Message must be displayed")
[scikits.learn.utils.testlogger_standard] Message must be displayed
>>> laconic_logger.info("Message mustn't be displayed")

>>> verbose_logger.warning("Message must be displayed")
[scikits.learn.utils.testlogger_verbose] Message must be displayed
>>> standard_logger.warning("Message must be displayed")
[scikits.learn.utils.testlogger_standard] Message must be displayed
>>> laconic_logger.warning("Message must be displayed")
[scikits.learn.utils.testlogger_laconic] Message must be displayed

Play with verbosity (absolute)
------------------------------
>>> laconic_logger.setLevel(DEBUG)
>>> laconic_logger.debug("Message must be displayed")
[scikits.learn.utils.testlogger_laconic] Message must be displayed
>>> laconic_logger.setLevel(PROGRESS)
>>> laconic_logger.setLevel(WARNING)


progress
========

Initialization
>>> for logger in (verbose_logger, standard_logger, laconic_logger):
...    logger.progress_every(1000)
>>> for logger in (verbose_logger, standard_logger, laconic_logger):
...    logger.dot_every(10)

progress_dot: Same semantic as print_msg
----------------------------------------

no verbosity: always output
>>> verbose_logger.progress_dot()
.
>>> standard_logger.progress_dot()
.
>>> laconic_logger.progress_dot()
.

verbosity=False -> DEBUG
>>> verbose_logger.progress_dot(verbosity=False)
.
>>> standard_logger.progress_dot(verbosity=False)
>>> laconic_logger.progress_dot(verbosity=False)

verbosity=True -> always print
>>> verbose_logger.progress_dot(verbosity=True)
.
>>> standard_logger.progress_dot(verbosity=True)
.
>>> laconic_logger.progress_dot(verbosity=True)
.

verbosity=PROGRESS for instance
>>> verbose_logger.progress_dot(verbosity=PROGRESS)
.
>>> standard_logger.progress_dot(verbosity=PROGRESS)
>>> laconic_logger.progress_dot(verbosity=PROGRESS)

progress_step
-------------
>>> for count in range(9):
...     verbose_logger.progress_step()
>>> verbose_logger.progress_step()
.
>>> for count in range(90):
...     verbose_logger.progress_step()
.........
>>> for count in range(2000):
>>>     verbose_logger.progress_dot() # doctest: +ELLIPSIS
...
[scikits.learn.utils.testlogger_verbose] Iteration 1000 done
...
[scikits.learn.utils.testlogger_verbose] Iteration 2000 done

>>> for count in range(2000): # doctest: -ELLIPSIS
>>>     standard_logger.progress_dot()
[scikits.learn.utils.testlogger_standard] Iteration 1000 done
[scikits.learn.utils.testlogger_standard] Iteration 2000 done

>>> for count in range(2000):
>>>     laconic_logger.progress_dot()

FROM HERE: questionnable ideas
>>> standard_logger.progress_complete()
[scikits.learn.utils.testlogger_standard] Successfully completed 2000 iterations
>>> standard_logger.percent_print_every(10) # every 10 percent, requires a target
>>> standard_logger.percent_target(1000) # the scale. rename function?
>>> for count in range(2000):
>>>     standard_logger.progress_step()
[scikits.learn.utils.testlogger_standard] 10%
[scikits.learn.utils.testlogger_standard] 20%
[scikits.learn.utils.testlogger_standard] 30%
[scikits.learn.utils.testlogger_standard] 40%
[scikits.learn.utils.testlogger_standard] 50%
[scikits.learn.utils.testlogger_standard] 60%
[scikits.learn.utils.testlogger_standard] 70%
[scikits.learn.utils.testlogger_standard] 80%
[scikits.learn.utils.testlogger_standard] 90%
[scikits.learn.utils.testlogger_standard] 100%
[scikits.learn.utils.testlogger_standard] 110%
[scikits.learn.utils.testlogger_standard] 120%
[scikits.learn.utils.testlogger_standard] 130%
[scikits.learn.utils.testlogger_standard] 140%
[scikits.learn.utils.testlogger_standard] 150%
[scikits.learn.utils.testlogger_standard] 160%
[scikits.learn.utils.testlogger_standard] 170%
[scikits.learn.utils.testlogger_standard] 180%
[scikits.learn.utils.testlogger_standard] 190%
[scikits.learn.utils.testlogger_standard] 200%


"""

