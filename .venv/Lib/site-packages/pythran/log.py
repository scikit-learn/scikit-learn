import logging

logger = logging.getLogger("pythran")
stream = logging.StreamHandler()

# Initialize logging
try:
    # Set a nice colored output
    from colorlog import ColoredFormatter
    formatter = ColoredFormatter(
        "%(log_color)s%(levelname)-8s%(reset)s %(blue)s%(message)s",
        log_colors={
            'DEBUG':    'cyan',
            'INFO':     'green',
            'WARNING':  'yellow',
            'ERROR':    'red',
            'CRITICAL': 'red',
        }
    )
except ImportError:
    # No color available, use default config
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    color_disabled = True
else:
    color_disabled = False

stream.setFormatter(formatter)
logger.addHandler(stream)

if color_disabled:
    logger.info("Disabling color, you really want to install colorlog.")

