import logging


logger = logging.getLogger()

logging.basicConfig(
    format="# %(asctime)s %(levelname)s %(name)s %(filename)s:%(lineno)s -- %(message)s\n",
    level=logging.INFO,
)
