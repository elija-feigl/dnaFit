import logging
from pathlib import Path


def get_resource(resources: str) -> Path:
    return Path(__file__).parent / "resources" / resources


def _init_logging():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s | [%(name)s] %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)


_init_logging()
