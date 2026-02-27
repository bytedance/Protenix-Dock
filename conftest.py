import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "e2e: end-to-end tests requiring AmberTools")
    config.addinivalue_line("markers", "slow: slow tests")
