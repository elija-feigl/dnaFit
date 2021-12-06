import pytest
import numpy as np

from dnaFit.core.utils import _norm

assert np.linalg.norm(_norm(np.random.rand(3))) == pytest.approx(1.)
