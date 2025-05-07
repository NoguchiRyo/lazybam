import pytest
import bampy


def test_sum_as_string():
    assert bampy.sum_as_string(1, 1) == "2"
