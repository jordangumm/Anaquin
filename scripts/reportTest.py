import unittest
from report import *

class TestRun(unittest.TestCase):
    def test1(self):
        self.assertTrue(len(run("whoami")) > 0)

if __name__ == '__main__':
    unittest.main()
